#include "dada_affinity.h"
#include "dada_client.h"
#include "dada_hdu.h"
#include "dada_def.h"

#include "ascii_header.h"

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include <fcntl.h>
#include <assert.h>
#include <math.h>

#include <sys/types.h>
#include <sys/stat.h>

typedef struct 
{
  dada_hdu_t *  hdu;
  key_t         key;
  uint64_t      block_size;
  uint64_t      bytes_written;
  uint64_t      sample_offset;
  unsigned      block_open;
  char *        curr_block;
} craft_dbtransposedb_hdu_t;

typedef struct {

  craft_dbtransposedb_hdu_t output;

  // number of bytes read
  uint64_t bytes_in;

  // number of bytes written
  uint64_t bytes_out;

  // verbose output
  int verbose;

  unsigned int nsig;
  unsigned int nchan;
  unsigned int npol; 
  unsigned int nbit;

  unsigned quit;

  char order[6];

  unsigned int reblock_factor;
  unsigned int pol_factor;
  unsigned int bit_factor;

  uint64_t nsamps_integrated;
  uint64_t nsamps_to_integrate;

  float * offsets;
  float * scales;
  double * sums;
  double * sums_sq;

  float digi_sigma;
  float digi_mean;
  float digi_scale;
  float digi_min;
  float digi_max;

} craft_dbtransposedb_t;

#define DADA_DBSUMDB_INIT { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }

int quit_threads = 0;

int64_t dbtransposedb_write_block_TFBP_to_BTF (dada_client_t *, void *, uint64_t, uint64_t);
void compute_scales_offsets_TFBP (craft_dbtransposedb_t *, float *, const uint64_t);

void usage()
{
  fprintf (stdout,
           "craft_dbtransposedb [options] in_key out_key\n"
           " -b core   bind processing to specified CPU core [default no binding]\n"
           " -n nsamps use nsamps in each channel and signal to determine scale factors [default 2048]\n"
           " -o nbit   requantise to nbit output: 8 or 16 [default 16]\n"
           " -r factor reblock samples by factor [default 64]\n"
           " -s        1 transfer, then exit\n"
           " -v        verbose mode\n"
           " in_key    DADA key for input data block\n"
           " out_key   DADA key for output data block\n");
}



/*! Function that opens the data transfer target */
int dbtransposedb_open (dada_client_t* client)
{
  // the craft_dbtransposedb specific data
  craft_dbtransposedb_t * ctx = (craft_dbtransposedb_t *) client->context;

  // status and error logging facilty
  multilog_t* log = client->log;

  if (ctx->verbose)
    multilog (log, LOG_INFO, "dbtransposedb_open()\n");

  char output_order[4];

  // header to copy from in to out
  char * header = 0;

  if (ctx->verbose)
    multilog (log, LOG_INFO, "open: HDU (key=%x) lock_write on HDU\n", ctx->output.key);

  if (dada_hdu_lock_write (ctx->output.hdu) < 0)
  {
    multilog (log, LOG_ERR, "cannot lock write DADA HDU (key=%x)\n", ctx->output.key);
    return -1;
  }

  // get the transfer size (if it is set)
  int64_t transfer_size = 0;
  ascii_header_get (client->header, "TRANSFER_SIZE", "%"PRIi64, &transfer_size);

  int nant;
  int nbeam;
  // get the number of antenna
  if (ascii_header_get (client->header, "NANT", "%u", &nant) != 1)
  {
    nant = 1;
  }

  if (ascii_header_get (client->header, "NBEAM", "%u", &nbeam) != 1)
  {
    nbeam = 1;
  }

  ctx->nsig = nant > nbeam ? nant : nbeam;
  if (ascii_header_get (client->header, "NBIT", "%u", &(ctx->nbit)) != 1)
  {
    multilog (log, LOG_ERR, "open: header with no NBIT\n");
    return -1;
  }

  if (ascii_header_get (client->header, "NPOL", "%u", &(ctx->npol)) != 1)
  {           
    multilog (log, LOG_ERR, "open: header with no NPOL\n");
    return -1;                
  }                             
  if (ctx->npol != 2)
  {
    multilog (log, LOG_ERR, "open: expected NPOL==2, but found %u\n", ctx->npol);
    return -1;
  }

  if (ascii_header_get (client->header, "NCHAN", "%u", &(ctx->nchan)) != 1)
  {           
    multilog (log, LOG_ERR, "open: header with no NCHAN\n");
    return -1;                
  }
 
  if (ascii_header_get (client->header, "DORDER", "%s", &(ctx->order)) != 1)
  {
    multilog (log, LOG_ERR, "open: header with no ORDER\n");
    return -1;
  }
  else
  {
    // for summing of data blocks we always want TF output mode
    multilog (log, LOG_INFO, "open: DORDER=%s\n", ctx->order);
    if (strcmp(ctx->order, "TFBP") != 0)
    {
      multilog (log, LOG_ERR, "open: input ORDER=%s is not supported\n", ctx->order);
      return -1;
    }
  }

  // Craft has TSAMP in seconds, convert to micro seconds
  double tsamp;
  if (ascii_header_get (client->header, "TSAMP", "%lf", &tsamp) != 1)
  {
    multilog (log, LOG_ERR, "open: header with no TSAMP\n");
    return -1;
  }

  // CRAFT has this at the top of the band, change to centre frequency
  double freq;
  if (ascii_header_get (client->header, "FREQ", "%lf", &freq) != 1)
  {
    multilog (log, LOG_ERR, "open: header with no FREQ\n");
    return -1;
  }

  // CRAFT uses channel bandwidth here, change to total bandwidth
  double bw;
  if (ascii_header_get (client->header, "BW", "%lf", &bw) != 1)
  {
    multilog (log, LOG_ERR, "open: header with no BW\n");
    return -1;
  }

  // CRAFT uses BLOCK_SIZE in place of RESOLUTION
  uint64_t resolution;
  if (ascii_header_get (client->header, "BLOCK_SIZE", "%lu", &resolution) != 1)
  {
    multilog (log, LOG_ERR, "open: header with no RESOLUTION\n");
    return -1;
  }
  multilog (log, LOG_INFO, "BLOCK_SIZE or RESOLUTION=%lu\n", resolution);

  char tmp[32];
  if (ascii_header_get (client->header, "UTC_START", "%s", tmp) == 1)
  {
    multilog (log, LOG_INFO, "open: UTC_START=%s\n", tmp);
  }
  else
  {
    multilog (log, LOG_INFO, "open: UTC_START=UNKNOWN\n");
  }

  // get the header from the input data block
  uint64_t header_size = ipcbuf_get_bufsz (client->header_block);
  multilog (log, LOG_INFO, "open: input header_size=%lu\n", header_size);
  multilog (log, LOG_INFO, "open: output header_size=%lu\n", ipcbuf_get_bufsz (ctx->output.hdu->header_block));

  // setup header for output HDU
  if (ctx->verbose)
    multilog (log, LOG_INFO, "open: writing HDU %x\n",  ctx->output.key);

  if (ctx->verbose)
    multilog (log, LOG_INFO, "open: enabling HDU %x\n", ctx->output.key);
  assert( header_size == ipcbuf_get_bufsz (ctx->output.hdu->header_block) );

  header = ipcbuf_get_next_write (ctx->output.hdu->header_block);
  if (!header) 
  {
    multilog (log, LOG_ERR, "open: could not get next header block\n");
    return -1;
  }

  unsigned nchansig = ctx->nchan * ctx->nsig;
  ctx->offsets = (float *) malloc (nchansig * sizeof(float));
  ctx->scales = (float *) malloc (nchansig * sizeof(float));
  ctx->sums = (double *) malloc (nchansig * sizeof(double));
  ctx->sums_sq = (double *) malloc (nchansig * sizeof(double));

  // ensure these are zeroed at the start of the observation
  bzero (ctx->sums, nchansig * sizeof(double));
  bzero (ctx->sums_sq, nchansig * sizeof(double));

  ctx->nsamps_integrated = 0;

  // copy the header from the in to the out
  memcpy ( header, client->header, header_size );

  unsigned ndim = 1;
  if (ascii_header_set (header, "NDIM", "%u", ndim) < 0)
  {
    multilog (log, LOG_ERR, "open: failed to write ndim=%u to header\n", ndim);
    return -1;
  }

  unsigned new_nbit = ctx->nbit / ctx->bit_factor;
  if (ascii_header_set (header, "NBIT", "%u", new_nbit) < 0)
  {
    multilog (log, LOG_ERR, "open: failed to write NBIT=%u to header\n", new_nbit);
    return -1;
  }

  double new_bw = bw * ctx->nchan;
  if (ascii_header_set (header, "BW", "%lf", new_bw) < 0)
  {
    multilog (log, LOG_ERR, "open: failed to write BW=%lf to header\n", new_bw);
    return -1;
  }

  // lower sideband
  double new_freq = freq + (new_bw / 2);
  if (ascii_header_set (header, "FREQ", "%lf", new_freq) < 0)
  {
    multilog (log, LOG_ERR, "open: failed to write FREQ=%lf to header\n", new_freq);
    return -1;
  }

  double new_tsamp = tsamp * 1000000;
  if (ascii_header_set (header, "TSAMP", "%lf", new_tsamp) < 0)
  {
    multilog (log, LOG_ERR, "open: failed to write TSAMP=%lf to header\n", new_tsamp);
    return -1;
  }

    multilog (log, LOG_INFO, "ctx->reblock_factor=%u ctx->bit_factor=%u ctx->pol_factor=%u\n", ctx->reblock_factor, ctx->bit_factor, ctx->pol_factor);
  unsigned out_factor = ctx->reblock_factor / (ctx->bit_factor * ctx->pol_factor);
  uint64_t new_resolution = resolution * out_factor;
  multilog (log, LOG_INFO, "resolution=%lu, out_factor=%u new resolution=%lu\n", resolution, out_factor, new_resolution);
  if (ascii_header_set (header, "RESOLUTION", "%lu", new_resolution) < 0)
  {
    multilog (log, LOG_ERR, "open: failed to write RESOLUTION=%lu to header\n", new_resolution);
    return -1;
  }

  unsigned new_npol = 1;
  if (ascii_header_set (header, "NPOL", "%u", new_npol) < 0)
  {
    multilog (log, LOG_ERR, "open: failed to write NPOL=%u to header\n", new_npol);
    return -1;
  }

  sprintf (output_order, "%s", "STF");
  if (ascii_header_set (header, "ORDER", "%s", output_order) < 0)
  {
    multilog (log, LOG_ERR, "open: failed to write ORDER=%s to header\n", output_order);
    return -1;
  }

  // mark the outgoing header as filled
  if (ipcbuf_mark_filled (ctx->output.hdu->header_block, header_size) < 0)  {
    multilog (log, LOG_ERR, "Could not mark filled Header Block\n");
    return -1;
  }
  if (ctx->verbose) 
    multilog (log, LOG_INFO, "open: HDU (key=%x) opened for writing\n", ctx->output.key);

  multilog (log, LOG_INFO, "open: transfer_size=%ld\n", transfer_size);
  client->transfer_bytes = transfer_size; 
  client->optimal_bytes = 64*1024*1024;

  ctx->bytes_in = 0;
  ctx->bytes_out = 0;
  client->header_transfer = 0;

  if (ctx->verbose) 
    multilog (log, LOG_INFO, "open: return 0\n");

  return 0;
}

int dbtransposedb_close (dada_client_t* client, uint64_t bytes_written)
{
  craft_dbtransposedb_t* ctx = (craft_dbtransposedb_t*) client->context;
  
  multilog_t* log = client->log;

  if (ctx->verbose)
    multilog (log, LOG_INFO, "close: bytes_in=%"PRIu64", bytes_out=%"PRIu64"\n",
                    ctx->bytes_in, ctx->bytes_out );

  // close the block if it is open
  if (ctx->output.block_open)
  {
    if (ctx->verbose)
      multilog (log, LOG_INFO, "close: ipcio_close_block_write bytes_written=%"PRIu64"\n");
    if (ipcio_close_block_write (ctx->output.hdu->data_block, ctx->output.bytes_written) < 0)
    {
      multilog (log, LOG_ERR, "dbtransposedb_close: ipcio_close_block_write failed\n");
      return -1;
    }
    ctx->output.block_open = 0;
    ctx->output.bytes_written = 0;
  }

  // unlock write on the datablock (end the transfer)
  if (ctx->verbose)
    multilog (log, LOG_INFO, "close: dada_hdu_unlock_write\n");

  if (dada_hdu_unlock_write (ctx->output.hdu) < 0)
  {
    multilog (log, LOG_ERR, "dbtransposedb_close: cannot unlock DADA HDU (key=%x)\n", ctx->output.key);
    return -1;
  }

  if (ctx->scales)
    free (ctx->scales);
  ctx->scales = 0;
  if (ctx->offsets)
    free (ctx->offsets);
  ctx->offsets = 0;
  if (ctx->sums)
    free (ctx->sums);
  ctx->sums = 0;
  if (ctx->sums_sq)
    free (ctx->sums_sq);
  ctx->sums_sq = 0;

  return 0;
}

/*! Pointer to the function that transfers data to/from the target */
int64_t dbtransposedb_write (dada_client_t* client, void* data, uint64_t data_size)
{
  craft_dbtransposedb_t* ctx = (craft_dbtransposedb_t*) client->context;

  multilog_t * log = client->log;

  if (ctx->verbose)
    multilog (log, LOG_INFO, "write: to_write=%"PRIu64"\n", data_size);

  // write dat to all data blocks
  ipcio_write (ctx->output.hdu->data_block, data, data_size);

  ctx->bytes_in += data_size;
  ctx->bytes_out += data_size;

  if (ctx->verbose)
    multilog (log, LOG_INFO, "write: read %"PRIu64", wrote %"PRIu64" bytes\n", data_size, data_size);
 
  return data_size;
}

int64_t dbtransposedb_write_block_TFBP_to_BTF (dada_client_t* client, void* in_data, uint64_t data_size, uint64_t block_id)
{
  craft_dbtransposedb_t* ctx = (craft_dbtransposedb_t*) client->context;

  multilog_t * log = client->log;

  if (ctx->verbose)
    multilog (log, LOG_INFO, "write_block_TFBP_to_BTF: data_size=%"PRIu64", block_id=%"PRIu64"\n",
              data_size, block_id);

  float * in = (float *) in_data;
  uint16_t * out16;
  uint8_t * out8;
 
  const uint64_t nsamp_in  = data_size / (ctx->nsig * ctx->nchan * ctx->npol * ctx->nbit / 8); 
  const uint64_t nsamp_out = nsamp_in * ctx->reblock_factor;

  uint64_t out_block_id;
  unsigned isig, isamp, ichan;

  if (ctx->verbose > 1)
    multilog (log, LOG_INFO, "write_block_TFBP_to_BTF: nsamp_in=%lu\n", nsamp_in);

  if (!ctx->output.block_open)
  {
    if (ctx->verbose > 1)
      multilog (log, LOG_INFO, "write_block_TFBP_to_BTF [%x] ipcio_open_block_write()\n", ctx->output.key);
    ctx->output.curr_block = ipcio_open_block_write(ctx->output.hdu->data_block, &out_block_id);
    if (!ctx->output.curr_block)
    {
      multilog (log, LOG_ERR, "write_block_TFBP_to_BTF [%x] ipcio_open_block_write failed %s\n", ctx->output.key, strerror(errno));
      return -1;
    }
    ctx->output.block_open = 1;
    ctx->output.bytes_written = 0;
    ctx->output.sample_offset = 0;
  }

  if (ctx->bit_factor == 2)
    out16 = (uint16_t *) ctx->output.curr_block;
  else if (ctx->bit_factor == 4)
    out8 = (uint8_t *) ctx->output.curr_block;

  // data strides for input and output blocks
  const uint64_t in_sig_stride = ctx->npol;
  const uint64_t in_chan_stride = ctx->nsig * in_sig_stride;
  const uint64_t in_samp_stride = ctx->nchan * in_chan_stride;
  const uint64_t out_samp_stride = ctx->nchan;
  const uint64_t out_sig_stride = nsamp_out * ctx->nchan;

  if (ctx->verbose > 1)
  {
    multilog (log, LOG_INFO, "write_block_TFBP_to_BTF: out_samp_stride=%lu\n", out_samp_stride);
    multilog (log, LOG_INFO, "write_block_TFBP_to_BTF: out_sig_stride=%lu\n", out_sig_stride);
    multilog (log, LOG_INFO, "write_block_TFBP_to_BTF: out_sample_offset=%lu\n", ctx->output.sample_offset);
  }

  // compute scales and offsets
  if (ctx->nsamps_integrated < ctx->nsamps_to_integrate)
  {
    if (ctx->verbose)
      multilog (log, LOG_INFO, "write_block_TFBP_to_BTF: computing scales and offsets from %lu to %lu\n",
                ctx->nsamps_integrated, ctx->nsamps_integrated + nsamp_in);
    compute_scales_offsets_TFBP (ctx, in, nsamp_in);
  }

#define NEW
#ifdef NEW 

  // loop and channel and signal first
  for (ichan=0; ichan<ctx->nchan; ichan++)
  {
    for (isig=0; isig<ctx->nsig; isig++)
    {
      // determine the required scale and offset
      const unsigned ichansig = ichan * ctx->nsig + isig;
      const float scale = ctx->scales[ichansig];
      const float mean = ctx->offsets[ichansig];
      const float min = ctx->digi_min;
      const float max = ctx->digi_max;
      const float digi_mean = ctx->digi_mean + 0.5;

      // compute the output offset sample for this channel and signal
      //uint64_t idx = (isig * in_sig_stride) + (ichan * in_chan_stride);
      uint64_t idx = ichansig * in_sig_stride;
      uint64_t odx = (isig * out_sig_stride) + (ctx->output.sample_offset * out_samp_stride) + ichan;

      for (isamp=0; isamp<nsamp_in; isamp++)
      {
        // sum the polarisations
        float pscr = in[idx] + in[idx+1];

        // compute the rescaled result
        float result = (pscr - mean) * scale + digi_mean;
        if (result < min)
          result = min;
        if (result > max)
          result = max;

        // write the output
        if (ctx->bit_factor == 2)
          out16[odx] = (uint16_t) result;
        else 
          out8[odx] = (uint8_t) result;

        idx += in_samp_stride;
        odx += out_samp_stride;
      }
    }
  }
  
#else
  // index on input ordering
  uint64_t idx = 0;
  unsigned ipol;
  for (isamp=0; isamp<nsamp_in; isamp++)
  {
    for (ichan=0; ichan<ctx->nchan; ichan++)
    {
      for (isig=0; isig<ctx->nsig; isig++)
      {
        float pscr = 0;
        for (ipol=0; ipol<ctx->npol; ipol++)
        {
          pscr += in[idx];
          idx++;
        }

        unsigned i = ichan * ctx->nsig + isig;
        const float combined_scale = ctx->scales[i] * ctx->digi_scale;
        const float mean = ctx->offsets[i];

        int result = (pscr - mean) * combined_scale + ctx->digi_mean + 0.5;

        if (result < ctx->digi_min)
          result = ctx->digi_min;
        if (result > ctx->digi_max)
          result = ctx->digi_max;
        
        const uint64_t odx = (isig * out_sig_stride) + ((isamp + ctx->output.sample_offset) * out_samp_stride) + ichan;

        if (ctx->bit_factor == 4)
          out8[odx] = (uint8_t) result;
        else if (ctx->bit_factor == 2)
          out16[odx] = (uint16_t) result;
      }
    }
  }
#endif
  

  uint64_t out_data_size = data_size / (ctx->pol_factor * ctx->bit_factor);
  ctx->output.bytes_written += out_data_size;
  ctx->output.sample_offset += nsamp_in;

  if (ctx->output.bytes_written > ctx->output.block_size)
    multilog (log, LOG_ERR, "write_block_TFBP_to_BTF [%x] output block overrun by "
              "%"PRIu64" bytes\n", ctx->output.key, ctx->output.bytes_written - ctx->output.block_size);

  if (ctx->verbose > 1)
    multilog (log, LOG_INFO, "write_block_TFBP_to_BTF [%x] bytes_written=%"PRIu64", "
              "block_size=%"PRIu64"\n", ctx->output.key, ctx->output.bytes_written, ctx->output.block_size);

  // check if the output block is now full
  if (ctx->output.bytes_written >= ctx->output.block_size)
  {
    if (ctx->verbose > 1)
      multilog (log, LOG_INFO, "write_block_TFBP_to_BTF [%x] block now full bytes_written=%"PRIu64", block_size=%"PRIu64"\n", ctx->output.key, ctx->output.bytes_written, ctx->output.block_size);

    // check if this is the end of data
    if (client->transfer_bytes && ((ctx->bytes_in + data_size) == client->transfer_bytes))
    {
      if (ctx->verbose)
        multilog (log, LOG_INFO, "write_block_TFBP_to_BTF [%x] update_block_write written=%"PRIu64"\n", ctx->output.key, ctx->output.bytes_written);
      if (ipcio_update_block_write (ctx->output.hdu->data_block, ctx->output.bytes_written) < 0)
      {
        multilog (log, LOG_ERR, "write_block_TFBP_to_BTF [%x] ipcio_update_block_write failed\n", ctx->output.key);
         return -1;
      }
    }
    else
    {
      if (ctx->verbose > 1)
        multilog (log, LOG_INFO, "write_block_TFBP_to_BTF [%x] close_block_write written=%"PRIu64"\n", ctx->output.key, ctx->output.bytes_written);
      if (ipcio_close_block_write (ctx->output.hdu->data_block, ctx->output.bytes_written) < 0)
      {
        multilog (log, LOG_ERR, "write_block_TFBP_to_BTF [%x] ipcio_close_block_write failed\n", ctx->output.key);
        return -1;
      }
    }
    ctx->output.block_open = 0;
    ctx->output.bytes_written = 0;
  }
  else
  {
    if (ctx->output.bytes_written == 0)
      ctx->output.bytes_written = 1;
  }

  ctx->bytes_in += data_size;
  ctx->bytes_out += out_data_size;

  if (ctx->verbose > 1)
    multilog (log, LOG_INFO, "write_block_TFBP_to_BTF read %"PRIu64", wrote %"PRIu64" bytes\n", data_size, out_data_size);

  return data_size;
}

/*
 * compute the offset (mean) and scale (1/stddev) from the input data
 */
void compute_scales_offsets_TFBP (craft_dbtransposedb_t * ctx, float * in, const uint64_t nsamp)
{
  uint64_t isamp, idx=0;
  unsigned ichan, isig, ipol;

  for (isamp=0; isamp<nsamp; isamp++)
  { 
    for (ichan=0; ichan<ctx->nchan; ichan++)
    { 
      for (isig=0; isig<ctx->nsig; isig++)
      { 
        const unsigned ichansig= ichan * ctx->nsig + isig;
        double val = 0;
        for (ipol=0; ipol<ctx->npol; ipol++)
        {
          val += (double) in[idx];
          idx++;
        }
        ctx->sums[ichansig] += val;
        ctx->sums_sq[ichansig] += (val * val);
      }
    }
  }
    
  // update the total number of samples instragred in the the sums and sums_sq
  ctx->nsamps_integrated += nsamp;

  // now update the offset and scale
  for (ichan=0; ichan<ctx->nchan; ichan++)
  { 
    for (isig=0; isig<ctx->nsig; isig++)
    { 
      const unsigned ichansig = ichan * ctx->nsig + isig;
        
      const double mean = ctx->sums[ichansig] / ctx->nsamps_integrated;
      const double mean_sq = ctx->sums_sq[ichansig] / ctx->nsamps_integrated;
      const double variance = mean_sq - (mean * mean);
        
      ctx->offsets[ichansig] = (float) mean;
      if (variance == 0) 
        ctx->scales[ichansig] = 1.0f;
      else
        ctx->scales[ichansig] = (float) (1.0 / sqrt(variance));

      // also include the conversion to the right output bitwidth
      ctx->scales[ichansig] *= ctx->digi_scale;
    } 
  }
}


int main (int argc, char **argv)
{
  craft_dbtransposedb_t dbtransposedb;

  dada_client_t* client = 0;

  dada_hdu_t* hdu;

  /* DADA Logger */
  multilog_t* log = 0;

  /* Flag set in verbose mode */
  char verbose = 0;

  // number of transfers
  unsigned single_transfer = 0;

  // reblocking factor
  int reblock_factor = 64;

  // processing affinity
  int core = -1;

  // output quantization
  int nbit_out = 16;

  // number of samples to use for scale factors
  int nsamps_to_integrate = 2048;

  // input data block HDU key
  key_t in_key = 0;

  int arg = 0;

  while ((arg=getopt(argc,argv,"b:n:o:r:sv")) != -1)
  {
    switch (arg) 
    {
      case 'b':
        core = atoi(optarg);
        break;

      case 'n':
        nsamps_to_integrate = atoi(optarg);
        break;

      case 'o':
        nbit_out = atoi(optarg);
        break;

      case 'r':
        reblock_factor = atoi(optarg);
        break;

      case 's':
        single_transfer = 1;
        break;

      case 'v':
        verbose++;
        break;
        
      default:
        usage ();
        return 0;
      
    }
  }

  if (core >= 0)
    if (dada_bind_thread_to_core(core) < 0)
       fprintf (stderr, "failed to bind to core %d\n", core);

  craft_dbtransposedb_t * ctx = &dbtransposedb;

  ctx->digi_sigma = 6;
  // 8 bit values from Jenet and Anderson / digifil
  if (nbit_out == 8)
  {
    ctx->bit_factor = 4;
    ctx->digi_mean = 127.5;
    ctx->digi_scale = ctx->digi_mean / ctx->digi_sigma;
    ctx->digi_min = 0;
    ctx->digi_max = 255;
  }
  else if (nbit_out == 16)
  {
    ctx->bit_factor = 2;
    ctx->digi_mean = 32768;
    ctx->digi_scale = 512.0 / ctx->digi_sigma;
    ctx->digi_min = 0;
    ctx->digi_max = 65536;
  }
  else
  {
    fprintf (stderr, "ERROR: only 8 and 16 bit output supported\n");
    usage();
    return EXIT_FAILURE;
  }

  ctx->pol_factor = 2; // npol input==2, output== 1
  ctx->verbose = verbose;
  ctx->reblock_factor = (unsigned) reblock_factor;
  ctx->nsamps_to_integrate = nsamps_to_integrate;

  int num_args = argc-optind;
  if (num_args != 2)
  {
    fprintf(stderr, "craft_dbtransposedb: 2 arguments required, %d provided\n", num_args);
    usage();
    exit(EXIT_FAILURE);
  } 

  if (verbose)
    fprintf (stderr, "parsing input key=%s\n", argv[optind]);
  if (sscanf (argv[optind], "%x", &in_key) != 1) {
    fprintf (stderr, "craft_dbtransposedb: could not parse in key from %s\n", argv[optind]);
    return EXIT_FAILURE;
  }

  // read output DADA key from command line arguments
  if (verbose)
    fprintf (stderr, "parsing output key %s\n", argv[optind+1]);
  if (sscanf (argv[optind+1], "%x", &(ctx->output.key)) != 1) {
    fprintf (stderr, "craft_dbtransposedb: could not parse out key from %s\n", argv[optind+1]);
    return EXIT_FAILURE;
  }

  log = multilog_open ("craft_dbtransposedb", 0);

  multilog_add (log, stderr);

  if (verbose)
    multilog (log, LOG_INFO, "main: creating in hdu\n");

  // setup input DADA buffer
  hdu = dada_hdu_create (log);
  dada_hdu_set_key (hdu, in_key);
  if (dada_hdu_connect (hdu) < 0)
  {
    fprintf (stderr, "craft_dbtransposedb: could not connect to input data block\n");
    return EXIT_FAILURE;
  }

  if (verbose)
    multilog (log, LOG_INFO, "main: lock read key=%x\n", in_key);
  if (dada_hdu_lock_read (hdu) < 0)
  {
    fprintf(stderr, "craft_dbtransposedb: could not lock read on input data block\n");
    return EXIT_FAILURE;
  }

  // get the block size of the DADA data block
  uint64_t block_size = ipcbuf_get_bufsz ( (ipcbuf_t *) hdu->data_block);

  // setup output data block
  ctx->output.hdu = dada_hdu_create (log);
  dada_hdu_set_key (ctx->output.hdu, ctx->output.key);
  if (dada_hdu_connect (ctx->output.hdu) < 0)
  {
    multilog (log, LOG_ERR, "cannot connect to DADA HDU (key=%x)\n", ctx->output.key);
    return -1;
  }
  ctx->output.curr_block = 0;
  ctx->output.bytes_written = 0;
  ctx->output.block_open = 0;
  ctx->output.block_size = ipcbuf_get_bufsz ( (ipcbuf_t *) ctx->output.hdu->data_block);

  unsigned out_scale_factor = ctx->reblock_factor / (ctx->bit_factor * ctx->pol_factor);
  if (verbose)
    multilog (log, LOG_INFO, "main: ctx->output.block_size=%"PRIu64"\n", ctx->output.block_size);
  if (block_size * out_scale_factor != ctx->output.block_size)
  {
    multilog (log, LOG_ERR, "Block size mismatch, input=%lu, reblock factor=%u, output should be %lu, but was %lu\n", block_size, out_scale_factor, block_size * out_scale_factor, ctx->output.block_size);
   return EXIT_FAILURE;
  }

  client = dada_client_create ();

  client->log           = log;
  client->data_block    = hdu->data_block;
  client->header_block  = hdu->header_block;
  client->open_function = dbtransposedb_open;
  client->io_function   = dbtransposedb_write;
  client->io_block_function = dbtransposedb_write_block_TFBP_to_BTF;
  client->close_function = dbtransposedb_close;
  client->direction      = dada_client_reader;

  client->context = &dbtransposedb;
  client->quiet = (verbose > 0) ? 0 : 1;

  while (!client->quit)
  {
    if (verbose)
      multilog (log, LOG_INFO, "main: dada_client_read()\n");

    if (dada_client_read (client) < 0)
      multilog (log, LOG_ERR, "Error during transfer\n");

    if (verbose)
      multilog (log, LOG_INFO, "main: dada_hdu_unlock_read()\n");

    if (dada_hdu_unlock_read (hdu) < 0)
    {
      multilog (log, LOG_ERR, "could not unlock read on hdu\n");
      return EXIT_FAILURE;
    }

    if (single_transfer || ctx->quit)
      client->quit = 1;

    if (!client->quit)
    {
      if (dada_hdu_lock_read (hdu) < 0)
      {
        multilog (log, LOG_ERR, "could not lock read on hdu\n");
        return EXIT_FAILURE;
      }
    }
  }

  if (dada_hdu_disconnect (hdu) < 0)
    return EXIT_FAILURE;

  return EXIT_SUCCESS;
}
