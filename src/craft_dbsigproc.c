//
// Write multi-beam data straight to a sigproc file
//

// PSRDada includes
#include "dada_affinity.h"
#include "dada_client.h"
#include "dada_hdu.h"
#include "dada_def.h"
#include "ascii_header.h"
#include "tmutil.h"

// craft heimdall includes
#include "sigproc/filterbank.h"
#include "craft_utils.h"

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
#include <emmintrin.h>
#include <stdint.h>

int64_t dbsigproc_write_block (dada_client_t *, void *, uint64_t, uint64_t);

void usage()
{
  fprintf (stdout,
           "craft_dbsigproc [options] start_beam total_nbeam\n"
           " -b <core>    bind compuation to CPU core\n"
           " -k <key>     DADA key for input\n"
           " -r <nbit>    requantise from 8-bit to n-bit\n"
           " -s           1 transfer, then exit\n"
           " -v           verbose mode\n"
           " start_beam   number of first beam in block\n"
           " total_nbeam  total number of beams in observation\n");
}

typedef struct {

  unsigned int start_beam;

  unsigned int total_nbeam;

  // number of bytes read
  uint64_t bytes_in;

  // verbose output
  int verbose;

  unsigned int nbeam;

  unsigned int nchan;

  unsigned int nbit;

  unsigned int ndim;

  unsigned int nbit_out;

  float freq;

  float bw;

  float tsamp;

  unsigned quit;

  char order[4];

  FILE ** fptrs;

  uint64_t block_size;

  void * out;

  uint64_t out_bufsz;

  char lsb;

} craft_dbsigproc_t;

#define DADA_DBSIGPROC_INIT { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, "", 0, 0, 0, 0, 0 }

/*! Function that opens the data transfer target */
int dbsigproc_open (dada_client_t* client)
{
  // the craft_dbsigproc specific data
  craft_dbsigproc_t* ctx = (craft_dbsigproc_t *) client->context;

  // status and error logging facilty
  multilog_t* log = client->log;

  // header to copy from in to out
  //char * header = 0;
 
  unsigned i = 0;

  if (ctx->verbose)
    multilog (log, LOG_INFO, "dbsigproc_open()\n");

  char utc_start_str[20];
  if (ascii_header_get (client->header, "UTC_START", "%s", utc_start_str) != 1)
  {
    multilog (log, LOG_INFO, "open: header with no UTC_START\n");
    return -1;
  }
  // convert to a unixtime
  const time_t utc_start = str2utctime (utc_start_str);
  struct tm * utc_t = gmtime (&utc_start);

  // start time of the observation in MJD
  tstart = crafthd_mjd_from_utc (utc_t);

  // get the transfer size (if it is set)
  int64_t transfer_size;
  if (ascii_header_get (client->header, "TRANSFER_SIZE", "%"PRIi64, &transfer_size) != 1)
  {
    transfer_size = 0;
  }

  int64_t file_size;
  if (ascii_header_get (client->header, "FILE_SIZE", "%"PRIi64, &file_size) != 1)
  {
    file_size = 0;
  }

  uint64_t obs_offset;
  if (ascii_header_get (client->header, "OBS_OFFSET", "%"PRIu64, &obs_offset) != 1)
  {
    multilog (log, LOG_ERR, "open: header with no OBS_OFFSET\n");
    return -1;
  }

  uint64_t resolution;
  if (ascii_header_get (client->header, "RESOLUTION", "%"PRIu64, &resolution) != 1)
  {
    multilog (log, LOG_WARNING, "open: header with no RESOLUTION\n");
    resolution = 0;
  }

  if (ascii_header_get (client->header, "NBEAM", "%u", &(ctx->nbeam)) != 1)
  {
    multilog (log, LOG_WARNING, "open: header with no NBEAM\n");
    ctx->nbeam = 1;
  }

  if (ascii_header_get (client->header, "NCHAN", "%u", &(ctx->nchan)) != 1)
  {           
    multilog (log, LOG_ERR, "open: header with no NCHAN\n");
    return -1;                
  }

  if (ascii_header_get (client->header, "FREQ", "%f", &(ctx->freq)) != 1)
  {
    multilog (log, LOG_ERR, "open: header with no FREQ\n");
    return -1;
  }

  if (ascii_header_get (client->header, "BW", "%f", &(ctx->bw)) != 1)
  {     
    multilog (log, LOG_ERR, "open: header with no BW\n");
    return -1;    
  }

  // check if lower side band
  if (ctx->bw < 0)
    ctx->lsb = 1;
  else
    ctx->lsb = 0;

  if (ascii_header_get (client->header, "TSAMP", "%f", &(ctx->tsamp)) != 1)
  {
    multilog (log, LOG_ERR, "open: header with no TSAMP\n");
    return -1;
  }

  if (ascii_header_get (client->header, "NBIT", "%u", &(ctx->nbit)) != 1)
  {
    multilog (log, LOG_ERR, "open: header with no NBIT\n");
    return -1;
  }
  if ((ctx->nbit != 8) && (ctx->nbit != 16))
  {
    multilog (log, LOG_ERR, "input NBIT=%u, must be 8 or 16\n", ctx->nbit);
    return -1;
  }
  //if (ctx->nbit == 32)
  //  ctx->nbit_out = 32;

  if (ascii_header_get (client->header, "NDIM", "%u", &(ctx->ndim)) != 1)
  {           
    multilog (log, LOG_ERR, "open: header with no NDIM\n");
    return -1;                
  }                             

  if (ascii_header_get (client->header, "ORDER", "%s", &(ctx->order)) != 1)
  {
    multilog (log, LOG_ERR, "open: header with no ORDER\n");
    return -1;
  }
  else
  {
    if (strcmp(ctx->order, "STF") != 0)
    {
      multilog (log, LOG_ERR, "require STF input ordering");
      return -1;
    }
  }
   
  // at this point we fill in the SIGPROC global variables (yuck)
  // source_name -> char[80]

  if (ascii_header_get (client->header, "SOURCE", "%s", source_name) != 1)
  {
    multilog (log, LOG_ERR, "open: header with no SOURCE\n");
    return -1;
  }

  machine_id = 0;   // 
  telescope_id = 0; // AKA Fake

  char buffer[128];
  if (ascii_header_get (client->header, "RA", "%s", buffer) != 1)
  {
    multilog (log, LOG_ERR, "open: header with no RAJ\n");
    return -1;
  }
  if (crafthd_hhmmss_to_sigproc (buffer, &src_raj) < 0)
  {
    multilog (log, LOG_ERR, "open: could not parse RAJ from %s\n", buffer);
    return -1;
  }

  if (ascii_header_get (client->header, "DEC", "%s", buffer) != 1)
  { 
    multilog (log, LOG_ERR, "open: header with no DEC\n");
    return -1;
  }
  if (crafthd_ddmmss_to_sigproc (buffer, &src_dej) < 0)
  {
    multilog (log, LOG_ERR, "open: could not parse DEC from %s\n", buffer);
    return -1;
  }

  az_start = 0;
  za_start = 0;

  // channel_bandwidth (convert to LSB if required)
  if (ctx->lsb)
    foff = ctx->bw / ctx->nchan;
  else
    foff = -1 * (ctx->bw / ctx->nchan);

  // centre frequnecy of channel (0) (highest channel)
  if (ctx->lsb)
    fch1 = ctx->freq - (ctx->bw/2) - (foff/2);
  else
    fch1 = ctx->freq + (ctx->bw/2) + (foff/2);

  nchans = ctx->nchan;

  // npol
  nifs = 1;

  // number of output bits
  obits = ctx->nbit_out;

  // sampling time
  tsamp = ctx->tsamp / 1e6;
  
  ifstream[0] = 'Y';

  ctx->fptrs = (FILE **) malloc (sizeof (FILE *) * ctx->nbeam);

  nbeams = ctx->total_nbeam;

  // open filterbank output files for each beam
  for (i=0; i<ctx->nbeam; i++)
  {
    ibeam = i;

    sprintf (buffer, "mkdir BEAM_%03d", ctx->start_beam + i);
    system (buffer);

    // write ascii header to file 
    sprintf (buffer, "BEAM_%03d/obs.header", ctx->start_beam + i);
    ctx->fptrs[i] = fopen (buffer, "w");
    size_t header_len = strlen (client->header);
    client->header[header_len] = '\n';
    header_len++;
    client->header[header_len] = '\0';
    fwrite (client->header, header_len, 1, ctx->fptrs[i]);
    fclose (ctx->fptrs[i]);

    // sigproc filename
    sprintf (buffer, "BEAM_%03d/%s.fil", ctx->start_beam + i, utc_start_str);

    if (ctx->verbose)
      multilog (log, LOG_INFO, "open: opening file %s\n", buffer);

    // open a file for each beam
    ctx->fptrs[i] = fopen (buffer, "w");

    // write the header to the file pointer
    filterbank_header (ctx->fptrs[i]);

    // TODO update the source positions for each beam!
  }

  ctx->out_bufsz = ctx->block_size / (ctx->nbeam * (ctx->nbit / ctx->nbit_out));
  ctx->out = malloc (ctx->out_bufsz);
	multilog (log, LOG_INFO, "open:ctx->out_bufsz=%lu\n", ctx->out_bufsz);

  client->transfer_bytes = transfer_size; 
  client->optimal_bytes = 64*1024*1024;

  ctx->bytes_in = 0;
  client->header_transfer = 0;

  return 0;
}

int dbsigproc_close (dada_client_t* client, uint64_t bytes_written)
{
  craft_dbsigproc_t* ctx = (craft_dbsigproc_t*) client->context;
  
  multilog_t* log = client->log;

  unsigned i = 0;

  if (ctx->verbose)
    multilog (log, LOG_INFO, "close: bytes_read=%"PRIu64"\n", ctx->bytes_in);

  if (ctx->fptrs)
  {
    for (i=0; i<ctx->nbeam; i++)
    {  
      if (ctx->fptrs[i])
      {
        if (ctx->verbose > 1)
          multilog (log, LOG_INFO, "close: closing file ptr for beam %d\n", i);
        fclose (ctx->fptrs[i]);
      }
      ctx->fptrs[i] = 0;
    }
    free (ctx->fptrs);
  }
  ctx->fptrs = 0;

  return 0;
}

/*! Pointer to the function that transfers data to/from the target */
int64_t dbsigproc_write (dada_client_t* client, void* data, uint64_t data_size)
{
  craft_dbsigproc_t* ctx = (craft_dbsigproc_t*) client->context;

  multilog_t * log = client->log;

  if (ctx->verbose)
    multilog (log, LOG_INFO, "write: to_write=%"PRIu64"\n", data_size);

  // write data to all data blocks

  ctx->bytes_in += data_size;

  if (ctx->verbose)
    multilog (log, LOG_INFO, "write: read %"PRIu64", wrote %"PRIu64" bytes\n", data_size, data_size);
 
  return data_size;
}

int64_t dbsigproc_write_block_STF_to_TF (dada_client_t* client, void* in_data, uint64_t in_data_size, uint64_t block_id)
{
  craft_dbsigproc_t* ctx = (craft_dbsigproc_t*) client->context;

  multilog_t * log = client->log;

  if (ctx->verbose > 1)
    multilog (log, LOG_INFO, "write_block_STF_to_TF: data_size=%"PRIu64", block_id=%"PRIu64"\n",
              in_data_size, block_id);

  const uint64_t out_data_size = in_data_size / ctx->nbeam;
  const uint64_t nsamp = out_data_size / (ctx->ndim * (ctx->nbit/8) * ctx->nchan);

  if (in_data_size < ctx->block_size)
  {
    multilog (log, LOG_INFO, "write_block_STF_to_TF: in_data_size[%lu] < in_block_size[%lu]\n",
              in_data_size, ctx->block_size);
    return in_data_size;
  }

  if (ctx->verbose > 1)
    multilog (log, LOG_INFO, "write_block_STF_to_TF: in_data_size=%"PRIu64", "
              "out_data_size=%"PRIu64" nsamp=%"PRIu64"\n", in_data_size, out_data_size, nsamp);

  unsigned i, ichan, isamp;

  if ((ctx->nbit == 8) && (ctx->nbit_out == 8))
  {
    uint8_t * in = (uint8_t *) in_data;
    uint8_t * out = (uint8_t *) ctx->out;
	  uint64_t idx = 0;
	  uint64_t odx = 0;

    for (i=0; i<ctx->nbeam; i++)
    {
	  	odx = 0;
      for (isamp=0; isamp<nsamp; isamp++)
	    {
      	for (ichan=0; ichan<ctx->nchan; ichan++)
      	{
        	out[odx] = in[idx];
					idx++;
					odx++;
        }
      }
      fwrite (out, ctx->out_bufsz, 1, ctx->fptrs[i]);
    }
  }
  else if ((ctx->nbit == 16) && (ctx->nbit_out == 16))
  {
    uint16_t * in = (uint16_t *) in_data;
    uint16_t * out = (uint16_t *) ctx->out;
    uint64_t idx = 0;
    uint64_t odx = 0;

    for (i=0; i<ctx->nbeam; i++)
    {
      odx = 0;
      for (isamp=0; isamp<nsamp; isamp++)
      {
        for (ichan=0; ichan<ctx->nchan; ichan++)
        {
          out[odx] = in[idx];
          idx++;
          odx++;
        }
      }
      fwrite (out, ctx->out_bufsz, 1, ctx->fptrs[i]);
    }
  }
  else
  {
    multilog (log, LOG_ERR, "unsupported bit width input=%d output=%d\n", ctx->nbit, ctx->nbit_out);
  }

  ctx->bytes_in += in_data_size;

  if (ctx->verbose > 1)
    multilog (log, LOG_INFO, "write_block_STF_to_TF read %"PRIu64" bytes\n", in_data_size);

  return in_data_size;
}



int main (int argc, char **argv)
{
  craft_dbsigproc_t dbsigproc = DADA_DBSIGPROC_INIT;

  dada_hdu_t* hdu = 0;

  dada_client_t* client = 0;

  /* DADA Logger */
  multilog_t* log = 0;

  /* Flag set in verbose mode */
  char verbose = 0;

  // number of transfers
  unsigned single_transfer = 0;

  // input data block HDU key
  key_t key = DADA_DEFAULT_BLOCK_KEY;

  int core = -1;

  int arg = 0;

  // by default output 8-bit data
  int nbit_out = 8;

  while ((arg=getopt(argc,argv,"b:k:hr:sv")) != -1)
  {
    switch (arg) 
    {
      case 'b':
        if (optarg)
        {
          core = atoi(optarg);
          break;
        }
        else
        {
          fprintf(stderr, "ERROR: -b requires argument\n");
          usage();
          return EXIT_FAILURE;
        }

      case 'k':
        if (optarg)
        {
          if (sscanf (optarg, "%x", &key) != 1) 
          {
            fprintf (stderr, "craft_dbsigproc: could not parse in key from %s\n", optarg);
            return EXIT_FAILURE;
          }
        }
        else
        {
          fprintf (stderr, "ERROR: -k requires argument\n");
          return EXIT_FAILURE;
        }
        break;

      case 'h':
        usage();
        return EXIT_SUCCESS;

      case 'r':
        nbit_out = atoi(optarg);
        if (!((nbit_out == 8) || (nbit_out == 16)))
        {
          fprintf (stderr, "ERROR: output nbit must be 8 or 16\n");
          return EXIT_FAILURE;
        }
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
      multilog(log, LOG_WARNING, "main: failed to bind to core %d\n", core);

  dbsigproc.verbose = verbose;
  dbsigproc.nbit_out = nbit_out;

  int num_args = argc - optind;
  if (num_args != 2)
  {
    fprintf(stderr, "craft_dbsigproc: 2 command line args required\n");
    usage();
    exit(EXIT_FAILURE);
  } 

  if (verbose)
    fprintf (stderr, "parsing start beam=%s\n", argv[optind]);
  if (sscanf (argv[optind], "%u", &(dbsigproc.start_beam)) != 1) {
    fprintf (stderr, "craft_dbsigproc: could not parse start beam from %s\n", argv[optind]);
    return EXIT_FAILURE;
  }

  if (verbose)
    fprintf (stderr, "parsing totai nbeam=%s\n", argv[optind+1]);
  if (sscanf (argv[optind+1], "%u", &(dbsigproc.total_nbeam)) != 1) {
    fprintf (stderr, "craft_dbsigproc: could not parse total nbeam from %s\n", argv[optind+1]);
    return EXIT_FAILURE;
  }

  log = multilog_open ("craft_dbsigproc", 0);
  multilog_add (log, stderr);

  if (verbose)
    multilog (log, LOG_INFO, "main: creating in hdu\n");

  // setup input DADA buffer
  hdu = dada_hdu_create (log);
  dada_hdu_set_key (hdu, key);
  if (dada_hdu_connect (hdu) < 0)
  {
    fprintf (stderr, "craft_dbsigproc: could not connect to input data block\n");
    return EXIT_FAILURE;
  }

  if (verbose)
    multilog (log, LOG_INFO, "main: lock read key=%x\n", key);
  if (dada_hdu_lock_read (hdu) < 0)
  {
    fprintf(stderr, "craft_dbsigproc: could not lock read on input data block\n");
    return EXIT_FAILURE;
  }

  // get the block size of the DADA data block
  dbsigproc.block_size = ipcbuf_get_bufsz ( (ipcbuf_t *) hdu->data_block);

  client = dada_client_create ();

  client->log           = log;
  client->data_block    = hdu->data_block;
  client->header_block  = hdu->header_block;
  client->open_function = dbsigproc_open;
  client->io_function   = dbsigproc_write;

  client->io_block_function = dbsigproc_write_block_STF_to_TF;

  client->close_function = dbsigproc_close;
  client->direction      = dada_client_reader;

  client->context = &dbsigproc;
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

    if (single_transfer || dbsigproc.quit)
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
