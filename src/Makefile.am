
noinst_LTLIBRARIES = libcrafthd.la

libcrafthd_la_SOURCES = craft_utils.c

libcrafthd_la_headers = craft_utils.h

bin_PROGRAMS = 

AM_CFLAGS = -I.

LDADD = libcrafthd.la

if HAVE_PSRDADA

bin_PROGRAMS += craft_dbtransposedb
AM_CFLAGS += @PSRDADA_CFLAGS@
LDADD += @PSRDADA_LIBS@

if HAVE_SOFA

bin_PROGRAMS += craft_dbsigproc

AM_CFLAGS += -I$(top_builddir)/sigproc \
  @SOFA_CFLAGS@

LDADD += $(top_builddir)/sigproc/libsigproc.la \
  @SOFA_LIBS@

endif

endif

