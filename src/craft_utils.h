/***************************************************************************
 *  
 *    Copyright (C) 2018 by Andrew Jameson
 *    Licensed under the Academic Free License version 2.1
 * 
 *****************************************************************************/

#ifndef __CRAFT_UTILS_H
#define __CRAFT_UTILS_H

#include "config.h"

#include "time.h"

#ifdef __cplusplus
extern "C" {
#endif

#ifdef HAVE_SOFA
// convert a struct tm * UTC into an MJD
double crafthd_mjd_from_utc (struct tm * utc);
#endif

// convert HH:MM:SS char string to sigproc double
int crafthd_hhmmss_to_sigproc (char * hhmmss, double * sigproc);

// convert DD:MM:SS char string to sigproc double
int crafthd_ddmmss_to_sigproc (char * ddmmss, double * sigproc);

#ifdef __cplusplus
}
#endif

#endif
