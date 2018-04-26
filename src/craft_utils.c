/***************************************************************************
 *  
 *    Copyright (C) 2018 by Andrew Jameson
 *    Licensed under the Academic Free License version 2.1
 * 
 ****************************************************************************/

#include "craft_utils.h"

#ifdef HAVE_SOFA
#include "sofa.h"
#endif

#include "stdlib.h"
#include "stdio.h"
#include "string.h"

#ifdef HAVE_SOFA
double crafthd_mjd_from_utc (struct tm* utc)
{
  int iy = 1900 + utc->tm_year;
  int im = 1 + utc->tm_mon;
  int id = (int) utc->tm_mday;

  double d, d1, d2;

  int rval = iauCal2jd ( iy, im, id, &d1, &d2 );
  if (rval != 0)
    return -1;

  int ihour = (int) utc->tm_hour;
  int imin  = (int) utc->tm_min;
  double sec = (double) utc->tm_sec;

  rval = iauTf2d ('+', ihour, imin, sec, &d );
  if (rval != 0)
    return -1;

  double mjd = d2 + d;

  return mjd;
}
#endif

/**
 *  Configure HH:MM:SS string to a sigproc double
 */
int crafthd_hhmmss_to_sigproc (char * hhmmss, double * sigproc)
{
  int ihour = 0;
  int imin = 0;
  double sec = 0;
  const char *sep = ":";
  char * saveptr;

  char * copy = (char *) malloc (strlen(hhmmss) + 1);
  strcpy (copy, hhmmss);

  char * str = strtok_r(copy, sep, &saveptr);
  if (str != NULL)
  {
    if (sscanf(str, "%d", &ihour) != 1)
      return -1;

    str = strtok_r(NULL, sep, &saveptr);
    if (str != NULL)
    {
      if (sscanf(str, "%d", &imin) != 1)
        return -1;

      str = strtok_r(NULL, sep, &saveptr);
      if (str != NULL)
      {
        if (sscanf(str, "%lf", &sec) != 1)
          return -1;
      }
    }
  }
  free (copy);

  //char s = '\0';
  if (ihour < 0)
  {
    ihour *= -1;
    //s = '-';
  }

  *sigproc = ((double)ihour*1e4  + (double)imin*1e2 + sec);
  return 0;
}

/**
 *  Configure DD:MM:SS string to a sigproc double
 */
int crafthd_ddmmss_to_sigproc (char * ddmmss, double * sigproc)
{
  int ideg = 0;
  int iamin = 0;
  double asec = 0;
  const char *sep = ":";
  char * saveptr;

  char * copy = (char *) malloc (strlen(ddmmss) + 1);
  strcpy (copy, ddmmss);

  char * str = strtok_r(ddmmss, sep, &saveptr);
  if (str != NULL)
  {
    if (sscanf(str, "%d", &ideg) != 1)
      return -1;

    str = strtok_r(NULL, sep, &saveptr);
    if (str != NULL)
    {
      if (sscanf(str, "%d", &iamin) != 1)
        return -1;

      str = strtok_r(NULL, sep, &saveptr);
      if (str != NULL)
      {
        if (sscanf(str, "%lf", &asec) != 1)
          return -1;
      }
    }
  }

  free (copy);

  if (ideg < 0)
    *sigproc = ((double) ideg*1e4 - (double) iamin*1e2) - asec;
  else
    *sigproc = ((double) ideg*1e4  + (double) iamin*1e2) + asec;

  return 0;
}
