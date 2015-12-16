/*
 * CPU timer
 */
#include <stdio.h>
#include <sys/time.h>

#include "timer.h"

static struct timeval start, stop;

void timer_start (void)
{
  gettimeofday(&start, (struct timezone *) NULL);
}

double timer_stop (void)
{
  gettimeofday( &stop, (struct timezone *) NULL);
  return ((double)(stop.tv_sec - start.tv_sec) + (double)((stop.tv_usec - start.tv_usec)) * 1.e-6);
}
