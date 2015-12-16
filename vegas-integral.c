#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>
#include "timer.c"
#include "timer.h"

double g (double *t, size_t dim, void *params);

void timer_start(void);
double timer_stop(void);

int main (void)
{
  double tvegas, tmonte;
  double res, err;
  double distmin = 1.001;
  double distmax = 4.;
  double dist = distmin;
  int np = 20;
  double nt = (distmax - distmin) / ((double) np - 1.);
  size_t dim = 6;
  double x1[] = { 0., 0., 0., 0., 0., 0. };
  double xu[] = { 1., 1., 1., 1., 1., 1. };

  double vegas[20], dipole[20], distance[20], monte[20];
  int f = 0;
  gsl_rng *r = gsl_rng_alloc (gsl_rng_taus2);
  unsigned long seed = 1UL;

  gsl_rng_set (r, seed);

  size_t calls = 1500000;

  //Vegas Integration

  gsl_monte_function G = { &g, dim, &dist };

  gsl_monte_vegas_state *sv = gsl_monte_vegas_alloc (dim);

  gsl_monte_vegas_init (sv);

  timer_start();

  do
    {

      gsl_monte_vegas_integrate (&G, x1, xu, dim, calls / 5, r, sv, &res,
				 &err);
      do
	{
	  gsl_monte_vegas_integrate (&G, x1, xu, dim, calls, r, sv, &res,
				     &err);
	  fflush (stdout);
	}
      while (fabs (gsl_monte_vegas_chisq (sv) - 1.) > .2);
      vegas[f] = res;
      f++;
      dist = dist + nt;
    }
  while (dist <= distmax);

  tvegas = timer_stop();

  gsl_monte_vegas_free (sv);

  // Montecarlo Integration

  double sum = 0.;
  double x[6];
  long i, j, l;

  l = 1E6;
  dist = distmin;
  
  timer_start();

  for (j = 0; j < np; j++)
    {
      sum = 0.;
      for (i = 0; i < l; i++)
	{
	  for (int k = 0; k < (int) dim; k++)
	    {
	      x[k] = gsl_rng_uniform (r);
	    }
	  sum += g (x, dim, &dist);
	}

      res = sum / (double) l;
      monte[j] = res;
      dist += nt;
    }

  tmonte =  timer_stop();

  gsl_rng_free (r);
  dist = distmin;

  for (int y = 0; y < np; y++)
    {
      dipole[y] = 2. / pow (dist, 3);
      distance[y] = dist;
      dist = dist + nt;
    }

  printf
    ("#    Dist               Vegas           Monte       Dipolapprox\n");
  for (int q = 0; q < np; q++)
    {
      double dd = distance[q];
      double vv = fabs (vegas[q]);
      double hh = fabs (monte[q]);
      double di = fabs (dipole[q]);
      printf ("   %.6f         %.6f       %.6f      %.6f\n", dd, vv, hh, di);
    }
  
  printf("\n\n");
  printf("Time for Vegas  = %f\n", tvegas);
  printf("Time for Monte  = %f\n", tmonte);
  printf("Speed up Factor = %f\n", tvegas/tmonte);

  return 0;
}
