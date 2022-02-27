////////////////////////////////////////////////////////////////////////////////////
//  Example program that shows how to use levmar in order to fit the four-
//  parameter exponential model
//  x_i,j = p[0]*exp( - ( (i-p[1])*(i-p[1]) + (j-p[2])*(j-p[2]) ) / 2 / sigma2) + p[3]
//  to a set of data measurements; example is based on a similar one from GSL.
//
//  Copyright (C) 2008-11  Manolis Lourakis (lourakis at ics forth gr)
//  Institute of Computer Science, Foundation for Research & Technology - Hellas
//  Heraklion, Crete, Greece.
//
//  Modified by P. Kestener (2021).
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
////////////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include <levmar.h>

#ifndef LM_DBL_PREC
#error Example program assumes that levmar has been compiled with double precision, see LM_DBL_PREC!
#endif


/* the following macros concern the initialization of a random number generator for adding noise */
#undef REPEATABLE_RANDOM /* #define this for repeatable program behavior across runs */
#define DBL_RAND_MAX (double)(RAND_MAX)

#ifdef _MSC_VER // MSVC
#include <process.h>
#define GETPID  _getpid
#elif defined(__GNUC__) // GCC
#include <sys/types.h>
#include <unistd.h>
#define GETPID  getpid
#else
#warning Do not know the name of the function returning the process id for your OS/compiler combination
#define GETPID  0
#endif /* _MSC_VER */

#ifdef REPEATABLE_RANDOM
#define INIT_RANDOM(seed) srandom(seed)
#else
#define INIT_RANDOM(seed) srandom((int)GETPID()) // seed unused
#endif

/* Gaussian noise with mean m and variance s, uses the Box-Muller transformation */
double gNoise(double m, double s)
{
  double r1, r2, val;

  r1=((double)random())/DBL_RAND_MAX;
  r2=((double)random())/DBL_RAND_MAX;

  val=sqrt(-2.0*log(r1))*cos(2.0*M_PI*r2);

  val=s*val+m;

  return val;
}

/* structure for passing user-supplied data to the objective function and its Jacobian */
struct xtradata
{

  int ni;
  int nj;
  double two_sigma2;

  /* more can be added here... */
};
typedef struct xtradata xtradata_t;

/* model to be fitted to measurements: x_i = p[0]*exp(-p[1]*i) + p[2], i=0...n-1 */
void gaussianFunc(double *p, double *x, int m, int n, void *data)
{
  xtradata_t *dat = (xtradata_t *) data;

  int idx = 0;
  for(int j=0; j<dat->nj; ++j)
  {
    for(int i=0; i<dat->ni; ++i)
    {
      x[idx] = p[0] * exp( - ( (i-p[1])*(i-p[1]) + (j-p[2])*(j-p[2]) ) / dat->two_sigma2) + p[3];
      idx++;
    }
  }

}

/* Jacobian of gaussianFunc() */
void jacGaussianFunc(double *p, double *jac, int m, int n, void *data)
{

  int i, j;
  xtradata_t *dat = (xtradata_t *) data;

  double sigma2 = dat->two_sigma2/2;

  /* fill Jacobian row by row */
  int idx = 0;
  int idx_jac = 0;

  for(int j=0; j<dat->nj; ++j)
  {
    for(int i=0; i<dat->ni; ++i)
    {
      double e = exp( - ( (i-p[1])*(i-p[1]) + (j-p[2])*(j-p[2]) ) / dat->two_sigma2);
      jac[idx_jac++] = e;
      jac[idx_jac++] = e * p[0] * (i-p[1])/sigma2;
      jac[idx_jac++] = e * p[0] * (j-p[2])/sigma2;
      jac[idx_jac++] = 1.0;
      idx++;
    }
  }
}

int main(int argc, char* argv)
{
  const int ni = 21;
  const int nj = 21;
  const int n=ni*nj, m=4; // 21^2 measurements, 4 parameters
  double p[m], x[n], opts[LM_OPTS_SZ], info[LM_INFO_SZ];
  int ret;

  double sigmaPSF = 0.66357;
  xtradata_t data;
  data.ni = ni;
  data.nj = nj;
  data.two_sigma2 = 2 * sigmaPSF * sigmaPSF;

  const double varTh = 1e-6;

  /* generate some measurement using the Gaussian model with
   * parameters (5.0, 10.35, 9.81, 0.4), corrupted with zero-mean
   * Gaussian noise of s=0.1
   */
  double p0 = 5.0;
  double p1 = 10.35;
  double p2 = 9.81;
  double p3 = 0.4;

  INIT_RANDOM(0);
  int idx = 0;
  for(int j=0; j<nj; ++j)
  {
    for(int i=0; i<ni; ++i)
    {
      x[idx] = (p0 * exp ( -( (i-p1)*(i-p1)+(j-p2)*(j-p2) )/data.two_sigma2) + p3) + gNoise(0.0, 0.1);
      idx++;
    }
  }

  /* initial parameters estimate: */
  p[0] = 4.0;
  p[1] = 10.0;
  p[2] = 10.0;
  p[3] = 0.0;

  /* optimization control parameters; passing to levmar NULL instead of opts reverts to defaults */
  opts[0] = LM_INIT_MU;
  opts[1] = 1E-15;
  opts[2] = 1E-15;
  opts[3] = 1E-20;
  opts[4] = LM_DIFF_DELTA; // relevant only if the finite difference Jacobian version is used

  /* invoke the optimization function */
  if (argc==1)
  {
    double *work, *covar;
    work = (double *) malloc((LM_DER_WORKSZ(m, n)+m*m)*sizeof(double));
    if(!work){
      fprintf(stderr, "memory allocation request failed in main()\n");
      exit(1);
    }
    covar=work+LM_DER_WORKSZ(m, n);

    clock_t t;
    t = clock();
    ret = dlevmar_der(gaussianFunc, jacGaussianFunc, p, x, NULL, varTh, m, n, 1000, opts, info, work, covar, (void *)&data); // with analytic Jacobian
    //ret=dlevmar_dif(expfunc, p, x, m, n, 1000, opts, info, NULL, NULL, (void *)&data); // without Jacobian
    t = clock() - t;
    double time_taken = ((double)t)/CLOCKS_PER_SEC; // in seconds

    printf("dlevmar_der() took %f micro-seconds to execute \n", 1e6*time_taken);

    printf("Levenberg-Marquardt LM_DER returned in %g iter, reason %g, sumsq %g [%g]\n", info[5], info[6], info[1], info[0]);
    printf("Best fit parameters: %.7g %.7g %.7g %.7g\n", p[0], p[1], p[2], p[3]);
  }
  else
  {
    double *work, *covar;
    work = (double *) malloc((LM_DIF_WORKSZ(m, n)+m*m)*sizeof(double));
    if(!work){
      fprintf(stderr, "memory allocation request failed in main()\n");
      exit(1);
    }
    covar=work+LM_DIF_WORKSZ(m, n);

    clock_t t;
    t = clock();
    ret = dlevmar_dif(gaussianFunc, p, x, NULL, varTh, m, n, 1000, opts, info, work, covar, (void *)&data); // without Jacobian
    t = clock() - t;
    double time_taken = ((double)t)/CLOCKS_PER_SEC; // in seconds

    printf("dlevmar_dif() took %f micro-seconds to execute \n", 1e6*time_taken);

    printf("Levenberg-Marquardt LM_DIF returned in %g iter, reason %g, sumsq %g [%g]\n", info[5], info[6], info[1], info[0]);
    printf("Best fit parameters: %.7g %.7g %.7g %.7g\n", p[0], p[1], p[2], p[3]);
  }

  return EXIT_SUCCESS;
}
