#include <cmath>
#include <cstdlib>
#include <cassert>
#include <cstdio>
#if defined (HAVE_GSL)
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_sf_gamma.h>

/* #ifdef __USE_ISOC99 */
/* #error */
/* #endif */
/*
#ifdef _BSD_SOURCE
#error
#endif
*/
void chisq2x2(const double *x, size_t size, double *chisq, double *chisqYates, double *p, double *pYates)
     /*chi-squared test on a 2x2 table with df=1 */
     /*See Sokal and Rohlf "Biometry", page 736 */
{
  double n = x[0]+x[1]+x[2]+x[3];
  double a = x[0];
  double b = x[1];
  double c = x[2];
  double d = x[3];
  assert(size==4);
  *chisq = (std::pow( a*d - b*c , 2. )*n)/( (a+b)*(c+d)*(a+c)*(b+d) );
  *chisqYates = (std::pow( fabs(a*d-b*c) - (n/2.) , 2 )*n )/( (a+b)*(c+d)*(a+c)*(b+d) );
  /*p-value is from an incomplete gamma function with parameters (0.5*DF,0.5+chi-squared)*/
  if(std::isfinite(*chisq))
    *p = gsl_sf_gamma_inc_Q(0.5,0.5*(*chisq)); 
  else
    *p = strtod("NAN",NULL);
  if(std::isfinite(*chisqYates))
    *pYates = gsl_sf_gamma_inc_Q(0.5,0.5*(*chisqYates)); 
  else
    *pYates = strtod("NAN",NULL);
}

#endif
