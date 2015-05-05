#if defined (HAVE_GSL)

#include <math.h>
#include <stdio.h>
#include <gsl/gsl_sf_gamma.h>

void g2(double *x, size_t size, double *G, double *correction, double *pG,double *pGcorr)
     /*G test for a 2x2 table w/William's correction*/
{
  double a = x[0];
  double b = x[1];
  double c = x[2];
  double d = x[3];
  double n = a+b+c+d;  
  double sumLnFreq = a*log(a)+b*log(b)+c*log(c)+d*log(d);
  double sumLnTots = (a+b)*log(a+b) + (c+d)*log(c+d) + (a+c)*log(a+c) + (b+d)*log(b+d);
  *G = 2.*(sumLnFreq - sumLnTots + n*log(n));/* the G statistic itself */
  if(G*<0)*G=0.;
  /* correction = William's correction */
  *correction =1.+( ((n/(a+b))+(n/(c+d))-1.)*((n/(a+c))+(n/(b+d))-1.))/(6*n);
  *pG = gsl_sf_gamma_inc_Q(0.5,0.5*(*G));/* probability is chi-squared with df = 1 */
  *pGcorr = gsl_sf_gamma_inc_Q(0.5,0.5*(*G)/(*correction));
}

#endif
