#ifdef __cplusplus
extern "C"
{
#include<cstdlib>
#else
#include<stdlib.h>
#endif

#if defined (HAVE_GSL)
  void g2(double *x, size_t size, double *G, double *correction, double *pG,double *pGcorr);
#endif

#ifdef __cplusplus
}
#endif
