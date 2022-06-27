

#include "siftgeo_and_nn.h"
#include <yael/vector.h>





float *siftgeo_to_fvecs (const point_t * corners, int n, int *d_out)
{

  if (n == 0)
    return NULL;

  int d = corners[0].dim;

  float *points = fvec_new (d * n);
  
  const pointset_t ps={
    n,
    corners};

  pointset_into_fvecs (&ps, points);

  *d_out = d;

  return points;
}


int pointset_into_fvecs (const pointset_t *ps,float *f) {

  int i, j;  
  
  if (ps->n == 0)
    return 0;

  int d = ps->pts[0].dim;

  for (i = 0; i < ps->n; i++) {
    const unsigned char *il = ps->pts[i].desc;
    float *ol = f + i * d;
    for (j = 0; j < d; j++)
      ol[j] = il[j];
  }
  
  return d;  
}


