
#ifndef IVFGEO_PQ_INCLUDED
#define IVFGEO_PQ_INCLUDED

#include <utils/ivfpq.h>
#include <siftgeo/siftgeo.h>

#include "ivfgeo.h"


typedef struct {
  unsigned int scale:NB_BITS_VAL_SCALE;  /*!< Quantized scale of the region */
  int angle:NB_BITS_VAL_ANGLE;           /*!< Quantized angle */
} PACKTHIS ivfgeo_pq_point_t;


typedef struct {

  ivfpq_t *ivfpq;
  
  int nim;
  
  int norm_type;                /*!< modify the default normalization scheme for the frequency vector normalization */
  int na_img;               /* allocates size of im-related arrays */
  float *norms;             /*!< norms associated with the images (allocated size:na_img) */
  int *im_ends;

  long tot_pts; 
  ivfgeo_pq_point_t *points;  
  long na_pts;

  int tfidf;                    /*!< define the kind of tf-idf scheme which is applied */
  float *tfidf_weights;     /*!< tf-idf weights (allocated size: nbvw) */


  /* query params (not stored)  */
  int ma;
  double ma_disratio;
  int k;
  int nthread;
  int dist_w_type;
  double sigma; /* may be used in dist_w_type */
  int wgc_type;                 /*!< method to extract maximum from the angle / scale histogram(s) */
  float *wgc_weights;       /*!< a-priori weighting of angle and scale differences */

  int burstiness;

  double scale_w;
  double * scale_weights;

} ivfgeo_pq_t;




/* ownership of ivf is transferred! */
ivfgeo_pq_t *ivfgeo_pq_new(ivfpq_t *ivfpq); 

/* used only on an empty ivf for add & merge */
ivfgeo_pq_t *ivfgeo_pq_dup(ivfgeo_pq_t * ivf);

/* ivf2 dealloc'ed on output */
void ivfgeo_pq_merge(ivfgeo_pq_t * ivf,ivfgeo_pq_t *ivf2);

void ivfgeo_pq_display(const ivfgeo_pq_t * ivf); 

/* returns label */
int ivfgeo_pq_add (ivfgeo_pq_t * ivf, const pointset_t * ps);


/* precomputed vw's */
int ivfgeo_pq_add_with_vw (ivfgeo_pq_t * ivf, const pointset_t * ps, const int *vw);



void ivfgeo_pq_query (const ivfgeo_pq_t * ivf,
                      const pointset_t * vw,
                      float *dists);

 
void ivfgeo_pq_compute_tfidf (ivfgeo_pq_t * ivf);
void ivfgeo_pq_compute_norms (ivfgeo_pq_t * ivf);
void ivfgeo_pq_compute_scale_weights(ivfgeo_pq_t * ivf); 
void ivfgeo_pq_set_wgc_type (ivfgeo_pq_t * ivf, int wgc_type);




void ivfgeo_pq_fwrite(const ivfgeo_pq_t * ivf, FILE *f);

ivfgeo_pq_t *ivfgeo_pq_fread(FILE *f,float *centroids);



#endif
