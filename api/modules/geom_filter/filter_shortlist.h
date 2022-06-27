/*---------------------------------------------------------------------------*/

#ifndef FILTER_SHORTLIST_H_INCLUDED
#define FILTER_SHORTLIST_H_INCLUDED



#include <siftgeo/siftgeo.h>



/* Minimalistic image descriptor                                             */
typedef pointset_t imdesc_t;

/*---------------------------------------------------------------------------*/
/* Small set of searchable images                                            */
/*---------------------------------------------------------------------------*/

typedef struct shortlist_t shortlist_t;

/* the imdesc_t*'s must remain valid during the lifespan of the shortlist_t */
shortlist_t* shortlist_new(imdesc_t **imdescs,int nimdesc);

void shortlist_delete(shortlist_t *);

/*---------------------------------------------------------------------------*/

/* match between a point on nthe query image and one from a database image */
typedef struct pointmatch_t {
  struct pointmatch_t *next;
  float score;  /* in [0,1] bigger = better */
  int qpt,dbpt; /* point numbers in their respective image descriptors */
} pointmatch_t;

typedef struct {

  /* internal */
  struct fast_alloc_t *fa_pointmatch;

  /* ---------- search results */

  /* linked list of point matches for this image */
  pointmatch_t *ptmatches;
  
  /* query in lowDim dimensions */
  int stage1_nmatch;

  /* in 128 dim */
  int stage2_nmatch;
  double stage2_votes;

  /* after geometric verification */
  int stage3_nmatch;
  double stage3_votes;

  /* new improved voting system */
  double final_votes;

  /* affine matrix (defined iff stage3_nmatch>0) */
  double affine[6];
  
  int* stage3_ids;
  
} imagematch_t;

imagematch_t* imagematches_new(int n);

  // added
int count_ptmatches(pointmatch_t *pm);  
  
/*---------------------------------------------------------------------------*/

/*! @brief Do exhaustive search
 *
 * method==0: exhaustive search of pts within thr (thr=0.2 typical for SIFT, 0.3 for CS-LBP) \n
 * method==1: take k nearest neighbours in each db image (k=(int)thr)  \n
 * method==2: take k nearest neighbours in union of db images (k=(int)thr)
 * 
 * return array of imagematch structures that must be free'd with free_imagematch
 */
imagematch_t *shortlist_match_points_exact(shortlist_t *,imdesc_t *query,
                                           int method,double thr,
                                           void (*peek_fun)(void *arg,double frac),
                                           void *peek_arg);

/*! @brief Query and db images should be vwgeo's */
imagematch_t *shortlist_match_points_vw(shortlist_t *,imdesc_t *query,
                                        int hamming_thresh,
                                        void (*peek_fun)(void *arg,double frac),
                                        void *peek_arg);


void imagematch_align_vw(imdesc_t *qim,
                         imdesc_t *dbim,
                         imagematch_t *imm,
                         int hamming_thresh);

/*! free n imagematch_t's */
void imagematches_delete(imagematch_t *,int n);

/*---------------------------------------------------------------------------*/

/* to record search stats */
typedef struct {
  int nbin;
  float *bins;  
  int *hamming_dist_bins;
  int hamming_dist_nbit;
} delta_a_stats_t;

/*! @brief lowehough parameters structure
 *
 *  Parameters for Lowe-Hough method. Default values are good for a
 *  320*240 image
 */
typedef struct {

  /* use an anisotropic scaling model in 2 dimensions instead of a
     similarity model during the Hough transform stage? */
  int use_2aff_model;

  /* layout of bin sizes:
     scale (in log) 
     angle 
     translation x 
     translation y 
  */
  double bin_sizes[4];

  /* layout of bin sizes:
     scaling factor in x (log)
     translation in x
     scaling factor in y (log)
     translation in y
  */
  double bin_sizes_2aff[4];

  /* don't handle more than this many bins */
  int max_nbin;

  /* min nb of matches before we try to estimate an affine transform (>=3) */
  int min_match_before_affine;

  /* when searching for agreeing point matches, allow this tolerance
     (in pixels) */
  double pos_error_in_affine;

  /* need at least this number of agreeing matches */ 
  int min_match_after_affine;  

  /* when verifying if points are not too close to each other, use
     this threshold (in pixels) */
  double distinct_tolerance;

  /* the bigger, the more diagnostic messages are written to stdout
     (0=silent) */
  int verbose;
  
  /* a priori weighting of deformations */

  /* enabled? */
  int weight_deformation;

  /* tolerances on parameters */
  double sigma_logscale;  /* scaling */
  double sigma_logar;     /* aspect ratio change */
  double sigma_a1;        /* rotation angle of horizontal axis */
  double sigma_a12;       /* difference of rotation of the vertical axis */
  int weight_portrait;    /* are pictures in portrait format weighted
                             the same as landscape? */
  double min_weight;      /* strip weights below this  */


  /* old weighting with hard thresholds (not used if weight_deformation) */

  /* absolute val of determinant of affine matrix should be in 
     [1/max_det,max_det] */
  double max_det;

  /* aspect ratio (small axis/big axis) of ellipse should be bigger
     than this */ 
  double min_ar;
 
  /* finally estimated affine matrix should not transform (0,0) more
     than this many bins away from the initial bin */
  double max_bin_deviation;
  
  delta_a_stats_t delta_a_stats;

} lowehough_parameters_t;

/*---------------------------------------------------------------------------*/

/*! @brief Intitialize structure with default parameters */
void lowehough_parameters_default(lowehough_parameters_t *params);

/*! @brief Filter point matches with Lowe-Hough geometrical filter. 
 *
 * params may be NULL (default values are used)
 */
void shortlist_filter_lowehough(shortlist_t *,imdesc_t *,imagematch_t *,
                                lowehough_parameters_t *params,
                                void (*peek_fun)(void *arg,double frac),
                                void *peek_arg);


/*---------------------------------------------------------------------------*/
/* Handling of image pairs (for video)                                       */
/*---------------------------------------------------------------------------*/



typedef struct image_pair_t {
  pointmatch_t *pm0; /*!< First point match of this pair in imm */
  double score;      /*!< Sum of scores of pts of this pair during match */
  int n_match;       /*!< Nb of matches */
} two_image_t;

typedef struct {
  float xmin,ymin,xmax,ymax;
} ptset_bbox_t;

/*---------------------------------------------------------------------------*/

int ptset_bbox_count_pts(ptset_bbox_t *bb,imdesc_t *imdesc);

/*---------------------------------------------------------------------------*/

/*! @brief List of image pairs
 *
 *  [(im0_0,im1_0),....,(im0_i,im1_i),....]
 */
typedef struct {
  
  two_image_t *pairs; 
  int n_pair;


  /* im0=catenated descriptors for im0_0....im0_i */ 
  imdesc_t *im0,*im1;  
  
  imagematch_t imm;

  ptset_bbox_t bbox0,bbox1;

} image_pairs_t;

/*---------------------------------------------------------------------------*/
/* im,im1 belong to caller, must remain valid until delete */

/*! @brief Create an empty image pairs list */
image_pairs_t *image_pairs_new();

/*! @brief Delete image pairs list */
void image_pairs_delete(image_pairs_t *ip);

/*! @brief Query and db images should be vwgeo's. Return index of match pair
 *  in pair_scores
 */
int image_pairs_add_match_vw(image_pairs_t *,imdesc_t *im0,imdesc_t *im1,int hamming_thresh);


/*! @brief run filter */
void image_pairs_filter_lowehough(image_pairs_t *,
                                  lowehough_parameters_t *params);





#endif



