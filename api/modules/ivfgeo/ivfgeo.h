/*---------------------------------------------------------------------------*/
/*! @defgroup ivfgeo Inverted file
 *
 *  @brief Defines the type and function when one want to use an inverted
 *  file in memory. The algorithms are described in: 
 *
 *  "Hamming embedding and weak geometric consistency for large scale image search"
 *  H. Jegou, M. Douze, C. Schmid, ECCV 2008
 */
/*---------------------------------------------------------------------------*/

#ifndef __invfilegeo_h
#define __invfilegeo_h

/*---------------------------------------------------------------------------*/

#include <siftgeo/siftgeo.h>

/*---------------------------------------------------------------------------*/
/* number of distinct scales and angles resulting from the previous quantization scheme */

extern const int ivfgeo_nb_scales;
extern const int ivfgeo_nb_angles;

/*---------------------------------------------------------------------------*/

/* define the different types of weighting schemes */
#define TFIDF_NO                0
#define TFIDF_OXFORD            1

/* log(N/N_w),where N is the total number of descriptors and N_w is
 * the number of descriptors assigned to visual word w
 */
#define TFIDF_NISTER            2

/* Number of bits used to represent the quantized values of scales and angles */
#define NB_BITS_VAL_ID          21  /* should be increased for bases of > 2M images */
#define NB_BITS_VAL_SCALE       5
#define NB_BITS_VAL_ANGLE       6
#define NB_BITS_VAL_BINSIGN     64  /* this value is a refinement of the position in a cluster */

/* Default log-scale quantization step */
#define LOG_SCALE_QUANTIZATION_STEP  1.2

/* wgc_type is a bitwise or of the 6 following constants */
/* 0) no geometry: ivfgeo_query_scale_angle does the same as ivfgeo_query */
#define AS_MAX_0                0

/* 3) histograms are smoothed by these sizes in angle and scale */
#define AS_MAX_NEIGH_SIZE(na,ns) ((ns)<<16 | (na)<<8)

/* 4) if != 0 before computing maximum, we multiply histogram bins
   with "Bayesian" coefficients */
#define AS_MAX_PREMULTIPLY      0xf0

/* weight only 0 */
#define AS_MAX_PREMULTIPLY_1    0xc0

/* weight 0 and +-pi/2 */
#define AS_MAX_PREMULTIPLY_3    0x60

/* weight all 4 1/4 turns */
#define AS_MAX_PREMULTIPLY_4    0x40

/*---------------------------------------------------------------------------*/
/*! @{ */
/*---------------------------------------------------------------------------*/

/*! @brief Type used for the distance values. */
typedef float distype_t;

/*! @brief The kind of normalization that is used. */
enum norm_type_e {
  norm_none,                    /*!< don't do normalization */
  norm_type_1,                  /*!< Manhantan distance */
  norm_type_2,                  /*!< Euclidean distance (default) */
  norm_type_1_sqrt,             /*!< root square of Manhantan distance */
};

/*---------------------------------------------------------------------------*/
/* Ivfgeo Structure definition                                               */
/*---------------------------------------------------------------------------*/

#ifndef SWIG
#define PACKTHIS __attribute__ ((__packed__))
#else
/* SWIG does not parse __attribute__ */
#define PACKTHIS 
#endif

/*! @brief The elementary data stored within an inverted file. */
typedef struct ivfgeo_elem_s {
  unsigned int id:NB_BITS_VAL_ID;        /*!< The identifier of the image */
  unsigned int scale:NB_BITS_VAL_SCALE;  /*!< Quantized scale of the region */
  int angle:NB_BITS_VAL_ANGLE;           /*!< Quantized angle */
  binsign_t binsign;                     /*!< Binary signature */
} PACKTHIS ivfgeo_elem_t;

/*---------------------------------------------------------------------------*/
/*! @} */
/*---------------------------------------------------------------------------*/

/* burstiness flags */
#define BURSTNORM_CANCEL_MM	1        /* strip multiple matches */
#define BURSTNORM_INTRA		2        /* apply intra-image burstiness downweighting */
#define BURSTNORM_INTER		0x10     /* idem for inter-image */

#define BURSTNORM_INTRA_FUNC_MASK 0xff00
#define BURSTNORM_INTER_FUNC_MASK 0xff0000

/*---------------------------------------------------------------------------*/
/*! @addtogroup ivfgeo
 *  @{ */
/*---------------------------------------------------------------------------*/

/*! @brief A structure representing the inverted file (loaded in memory). */
typedef struct ivfgeo_s {
  int nbvw;                     /*!< number of visual words */
  int nbvec;                    /*!< number of images stored */
  int nbvec_max;                /*!< maximum number of vectors stored in the inverted file */
  int norm_type;                /*!< modify the default normalization scheme for the frequency vector normalization */

  ivfgeo_elem_t **elems;        /*!< the data corresponding to the frequency vectors */
  int *segsize;                 /*!< the amount (number of elements) of data allocated for a given vw */
  int *nbelems;                 /*!< the number of elements stored for a given vw */

  distype_t *norms;             /*!< norms associated with the images (allocated size:nbvec_max) */
  int *labels;                  /*!< some labels that are associated with the images (optional) */

  int tfidf;                    /*!< define the kind of tf-idf scheme which is applied */
  distype_t *tfidf_weights;     /*!< tf-idf weights (allocated size: nbvw) */


  int he_dist_w;                /*!< !=0: weight bins is a function of the Hamming distance 
                                      =1: weight based on Shannon Information content 
                                      =2: Gaussian weighting */

  int he_thres;                 /*!< Hamming comparison of binary signatures threshold */
  double * he_dis_weights;      /*!< weights applied for the different values of the Hamming distance */

  int wgc_type;                 /*!< method to extract maximum from the angle / scale histogram(s) */
  distype_t *wgc_weights;       /*!< a-priori weighting of angle and scale differences */
  
  double scale_w;               /*!< weighting applied to favor bigger scales */
  double * scale_weights;       /*!< weights applied for the different scales */

  int burstiness;               /*!< !=0: burstiness receives a particular treatement
                                 *   burstiness is a binary or the BURSTNORM_* flags */

  int mm_fd;                    /*!< memory-mapped file */
  size_t mm_length;             /*!< length of mapped data */
  void *mm_base;                /*!< base pointer */
} ivfgeo_t;


/*---------------------------------------------------------------------------*/
/* Matching elements                                                         */
/*---------------------------------------------------------------------------*/

#define NB_BITS_VAL_DQID 23

/*! @brief The elementary data of the matching structure. */
typedef struct ivfgeo_match_elem_s {
  float score;                           /*!< The match strength */
  unsigned int ds: NB_BITS_VAL_SCALE+1;  /*!< Quantized scale of the region */
  unsigned int da: NB_BITS_VAL_ANGLE;    /*!< Quantized angle */
  unsigned int id : NB_BITS_VAL_ID;      /*!< The identifier of the image   */
  unsigned int dqid : NB_BITS_VAL_DQID;  /*!< The identifier of the vector within the query image */
} PACKTHIS ivfgeo_match_elem_t;

/*! @brief Data of the matching structure. */
typedef struct ivfgeo_match_s {
  long n;                           /*!< number of matches found */
  long nmax;                        /*!< allocated number of matches */
  ivfgeo_match_elem_t * matches;    /*!< the matches */
} ivfgeo_match_t;

/*---------------------------------------------------------------------------*/
/* Ivfgeo interface                                                          */
/*---------------------------------------------------------------------------*/

/*! @brief Create a void inverted file in memory.
 *
 *  The maxn value indicated the maximum of frequency vectors stored in the inverted file
 *  without re-allocation.\n
 *  Note that if an ivfgeo_add is performed, maxn is automatically 
 *  increased to accept new frequency vectors (geometric re-allocation).    
 */
ivfgeo_t *ivfgeo_new (int nbvw, int maxn);

void ivfgeo_delete (ivfgeo_t * ivf);

/*! @brief Add a -sparse- vector an existing inverted file.
 *
 *  The label is the number assigned to this vector when querying.\n 
 *  By default, if set to -1, it is set to the same value as the internal identifier.\n
 *  The function returns an integer which is the internal identifier of the image (-1 if failed).\n
 *  Note that this internal identifier may be different from the label stored in the set
 */
int ivfgeo_add (ivfgeo_t * ivf, const vwgeoset_t * vw, int label);

/*! @brief Sanity checks on bounds for debugging */
void ivfgeo_check(ivfgeo_t * ivf);

/*! @brief Useful for resume */
int ivfgeo_max_label (const ivfgeo_t * ivf);

/*! @brief Query vw.
 *
 *  Returns NULL if vw is empty.\n
 *  elif n<0: \n
 *    returns NULL and sets *retdis to the array of scores indexed by label \n
 *  else: \n
 *    returns imno and dis=*retdis such that, for 0<=i<min(n,ivf->nbvec) \n
 *    imno[i] is the i^th closest vector to vw \n
 *    dis[i] is the score of image imno[i] \n
 */
int *ivfgeo_query (const ivfgeo_t * ivf,
                   const vwgeoset_t * vw, 
		   distype_t ** retdis, int n);

/*! @brief ivfgeo_query() with a peek function to monitor the progress of the computation */
int *ivfgeo_query_peek (const ivfgeo_t * ivf,
                        const vwgeoset_t * vw,
                        distype_t ** retdis, int n,
                        void (*peek_fun) (void *arg, double frac),
                        void *peek_arg);

/*! @brief Open an inverted file and load it into memory.
 *
 *  fmt=0 standard \n
 *  fmt=1 memory image (ivfgeo_add not possible) \n
 *  fmt=2 load memory image (ivfgeo_add or ivfgeo_merge possible) \n
 *  fmt=3 memory image, scan all to make sure it is in RAM
 */
ivfgeo_t *ivfgeo_read (const char *filename, int fmt);

/*! @brief Merge ivf2 into ivf1.
 *
 *  The second file is deallocated on exit. \n
 *  Labels in ivf2 are increased by label2ofs
 */
void ivfgeo_merge (ivfgeo_t * ivf1, ivfgeo_t * ivf2, int label2ofs);

/*! remove entries with label>=maxlabel */
void ivfgeo_crop (ivfgeo_t * ivf, int maxlabel);

/*! keep maximum one entry per image for each visual word */
void ivfgeo_filter_duplicate_vw (ivfgeo_t * ivf);

/*! mask binary signatures */
void ivfgeo_mask_binsigns (ivfgeo_t * ivf,unsigned long long mask);


/*! @brief Write an inverted file on disk.
 *
 *  fmt=0 standard \n
 *  fmt=1 store flat dump of memory structures
 */
void ivfgeo_write (const char *filename, const ivfgeo_t * ivf, int fmt);

/*! @brief Also computes weighting tables */
void ivfgeo_set_wgc_type (ivfgeo_t * ivf, int);

/*! @brief Compute the norms associated with the elements of the inverted file.
 *
 *  This function should not be used, except is one decide to change the norm
 *  used when querying
 */
void ivfgeo_compute_norms (ivfgeo_t * ivf);

/*! @brief Compute the tf-idf weights according to the criterion ivf->tfidf */
void ivfgeo_compute_tfidf (ivfgeo_t * ivf);

/*! @brief Weighting of scales: bigger scales get higher weights.
 *
 *  The weighting depends on the variable ivf->scale_w
 */
void ivfgeo_compute_scale_weights(ivfgeo_t * ivf);

/*! @brief Weighting of hamming distances: smaller distances get higher weights */
void ivfgeo_compute_he_weights(ivfgeo_t * ivf);

/*! @brief Display the contents of an inverted file  */
void ivfgeo_display (const ivfgeo_t * ivf);

/*! @brief Display a summary of the contents of an inverted file */
void ivfgeo_summary (const ivfgeo_t * ivf);

/*---------------------------------------------------------------------------*/
/* Matching Structure Interface                                              */
/*---------------------------------------------------------------------------*/

/*! @brief Create the structure that will receive the matches */
ivfgeo_match_t * ivfgeo_match_new (int maxn);

/*! @brief Remove the structure */
void ivfgeo_match_delete (ivfgeo_match_t * matches);

/*! @brief Display the content of the structure */
void ivfgeo_match_display (const ivfgeo_match_t * matches);

/*---------------------------------------------------------------------------*/
/* Statistical Information                                                   */
/*---------------------------------------------------------------------------*/

/*! @brief Count the total number of elements (descriptors) in the inverted file */
int ivfgeo_count_nbelems (const ivfgeo_t * ivf);

/*! @brief Compute the "unbalanced factor" of the inverted file */
double ivfgeo_unbalanced_factor (const ivfgeo_t * ivf);

/*! @brief Display the contents of an inverted file */
int * ivfgeo_compute_burstiness_hist (const ivfgeo_t * ivf, int maxcount);

/*! @brief Return an histogram of the difference between images id */
int * ivfgeo_histo_diff_im_id (const ivfgeo_t * ivf);


/*---------------------------------------------------------------------------*/
/* Exported utility functions                                                */
/*---------------------------------------------------------------------------*/

int binsign_hamming (unsigned long long bs1, unsigned long long bs2);

/*! compute all Hamming distances between 2 sets of vectors 
 *
 *  a is na * sizeof_code
 *  b is nb * sizeof_code
 *  dist2[ i + na * j ] = || a(i,:)-b(j,:) || \n
 */
void binsign_cross_distances(int sizeof_code,int na,int nb,
                             const char *a,
                             const char *b,
                             float *dists);



int quantize_scale (geom_t * geom);

int quantize_angle (geom_t * geom);

void ivfgeo_compute_scale_weights_tab(int nb_scales,double scale_w,double **scale_weights_out);

void ivfgeo_compute_wgc_tab (int nb_angles,int nb_scales,int wgc_type,
                             float **wgc_weights_out);

distype_t synth_wgc_tab (distype_t * tab, int synth_type,
                         distype_t * weights);

/*---------------------------------------------------------------------------*/
/* Query stats (not MT-safe)                                                 */
/*---------------------------------------------------------------------------*/

/*! @brief Query stats (not MT-safe). */
typedef struct query_stats_s {
  int n_images;
  int n_points; 

  long n_visited_elements;
  long n_he_success;
  
  long n_mm_cancel;

  float *as_hist;

  int elem_size;
  int match_elem_size;
  
} query_stats_t;

extern query_stats_t query_stats;

/*---------------------------------------------------------------------------*/
/*! @} */
/*---------------------------------------------------------------------------*/

#endif

/*---------------------------------------------------------------------------*/

