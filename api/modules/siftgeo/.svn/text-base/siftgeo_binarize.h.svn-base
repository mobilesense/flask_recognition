/*---------------------------------------------------------------------------*/

#ifndef __siftgeo_binarize_h
#define __siftgeo_binarize_h

/*---------------------------------------------------------------------------*/

#include <siftgeo/siftgeo.h>

/*---------------------------------------------------------------------------*/
/*! @addtogroup siftgeo
 *  @{
 */
/*---------------------------------------------------------------------------*/

/*! @brief Parameters used to compute binary signatures */
typedef struct {
  int d;              /*!< dimension */
  int nproj;          /*!< nb of projections */
  int nbvw;           /*!< nb of visual words for which we have medians */
  float *p;           /*!< projection matrix, dimension nproj*d */
  float *medians;     /*!< thresholds to compare against, dimension nbvw*nproj */
} siftgeo_binarize_t;

/*---------------------------------------------------------------------------*/

/*! @brief Create a new binarization structure */
siftgeo_binarize_t * siftgeo_binarize_new (int n, int d, int nbvw, int nproj, 
					   float * points, 
					   int * clust_assign);

/*! @brief Generate a projection matrix using the principal subspace of
 * difference between the vectors and their centroids */
siftgeo_binarize_t * siftgeo_binarize_pca_new (int n, int d, int nbvw, 
					       int nproj, int whiten,
					       float * points, 
					       int * clust_assign, 
					       float * centroids);


/*! @brief Fillin the binsign in the vwgeos */
void siftgeo_binarize (siftgeo_binarize_t *,pointset_t *siftgeos,pointset_t *vwgeos);

/*! @brief Fillin the binsign in the vwgeos, but the points are in a d * vwgeos->n
 *  float array (to speed up coordinate transforms) 
 *
 *  k = # of instances of point for a given vwgeo
 */
void siftgeo_binarize_ffq (siftgeo_binarize_t *,const float *points,pointset_t *vwgeos,int k);


/*! limit nb of projections */
void siftgeo_binarize_crop_nproj(siftgeo_binarize_t *,int new_nproj);

/*! Fillin a set of binsigns. 
 *
 * for i=0:n-1
 * 
 * points[i*d:(i+1)*d-1] are the coordinates of point i
 * assign[i*k:(i+1)*k-1] gives the assignements of point i
 * 
 * on output, fill in 
 * 
 * binsigns[i*k*nl:(i+1)*k*nl-1] where nl=ceil(nproj/64)
 *  
 */
void siftgeo_binarize_ffq_table (siftgeo_binarize_t *,
                                 int np,int k,
                                 const float *points,const int *assign,
                                 unsigned long long *binsigns);


/*! @brief Delete structure */
void siftgeo_binarize_delete(siftgeo_binarize_t *);

/*! @brief Fillin the medians for the given visual word using these points 
 * 
 * vw: index of the vw to be filled in         \n
 * points: learning points (size npt * sb->d)  \n
 * npt: nb of learning points
 */
void siftgeo_binarize_fill_medians(siftgeo_binarize_t *sb,int vw,float *points,int npt);

/*! @brief Write description to stdout (verbose) */
void siftgeo_binarize_display(siftgeo_binarize_t *bs);

/*! @brief Read params */
siftgeo_binarize_t * read_sb_common (FILE *f);

/*! @brief Read params */
siftgeo_binarize_t *siftgeo_binarize_read(const char*fname);

/*! @brief Write params to disk */
void write_sb_common (FILE *f, const siftgeo_binarize_t *sb) ;

/*! @brief Write params to disk */
void siftgeo_binarize_write(const char*fname, const siftgeo_binarize_t *sb);

/*! @brief Compute the binary signatures of a set of descriptors
 *
 * p           the matrix used for the projection of descriptors \n
 * nbproj      number of distinct random projection = number of bits in the signature \n
 * siftgeo     the descriptors \n
 * nbdes       the number of descriptors
 */
binsign_t * siftgeo_binarize_full (int nproj,int nbvw,float *p,
                                   point_t * siftgeo, point_t * vwgeo, 
				   int nbdes, float *medians);

/*! @brief Display a signature */
void siftgeo_binarize_binsign_display (binsign_t bs); 

/*! @brief Use the rotation matrix of sb to project the points */
float * he_transform_points(siftgeo_binarize_t *sb, const float *points, int npt);

/*---------------------------------------------------------------------------*/
/*! @} */
/*---------------------------------------------------------------------------*/

#endif

/*---------------------------------------------------------------------------*/

