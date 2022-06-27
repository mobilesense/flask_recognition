/*---------------------------------------------------------------------------*/
/*! @defgroup siftgeo Siftgeo
 *
 *  @brief Defines the type and function when one want to use a siftgeo
 *  structure
 *
 *  @{
 */
/*---------------------------------------------------------------------------*/

#ifndef SIFTGEO_H_INCLUDED
#define SIFTGEO_H_INCLUDED

/*---------------------------------------------------------------------------*/

#include <stdio.h>

/*---------------------------------------------------------------------------*/
/* Generic SIFTGEO & VWGEO object                                            */
/*---------------------------------------------------------------------------*/

/*! @brief Type for binary signatures. */
typedef long long unsigned int binsign_t;

/*! @brief The structure containing the interest point values */
typedef struct {
  float x, y, scale, angle, mi11, mi12, mi21, mi22, cornerness;
} geom_t;

/*! @brief Interest point structure */
typedef struct {
  geom_t geom;			/*!< Interest point value */
  union {
    int  dim;			/*!< Dimension of the descriptor (for siftgeo) */
    int  vw;			/*!< visual word (for vwgeo) */
  };
  unsigned char *desc;		/*!< Descriptor data */

  binsign_t binsign;		/*!< Binary signature */
} point_t;

/*! @brief A structure to store a set of interest points */
typedef struct pointset_s {
  int  n;
  point_t *pts;
} pointset_t;

/*! @brief A structure to store a visual word */
typedef pointset_t vwgeoset_t;



/*---------------------------------------------------------------------------*/
/* General purpose functions                                                 */
/*---------------------------------------------------------------------------*/

/*! @brief All fields set to 0 */
pointset_t *pointset_alloc (int npt);

pointset_t *pointset_new ();

/*! all other fields are 0 */
pointset_t *pointset_from_vw (int *vw, int n);

int *vw_from_pointset(const pointset_t *); 

/*! @brief Delete pointset */
void pointset_delete (pointset_t * vgs);

/*! @brief Deallocates a set of point_t */
void delete_points (point_t * corners, int nc);


/*! @brief Compute the norm of a given set of vwgeo */
double vwgeoset_norm1 (const vwgeoset_t * vwg);

/*! @brief Should be ordered by vw on input */
double vwgeoset_norm2 (const vwgeoset_t * vwg);


/*! @brief Make new set of vwgeo's with geom_t from the descriptors and vw
 * from vws
 *
 * k=# of vw per point_t
 */
point_t *siftgeo_to_vwgeo (point_t * des, int n, int *vws, int k);


/*! @brief Appends points from ps2 to ps.
 *
 *  On output, ps2.n=0
 */
void pointset_append (pointset_t * ps, pointset_t * ps2);

/*! @brief Transform points with affine transform
 *
 * [aff[0] aff[1]; aff[3] aff[4]]*[x;y]+[aff[2]; aff[5]]\n\n
 *
 * and multiply scale of points with scalef
 */
void pointset_affine_transform (pointset_t * ps, double aff[6],
				double scalef);

/*! @brief Duplicate a pointset */
pointset_t *pointset_dup (pointset_t * ps);



/*---------------------------------------------------------------------------*/
/* I/O                                                                       */
/*---------------------------------------------------------------------------*/
/* point formats: \n 
 * 0 = siftgeo    \n
 * 1 = vwgeo      \n
 * 2 = vwsgeo = vwgeo + binary signature \n
 * 3 = siftbin, the binary format used by Lowe's description code \n
 */

/*! @brief Reads one corner form a siftgeo or vwgeo file
 *
 * returns:       \n
 * 1 for success  \n
 * 0 for eof      \n
 * <0 for failure \n
 */
int  read_point_t (FILE * f, point_t * corner, int fmt);

/*! @brief Writes one corner form a siftgeo or vwgeo file */
int  write_point_t (FILE * f, point_t * corner, int fmt);

/*! @brief Count points from a seekable file */
int  count_points (const char *fname, int fmt);

/*! @brief Count points from a seekable file.
 * if fmt is 0 or 3, the dimension must be determined. If you know it,
 * set *d_io, else set *d_io to -1. The function will read it. If you
 * don't care, just pass NULL.
*/
int pointset_file_size_and_dim (const char *fname, int fmt, int *d_io);


/*! @brief Fills *corners_out.
 *
 * Returns nb of corners or -1 for failure
 */
int  read_points (const char *fname, point_t ** corners_out, int fmt);

/*! @brief Adds corners to *corners_io that already contains nc corners.
 *
 * Returns new value of nc or <0 for failure
 */
int  read_points_add (const char *fname,
		      point_t ** corners_io, int nc, int fmt);

/*! @brief Add point i iff mask[i]!=0 */
int  read_points_add_with_mask (const char *fname,
				point_t ** corners_io,
				int nc, int fmt, int *mask);

/*! @brief Read from already opened FILE* */
int  read_points_file (FILE * f, point_t ** corners_out, int fmt);

/*! @brief Read from string (same layout as in a file) */
int  read_points_string (const char *string, int string_length,
			 point_t ** corners_out, int fmt);

/*! @brief Writes corners to file
 *
 *  return 0 for success
 */
int  write_points (const char *fname, point_t * corners, int nc, int fmt);

/*! @brief Writes corners to string
 *
 *  return 0 for success
 */
int  write_points_string (char *string, int string_length, point_t * corners,
			  int nc, int fmt);

/*! @brief Writes corners to FILE*
 *
 *  return 0 for success
 */
int  write_points_file (FILE *, point_t * corners, int nc, int fmt);


/*! @brief Read a pointset from a file */
pointset_t *pointset_read (char *fname, int fmt);

/*! @brief Read at most n points from a file  */
pointset_t *pointset_read_file_max (FILE *f, int n, int fmt);


/*! @brief Read a pointset list from a file and sort it according to cornerness */
pointset_t *pointset_read_cornerness (char *fname, double min_cornerness,
				      int fmt);

/*! @brief Write pointset list to a file */
void pointset_write (char *fname, pointset_t *, int vwgeo);


/*! @brief Display a set of descriptors */
void display_points (const point_t * des, int n);



/*---------------------------------------------------------------------------*/
/* Sort and filtering functions                                              */
/*---------------------------------------------------------------------------*/

/*! @brief This function sort all the vwgeo set such that the indexes 
 *  of the visual word are in increasing order
 */
void vwgeoset_sort (vwgeoset_t * vgs);


/*! @brief Filter the vwgeo set such that only the n vwgeo with highest
 *  cornerness are kept */
void vwgeoset_filter_n_cornerness_max (vwgeoset_t * vgs, int n);


/*! @brief Same as siftgeo_to_vwgeo, but keep only the points assigned to centroids
 *  which are at a similar distance (disratio*dismin) from the nearest centroid
 */
void vwgeoset_filter_ma (vwgeoset_t * vgs, int k, float *vw_centroids_dis,
			 double disratio);

/*! @brief Same as vwgeoset_filter_ma, but report how many NN have been kept for each group of k
 */
void vwgeoset_filter_ma_nkeep (vwgeoset_t * vgs, int k, float *vw_centroids_dis,
                               double disratio,int *nkeep);

/*! @brief Same as vwgeoset_filter_ma_nkeep, but use a hard threshold on distances
 */
void vwgeoset_filter_thresh_nkeep (vwgeoset_t * vgs, int k, float *vw_centroids_dis,
                                   double thresh,int *nkeep);


/*! @brief remove the multiple occurence of the same visual words and keep only one,
  producing a binary bag-of-features representations                         */
void vwgeoset_filter_duplicate_vw (vwgeoset_t * vgs);


/*! @brief Sort by cornerness */
void pointset_sort_by_cornerness (pointset_t * vgs);


/*! @brief Sort the elements of vgs using the order indicated in order */
void pointset_sort_by_permutation (pointset_t * vgs, const int * permutation);


/*! @brief Keep points inside/outside a bounding box */
void pointset_crop (pointset_t * ps, float xmin, float ymin, float xmax,
		    float ymax, int keep_outside);

/*! keep points inside/outside a polygon
 * 
 *  Keep points such that \n
 *  [ coeffs[3*i] coeffs[3*i+1] coeffs[3*i+2] ] * [x y 1]' > 0 \n
 *  for i=0..n-1
 */
void pointset_crop_polygon (pointset_t * ps, double *coeffs, int n);


/*! @brief crop n points (keep the first n ones) */
void pointset_crop_n (pointset_t * ps, int n);


/*! @brief Randomly keep some points with probability */
void pointset_filter_random (pointset_t * vgs, float pkeep);


/*! @brief keep points with cornerness>=cornerness_min */
void pointset_filter_cornerness (pointset_t *ps, float cornerness_min);

/*! @brief mask binary signature */
void vwgeoset_mask_binsign (pointset_t * ps, binsign_t mask);


/*---------------------------------------------------------------------------*/
/*! @} */
/*---------------------------------------------------------------------------*/

#endif

/*---------------------------------------------------------------------------*/
