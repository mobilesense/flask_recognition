#ifndef __hkm_h
#define __hkm_h


/*--------------------------------------------------------------
 * hierarchical clustering
 --------------------------------------------------------------*/

/* the structure used for the quantization */
typedef struct hkm_s {
  int nlevel;            /* number of levels */
  int bf;                /* the branching factor */
  int k;                 /* the number of leaves (bf^nlevel) */
  int d;                 /* dimension of the input vectors */
  float ** centroids;    /* centroids for all levels */
  float ** balfact;      /* multiplicative factor associated with the centroid assigment */

} hkm_t;


/* create/delete a hierarchical quantizer structure.
   nlevel is the number of levels in the tree and bf the branching factor */
hkm_t * hkm_learn (int n, int d, int nlevel, int bf, 
		   float *points, int nb_iter_max, 
		   int balmode, 
		   int ** clust_assign_out);
		 
void hkm_delete (hkm_t * hkm);

/* Quantization usign the hierarchical clustering */
int * hkm_quantize (const hkm_t * hkm, const float * v, int npt);

/* I/O function for hkm */
void hkm_write (const char * filename, const hkm_t * hkm);
hkm_t * hkm_read (const char * filename);

/* retrieve the centroids from a particular level */
float * hkm_get_centroids (const hkm_t * hkm, int l, int no);


/* compute the multiplicative factor such that the clusters are almost balanced */
void compute_balance_factors (const float * centroids, const float * v, 
			      int k, int n, int d, float * balfact);


#endif
