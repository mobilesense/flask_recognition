#ifndef __siftgeo_agg_h
#define __siftgeo_agg_h

#include "siftgeo.h"


/* Define a structure that contains a set of sparse functions used in BOF aggregation */
typedef struct {
  int n;                        /* number of atoms */
  int d;                        /* dimension of the sparse vectors which are stored */
  int *nz;                      /* number of non-zeros positions for each atom */
  int **idx;                    /* positions of the non-zeros values */
  float **val;                  /* values associated with the idx positions */
} atoms_t;


/*---------------------------------------------------------------------------*/
/* Elementary operations on atoms                                            */
/*---------------------------------------------------------------------------*/

atoms_t *atoms_new (int n, int d);
void atoms_free (atoms_t * atoms);

/* Add an atom to the structure (in position no) */
void atoms_add_spvec (atoms_t * atoms, int no, int *idx, float *val, int n);

atoms_t *atoms_read (const char *filename);
void atoms_write (const char *filename, const atoms_t * atoms);

void atoms_display (const atoms_t * atoms, int verboselevel);


float *atoms_mul_fvec (const atoms_t * atoms, float * v);

/* keep only vectors between n0 and n1 */
void atoms_crop(atoms_t *atoms,int n0,int n1);

void atoms_sort_indices(atoms_t * atoms);

/*---------------------------------------------------------------------------*/
/* Miscellaneous                                                             */
/*---------------------------------------------------------------------------*/

/*! @brief Compute the distances between a descriptor and a set of 
  descriptors. The output parameter dis should be allocated.       */
void point_to_points_dis (const point_t * set, const point_t * des, int n,
                          float *dis);


/*! convert a set of point into a bag-of-features vector (not normalized) */
int vwsgeo_to_bof (vwgeoset_t * vgs, int **idx_out, float **v_out);

/*
void atoms_project_pointset(const atoms_t * atoms,const pointset_t *ps,
                            int nkeep,int *atom_no,float *dotprod);
*/

/*---------------------------------------------------------------------------*/
/* Construction of special atoms                                             */
/*---------------------------------------------------------------------------*/


/*! Atoms consisting of "and" combinations of Hadamard functions 
 *
 * d must be a power of 2. The resulting atoms structure will have n==d */
atoms_t *atoms_new_hadamard_inter (int d, int l, int *h_indices_out);

/*! Random atoms with supports of size nsup that do intersect at most once */
atoms_t *atoms_new_random_by_support (int n, int d, int nsup);

atoms_t *atoms_new_random_by_sparsity_rate (int n, int d, float rate);

atoms_t * atoms_new_aggregate_components (int n, int d, int nz);


float * atoms_signature (const atoms_t * atoms, vwgeoset_t * vgs, const float * tfidf);

float * atoms_signature_fvec (const atoms_t * atoms, const float*bof, const float * tfidf);

int * hesp_signature (const atoms_t * atoms, vwgeoset_t * vgs, const float * tfidf);


#endif
