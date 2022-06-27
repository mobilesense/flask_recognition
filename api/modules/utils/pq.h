/* product quantizer and corresponding search methods */

#ifndef __pq_h
#define __pq_h

#include <stdio.h>

typedef float real;

#define PQMODE_STD   0
#define PQMODE_8x8   1

/* the structure that store the prod quantizer parameters */
typedef struct pq_s {
  int d;                  /* the dimension of the vector to represent */
  int nsq;                /* nb of subquantizers used in parallel */
  int nsq_bits;           /* nb of bits allocated to each subquantizer */
  int ksq;                /* nb centroid per subquantizer */
  int mode;               /* optimization mode */
  real ** centroids;      /* centroids[i] gives subquantizer i centroid*s* */
  real ** centroids_dis;  /* the distance between the centroids */
  real ** sberr;          /* sberr[i][j]:  quantiz error for (sb i, cell j) */
  int * extract_map;      /* permutation to define how a vector is mapped to subvectors */ 
} pq_t;


/* Generic type used for pqcode_t */
typedef unsigned char pqcode_t; 



/* [*** pqcodes_t functions ***] */

/* alloc/free array of pq codes */
pqcode_t * pqcode_new (const pq_t * pq, int n);

pqcode_t * pqcode_realloc (const pq_t * pq, pqcode_t *code, int n);

long pqcode_sizeof (const pq_t * pq);

/* return pointer on i^th code of the array */
pqcode_t *pqcode_offset (const pq_t * pq, const pqcode_t *code, int i);

void pqcode_free (pqcode_t * codes);



/*--- Input/Output ---*/

void pqcode_fwrite (FILE * f, const pq_t * pq, const pqcode_t * codes, int n);
void pqcode_write (const char * filename, const pq_t * pq, const pqcode_t * codes, int n);

pqcode_t * pqcode_new_fread (FILE * f, const pq_t * pq, int n);
pqcode_t * pqcode_new_read (const char * filename, const pq_t * pq);

void pqcode_display (const pq_t * pq, const pqcode_t * codes, int n);


/* [*** pqcodes_t functions ***] */

/*--- General functions ---*/

/* learn the structure from a vector set */
pq_t * pq_new (int d, int nsq, int nsq_bits);

/* free memory */
void pq_free (pq_t * pq);

/* learn the quantizer on a set of n vectors v */ 
pq_t * pq_learn (int d, int nsq, int nsq_bits, 
		 const real * v, int n, int nb_iter_max);


/* learn with a given extract map */
pq_t *pq_learn_with_extract_map (int d, int nsq, int nsq_bits,
                                 const real * v, const int *extract_map,
                                 int n, int nb_iter_max); 

/*--- Input/Output ---*/

/* print the product quantizer */
void pq_display (const pq_t * pq, int verbose);

/* write down the parameters */
void pq_fwrite (FILE * fe, const pq_t * pq);
void pq_write (const char * fname, const pq_t * pq);

/* read the parameters from a file */
pq_t * pq_fread (FILE * f);
pq_t * pq_read (const char * fname);

/*--- encode/decode product quantizer representations ---*/

/* encode the input vectors */
void pq_encode (const pq_t * pq, const float * v, int n, pqcode_t * codes_out);


/* decode the input codes */
void pq_decode (const pq_t * pq, const pqcode_t * codes, int n, float * v_out);

/* compute the residual vector associated with a vector */
void pq_residual (const pq_t * pq, const real * v, int n, real * resi);

/* the same but with known pre-computed codes */
void pq_residual_from_codes (const pq_t * pq, const real * v, const pqcode_t * codes, 
			     int n, real * resi);


/*--- Expected distance and search: the symmetric case ---*/

/* compute the expected square distance between an encoded query code1 and n encoded queries codes2. 
   The distances are stored in dis_out (externaly allocated) 
*/
void pq_L2sq_cor (const pq_t * pq, const pqcode_t * codes1, 
		  const pqcode_t * codes2, int n, real * dis_out);

/* The same but with no correcting term for the quantizer (better for NN search) */
void pq_L2sq (const pq_t * pq, const pqcode_t * codes1, 
	      const pqcode_t * codes2, int n, real * dis_out);

/* find the k-nearest neighbors for n1 coded points within a set of n2 vectors.
   Also output the distance, but only if nn_dis is allocaed (!= NULL) */
void pq_nns_cor (const pq_t * pq, const pqcode_t * codes1, 
		 const pqcode_t * codes2, int n1, int n2, int k, 
		 int * nns_out, real * nns_dis_out);

void pq_nns (const pq_t * pq, const pqcode_t * codes1, 
	     const pqcode_t * codes2, int n1, int n2, int k, 
	     int * nns_out, real * nns_dis_out);


/*--- Expected distance and search: asymetric case (query not quantized) */


/* compute the expected distance quantized query and points */
void pq_asym_L2sq_cor (const pq_t * pq, const real * v, 
		       const pqcode_t * codes, int n, real * dis_out);

/* biaised version (preferred) */
void pq_asym_L2sq (const pq_t * pq, const real * vquery, 
		   const pqcode_t * codesb, int n, real * dis_out);

/* find the k-nearest neighbors for a set of queries within a set of n2 vectors, 
   using the asymetric method: query is a float vector, and the base are pq-codes */
void pq_asym_nns_cor (const pq_t * pq, 
		      const real * v, const pqcode_t * codes2, 
		      int n1, int n2, int k, 
		      int * nns, real * nns_dis);

void pq_asym_nns (const pq_t * pq, 
		  const real * v, const pqcode_t * codesb, 
		  int n1, int n2, int k, 
		  int * nns_out, real * nns_dis_out);

/* here, the query is the code and the base a set of float vectors */
void pq_asym_fc_nns (const pq_t * pq, const pqcode_t * codesq, 
		     const real * vbase, 
		     int n1, int n2, int k, 
		     int * nns, real * nns_dis);


/*--- k-means clustering using codes ---*/

float * pq_kmeans_clustering (const pq_t * pq, const pqcode_t * codes, 
			      int n, int k, int nb_iter_max, 
			      int * clust_assign_out);

/* threaded NN search */

void pq_nns_thread (const pq_t * pq, const pqcode_t * codes1, 
                    const pqcode_t * codes2, int n1, int n2, int k, int nthread,
                    int * nns_out, real * nns_dis_out);


void pq_asym_nns_thread (const pq_t * pq, 
                         const real * v, const pqcode_t * codesb, 
                         int n1, int n2, int k, int nthread,
                         int * nns_out, real * nns_dis_out);



#endif
