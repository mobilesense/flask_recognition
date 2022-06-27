/* inverted file with geometric information and product quantizer-indexed data  */

#ifndef __ivfpq_h
#define __ivfpq_h

#include "pq.h"


/*---------------- Ivfpq Structure definition ----------------*/

/* Structure representing the inverted file (loaded in memory)) */

#define IVFPQ_FLAG_SYM            0
#define IVFPQ_FLAG_ASYM           1 

#define IVFPQ_FLAG_ABSOLUTE       0
#define IVFPQ_FLAG_RELATIVE       2 

#define IVFPQ_FLAG_CENTROID_REAL  0 
#define IVFPQ_FLAG_CENTROID_PQ    4


typedef struct ivfpq_t {
 
  pq_t *pq;                     /* product quantizer associated with the codes */

  int nbvw;                     /* number of visual words */
  
  int *nbelems;                 /* the number of elements stored for a given vw */
  int *segsize;                 /* the amount (number of elements) of data allocated for a given vw */

  int **labels;                 /* indexed vectors */
  pqcode_t **codes;             /* codes corresponding to the indexed vectors */

  int flags;

  float *centroids;             /* centroids used to build the ivf */

  /* version with product quantizer for the centroids */
  pqcode_t * centroids_codes;   /* product-quantized version of the centroids */
  pq_t * pq_assign;  
} ivfpq_t;



/*-------------------- Low-level interface (VW's are pre-processed) --------------------*/

/* create a void inverted file in memory. The pq will be deleted when
   the ivfpq_t is dealloc'ed */
ivfpq_t *ivfpq_new (int nbvw);

/* suppress the inverted file */
void ivfpq_delete (ivfpq_t * ivf);

/* Add n vectors to the inverted file. val is a n * pq->d array */
void ivfpq_add_vw (ivfpq_t * ivf, const int *vw, const int *labels, const float * val, int n);

/* low-level query

   Fills arrays labels and vals with the at most k nearest neighbours
   of val from the union of cells vw[0:n]. 

   if multiple_val, val is a n*d matrix with a different query vector for each cell. 
 */
int ivfpq_query_vw (const ivfpq_t * ivf, const int *vw, int n, 
		    const float * val, int k, int multiple_val,
                    int *labels,float *dists);

/* Count the total number of elements (descriptors) in the inverted file */
int ivfpq_count_nbelems (const ivfpq_t * ivf);

double ivfpq_unbalanced_factor(const ivfpq_t *ivf);

/* dup & merge to implement parallel adds.*/
ivfpq_t *ivfpq_dup(ivfpq_t *ivf);

/* ivf2 deallocated on output*/
void ivfpq_merge(ivfpq_t *ivf,ivfpq_t *ivf2);


void ivfpq_display(const ivfpq_t *ivf); 


/*-------------------- with centroids --------------------*/


/* train the pq & build the ivf */
ivfpq_t *ivfpq_new_learn(int nbvw, const float *centroids,
                         int d, int nsq, int nsq_bits,
                         const real * v, int n, int nb_iter_max,
                         int flags, int nt);


ivfpq_t *ivfpq_new_learn_with_extract_map (int nbvw, const float *centroids,
                                           int d, int nsq, int nsq_bits,
                                           const real * v, const int* extract_map, int n, int nb_iter_max,
                                           int flags, int nt);


/* use only for parallel add's */
ivfpq_t *ivfpq_dup(ivfpq_t *ivf);
void ivfpq_merge(ivfpq_t *ivf,ivfpq_t *ivf2);

/* add n vectors */
void ivfpq_add (ivfpq_t * ivf, const int *labels, const float * val, int n, int nt);


void ivfpq_add_with_vw (ivfpq_t * ivf, const int *labels, const float * val, int n, const int *vw, int nt);

int *ivfpq_add_and_get_vw (ivfpq_t * ivf, const int *labels, const float * val, int n, int nt);


/* 
query n vectors return k results per vector

for query i, the results and associated distances are in 

 labels[k*i : (k+1)*i]  
 dists[k*i : (k+1)*i] 

 if there are not enough results, labels for missing results are set to -1

*/

void ivfpq_query (const ivfpq_t * ivf, const float * val, int n, 
                  int ma,double ma_disratio, int k, 
                  int nt, 
                  int *labels, float *dists); 

/* 
in addition to the previous, this returns 

vw[i*ma : (i+1)*ma]: visual words corresponding to the query (-1 for absent)

for query i, visual word j:

let 
  i0=label_map[i*ma+j-1]
  i1=label_map[i*ma+j]

then 

  labels[k*i+i0 : k*i+i1]: labels corresponding to visual word j.

*/
void ivfpq_query_and_get_vw (const ivfpq_t * ivf, const float *val, int n,
                             int ma, double ma_disratio,
                             int k, int nt, 
                             int *labels, float *dists, 
                             int *vw,int *label_map);

#if 0

/* combine results from several ivfs. The dimension of val is n*d, where 
   d=sum(ivfs[i]->pq->d,i=0..n-1)
*/
void ivfpqs_query(const ivfpq_t ** ivfs, int ndb,
                  const float *val, int n,
                  int ma,double ma_disratio, 
                  int k,  
                  int *labels, float *dists); 


/* compute squared distance quantiles for points belonging to the same centroid */

void ivfpq_dist_quantiles(const float *centroids,
                          const float *val, int n, 
                          int nq,
                          float *quantiles);

#endif
                          

/*-------------------- I/O --------------------*/

/* pq is always written and loaded along with the ivf */
void ivfpq_fwrite(FILE *f, const ivfpq_t *ivf);


ivfpq_t *ivfpq_fread(FILE *f, const float *centroids);


typedef struct {
  long nvisited;
} ivfpq_query_stats_t; 

extern ivfpq_query_stats_t ivfpq_query_stats;


#endif
