/*---------------------------------------------------------------------------*/

#ifndef __annvw_h
#define __annvw_h


/*---------------------------------------------------------------------------*/
/* Approximate cluster assignement                                           */
/*---------------------------------------------------------------------------*/

/*! @brief The 2-level structure that receives NN information */
typedef struct ann_vw_t {
  int nbvw_coarse;
  int nbvw;
  int d;

  float * nodes;         /* level 2 (switching/coarse level) */
  const float * leaves;  /* level 1 (base/initial level) */
                         /* edges[i][j] for 0<=i<nbvw_coarse and 0<=j<n_edge[i]
                            is an arc between level-1 and level-2 nodes */
  int * n_edge; 
  int ** edges;

  float ** leave_blocks; /* block computations */

} ann_vw_t;

/*---------------------------------------------------------------------------*/

/*! @brief Load from a couple of separate files */
ann_vw_t * load_ann_vw (const char * fclustcoarse_name, const char * fedges_name, 
			int nbvw, float * centroids);



void free_ann_vw(ann_vw_t *);

void ann_vw_describe(const ann_vw_t *);

void ann_vw_make_blocks(ann_vw_t *);


/*! @brief Assignment */
void quantize_codebook_annvw (const ann_vw_t * ann, int npt,
                              const float * points, int * vw,
                              void (*peek_fun) (void *arg, double frac),
                              void *peek_arg);

/*! @brief Assignement sliced in threads */
void quantize_codebook_annvw_thread (const ann_vw_t * ann, int npt,
                                     const float * points, int * vw,
                                     int n_thread,
                                     void (*peek_fun) (void *arg, double frac),
                                     void *peek_arg);


/*! @brief Multiple (k) assignement.
 *
 * Returns the set of distances to the centroids.
 */
float* quantize_codebook_annvw_multiple (const ann_vw_t * ann, int npt, float * points,
                                         int k, int * vw);

float* quantize_codebook_annvw_multiple_thread (const ann_vw_t * ann, int npt, float * points,
                                                int k, int * vw,int n_thread);


/* full call n_scaleprod_out returns nb of scalar products */
void quantize_codebook_annvw_full_thread (const ann_vw_t * ann, int npt,
                                          const float *points, int knn, int knn_coarse, 
                                          int *vw, float *vwdis, 
                                          int *n_scalprod_out,
                                          int n_thread,
                                          void (*peek_fun) (void *arg, double frac),
                                          void *peek_arg);


/*****************************************************************************
 * learning part of ANN
 **/

/* compute coarse centroids. 
 * centroids should remain valid until the ann_vw_t is deleted */
struct ann_vw_t * ann_vw_new(int d,int nbvw, const float * centroids,int nbvw_coarse,int max_iter,
                             int **clust_assign_out);

/* perform assignment on level 1 (leaves) or 2 (coarse) for learning */
void ann_vw_learn_assign(const struct ann_vw_t *,int level,float *points,int n,int *clust_assign);

/* fill in edges table from co-occurrences */
void ann_vw_make_edges(struct ann_vw_t *,const int *clust_assign1,const int *clust_assign2,int n,int thresh);

/* write centroids and edges */
void ann_vw_write(const struct ann_vw_t *,const char * fclustcoarse_name, const char * fedges_name);


#endif

/*---------------------------------------------------------------------------*/

