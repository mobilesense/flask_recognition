/*---------------------------------------------------------------------------*/

#ifndef NN_H_INCLUDED
#define NN_H_INCLUDED

/*---------------------------------------------------------------------------*/
/*! @addtogroup utils
 *  @{
 */
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/* Cluster assignment                                                        */
/*---------------------------------------------------------------------------*/





/*! @brief Quantizes the vectors coords according to the codebook 
 * clusters,  and sets the set of corresponding indexes in the index vector 
 * vw accordingly  
 * 
 *     codebook is nclust*d   
 *     coords is npt*d        
 *     vws corresponding to i are vw[i*k] to vw[i*k+k-1].
 *
 * all stored in row-major order
 *
 *     k:           nb of multiple assignments 
 *     vw_weigths:  multiply squared distances by this for each centroid
 *     vwdis2:      squared distances for each point
 * 
 * pek_fun needs not to be reentrant 
 * falls back on mono-thread version if task too small.
 */
void quantize_codebook_full (int npt, int nclust, int d, int k,
                             const float *codebook, const float *coords,
                             const float *vw_weights,
                             int *vw, float *vwdis2,                                             
                             void (*peek_fun) (void *arg,double frac),
                             void *peek_arg);

/* multi-threaded version */

void quantize_codebook_full_thread (int npt, int nclust, int d, int k,
                                    const float *codebook, const float *coords,
                                    const float *vw_weights,
                                    int *vw, float *vwdis2,
                                    int n_thread,
                                    void (*peek_fun) (void *arg,double frac),
                                    void *peek_arg);


/* next functions are simplified calls of the previous */

void quantize_codebook (int npt, int nclust, int d,
                        const float *codebook, const float *coords, int *vw,
                        void (*peek_fun) (void *arg, double frac),
                        void *peek_arg);


/*! @brief Threaded version.  */
void quantize_codebook_thread (int npt, int nclust, int d,
                               const float *codebook, const float *coords, int *vw,
                               int n_thread,
                               void (*peek_fun) (void *arg, double frac),
                               void *peek_arg);

/*! @brief function returns the set of distances to centroids (alloc'ed with malloc)
 */
float* quantize_codebook_multiple (int npt, int nclust, int d, int k,
                                   const float *codebook, const float *coords, int *vw,
                                   void (*peek_fun) (void *arg, double frac),
                                   void *peek_arg);


float* quantize_codebook_multiple_thread (int npt, int nclust, int d, int k,
                                          const float *codebook, const float *coords, int *vw,
                                          int n_thread,
                                          void (*peek_fun) (void *arg, double frac),
                                          void *peek_arg);





/*! @brief Low-level function to compute all distances between 2 sets of vectors 
 *
 *  a is na*d \n
 *  b is nb*d \n
 *  dist2[i+na*j] = || a(i,:)-b(j,:) ||^2 \n
 *  uses BLAS if available
 */
void compute_cross_distances (int d, int na, int nb,
                              const float *a,
                              const float *b, float *dist2);

void compute_cross_distances_nonpacked (int d, int na, int nb,
                                        const float *a, int lda,
                                        const float *b, int ldb, 
                                        float *dist2);


/*! @ alternative distances. 
 *
 * distance_type==1: L1, ==3: symmetric chi^2  
 */

void compute_cross_distances_alt (int distance_type, int d, int na, int nb,
                                  const float *a,
                                  const float *b, float *dist2);


/* compute_tasks creates nthread threads that call task_fun n times 
 * with arguments:
 *   arg=task_arg
 *   tid=identifier of the thread in 0..nthread-1
 *   i=call number in 0..n-1
 */

void compute_tasks (int n, int nthread,
                    void (*task_fun) (void *arg, int tid, int i),
                    void *task_arg);


#endif

/*---------------------------------------------------------------------------*/

