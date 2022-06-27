
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <string.h>

#include "vector.h"
#include "nn.h"
#include "sorting.h"




#define NEWA(type,n) (type*)malloc((n)*sizeof(type))


/**********************************************************************************
 * Distance functions  
 */


/*------------------ Blas subroutine ------------------*/

#define real float
#define integer int

int sgemm_ (char *transa, char *transb, integer * m, integer *
            n, integer * k, real * alpha, const real * a, integer * lda,
            const real * b, integer * ldb, real * beta, real * c__,
            integer * ldc);


int sgemv_(char *trans, integer *m, integer *n, real *alpha, 
           const real *a, integer *lda, const real *x, integer *incx, real *beta, real *y, 
           integer *incy);

#undef real
#undef integer


/*
 * computes dist2 := dist2 - 2 * descs * clusters' 
 * where 
 *   dist2    is ndesc-by-nclust 
 *   clusters is nclust-by-d
 *   descs    is ndesc-by-d
 * (all matrices stored by lines, à la C, and packed)
 */
static void add_matmul (int d, int na, int nb,
                        const float *a, int lda, 
                        const float *b, int ldb,
                        float *dist2)
{

  float minus_two = -2;
  float one = 1;

  sgemm_ ("Transposed", "Not trans", &na, &nb, &d,
          &minus_two, a, &lda, b, &ldb, &one, dist2, &na);
}


/* computes all distances between a line of a and a line of b. 
 *   a(na,d) by lines
 *   b(nb,d) by lines
 *  dist2[i+na*j] = || a(i,:)-b(j,:) ||^2
 */
void compute_cross_distances (int d, int na, int nb,
                              const float *a, const float *b,  
                              float *dist2) 
{
  compute_cross_distances_nonpacked (d, na, nb, a, d, b, d, dist2);
}

void compute_cross_distances_nonpacked (int d, int na, int nb,
                                        const float *a, int lda,
                                        const float *b, int ldb, 
                                        float *dist2)
{
  long i, j;
  float *sum_c2 = (float *) malloc (sizeof (float) * na);

  for (i = 0; i < na; i++) {
    float s = 0;
    const float *cl = a + lda * i;
    for (j = 0; j < d; j++)
      s += cl[j] * cl[j];
    sum_c2[i] = s;
  }

  for (i = 0; i < nb; i++) {
    double sum_d2 = 0;
    const float *dl = b + ldb * i;
    for (j = 0; j < d; j++)
      sum_d2 += dl[j] * dl[j];
    float *d2l = dist2 + i * na;
    for (j = 0; j < na; j++)
      d2l[j] = sum_d2 + sum_c2[j];
  }

  add_matmul (d, na, nb, a, lda, b, ldb, dist2);

  free (sum_c2);
}



static double sqr (double x)
{
  return x * x;
}

/* alternative distance functions (not optimized) */

void compute_cross_distances_alt (int distance_type, int d, int na, int nb,
                                  const float *a,
                                  const float *b, float *dist2) {
  int i,j,k;

  for(i=0;i<na;i++) for(j=0;j<nb;j++) {
    double sum=0;
    if(distance_type==1) 
      for(k=0;k<d;k++) sum+=fabs(a[k+i*d]-b[k+j*d]);
    else if(distance_type==2) 
      for(k=0;k<d;k++) sum+=sqr(a[k+i*d]-b[k+j*d]);
    else if(distance_type==3) 
      for(k=0;k<d;k++) {
        float av=a[k+i*d],bv=b[k+j*d];
        sum+=av+bv==0 ? 0 : sqr(av+bv)/(av+bv);
      }      

    dist2[i+j*na]=sum;
  }

}


/**********************************************************************************
 * Elementary cluster assignment 
 */


/*
 * Computations are done by blocks (nice for cache access, distance matrix must fit in mem and allows beautiful progress bar)
 * blocks are BLOCK_NPT * BLOCK_CLUST
 */

#define BLOCK_N1 8192

#define BLOCK_N2 8192


#define BLOCK_NPT BLOCK_N1
#define BLOCK_CLUST BLOCK_N2




#define MIN(a,b) ((a)<(b) ? (a) : (b))

/* This function quantizes the vectors coords according to the codebook clusters, 
   and sets the corresponding indexes in the index vector vw accordingly    
 */

/* n1 = pts */
static void quantize_codebook_single_full (int n1, int n2, int d,
                                           const float *mat2, const float *mat1, 
                                           const float *vw_weights,                             
                                           int *vw, float *vwdis,
                                           void (*peek_fun) (void *arg, double frac),
                                           void *peek_arg)
{
  int step1 = MIN (n1, BLOCK_N1), step2 = MIN (n2, BLOCK_N2);

  float *dists = fvec_new (step1 * step2);

  /* divide the dataset into sub-blocks to:
   * - not make a too big dists2 output array 
   * - call peek_fun from time to time
   */
  
  long i1,i2,j1,j2;
  for (i1 = 0; i1 < n1; i1 += step1) {  

    int m1 = MIN (step1, n1 - i1);

    /* clear mins */

    for (j1 = 0; j1 < m1; j1++) {
      vw[j1+i1]=-1;
      vwdis[j1+i1]=1e30;
    }

    for (i2 = 0; i2 < n2 ; i2 += step2) {     
      
      int m2 = MIN (step2, n2 - i2);
      
      compute_cross_distances (d, m2, m1, mat2+i2*d, mat1+i1*d, dists);

      if(vw_weights) {
        for(j1=0;j1<m1;j1++) for (j2 = 0; j2 < m2; j2++)
          dists[j1 * m2 + j2] *= vw_weights[j2 + i2];        
      }

      /* update mins */

      for(j1=0;j1<m1;j1++) {
        float *dline=dists+j1*m2;
        
        int imin=vw[i1+j1];
        float dmin=vwdis[i1+j1];

        for(j2=0;j2<m2;j2++) 
          if(dline[j2]<dmin) {
            imin=j2+i2;
            dmin=dline[j2];
          }
          
        vw[i1+j1]=imin;
        vwdis[i1+j1]=dmin;

      }      

    }  

    if (peek_fun)
      (*peek_fun) (peek_arg, i1 / (double) n1);

  }

  free (dists);

  if (peek_fun)
    (*peek_fun) (peek_arg, 1);
}




void quantize_codebook_full (int n1, int n2, int d, int k,
                             const float *mat2, const float *mat1,
                             const float *vw_weights,
                             int *vw, float *vwdis,                                             
                             void (*peek_fun) (void *arg,double frac),
                             void *peek_arg)
{
  assert (k <= n2);

  if(k==1) {
    quantize_codebook_single_full(n1, n2, d, mat2, mat1, vw_weights, vw, vwdis, 
                                  peek_fun, peek_arg);
    return;
  }

  
  int step1 = MIN (n1, BLOCK_N1), step2 = MIN (n2, BLOCK_N2);

  float *dists = fvec_new (step1 * step2);


  /* allocate all heaps at once */
  long oneh = sizeof (maxheap_t) + k * sizeof (heap_entry_t);
  // oneh=(oneh+7) & ~7; /* round up to 8 bytes */
  char *minbuf = malloc ((oneh) * step1);

#define MINS(i) ((maxheap_t*)(minbuf + oneh * i))
  
  
  long i1,i2,j1,j2;
  for (i1 = 0; i1 < n1; i1 += step1) {  

    int m1 = MIN (step1, n1 - i1);

    /* clear mins */
    for (j1 = 0; j1 < m1; j1++) {
      MINS(j1)->n = k;
      MINS(j1)->i = 0;
    }
    

    for (i2 = 0; i2 < n2 ; i2 += step2) {     
      
      int m2 = MIN (step2, n2 - i2);
      
      compute_cross_distances (d, m2, m1, mat2+i2*d, mat1+i1*d, dists);

      if(vw_weights) {
        for(j1=0;j1<m1;j1++) for (j2 = 0; j2 < m2; j2++)
          dists[j1 * m2 + j2] *= vw_weights[j2 + i2];        
      }

      /* update mins */

      for(j1=0;j1<m1;j1++) {
        float *dline=dists+j1*m2; 
        maxheap_add_multiple(MINS(j1),i2,m2,dline);
      }      

    }  

    for (j1 = 0; j1 < m1; j1++) {
      maxheap_t *mh = MINS(j1);
      assert (mh->i == k);
      
      for (j2 = 0; j2 < k; j2++) {
        vw[(i1+j1) * k + j2] = mh->elts[j2].label;
        vwdis[(i1+j1) * k + j2] = mh->elts[j2].val;
      }
    }

    if (peek_fun)
      (*peek_fun) (peek_arg, i1 / (double) n1);

  }

#undef MINS
  free (minbuf);
  free(dists);

  if (peek_fun)
    (*peek_fun) (peek_arg, 1);

}

/**********************************************************************************
 * Simple call versions
 */


void quantize_codebook (int npt, int nclust, int d,
                        const float *codebook, const float *coords, int *vw,
                        void (*peek_fun) (void *arg, double frac),
                        void *peek_arg) {
  
  /* The distances to centroids that will be returned */
  float *vwdis = fvec_new(npt);
  
  quantize_codebook_full (npt, nclust, d, 1, codebook, coords, NULL, vw, vwdis, peek_fun, peek_arg);

  free(vwdis);
}

float *quantize_codebook_multiple (int npt, int nclust, int d, int k,
                                   const float *codebook, const float *coords,
                                   int *vw, void (*peek_fun) (void *arg,
                                                              double frac),
                                   void *peek_arg)
{
  /* The distances to centroids that will be returned */
  float *vwdis = fvec_new(npt * k);

  quantize_codebook_full (npt, nclust, d, k, codebook, coords, NULL, vw, vwdis, peek_fun, peek_arg);
  
  return vwdis;
}



/**********************************************************************************
 * Threaded versions
 */

#include <pthread.h>

/* a common function dispatches the calls */
typedef struct {

  /* input */

  /* for normal assignement */
  int nclust, d, k;
  const float *codebook;

/*   /\* for ANN *\/ */
/*   int k_coarse; */
/*   const ann_vw_t *ann; */
/*   int *n_scalprod; */

  int npt;
  const float *points;
  const float *vw_weights;

  /* output */
  int *vw;
  float *vwdis;

  /* bookkeeping */
  int n_thread;
  void (*peek_fun) (void *arg, double frac);
  void *peek_arg;
  pthread_mutex_t peek_mutex;
} quantize_codebook_input_t;

static void quantize_codebook_peek (void *arg, double frac)
{
  quantize_codebook_input_t *t = arg;
  if (pthread_mutex_trylock (&t->peek_mutex) != 0)      /* don't insist if another thread is peeking */
    return;
  (*t->peek_fun) (t->peek_arg, frac);
  pthread_mutex_unlock (&t->peek_mutex);
}


static void quantize_codebook_task (void *arg, int tid, int i)
{
  quantize_codebook_input_t *t = arg;

  long n0 = t->npt * (long)i / t->n_thread;
  long n1 = t->npt * (long)(i + 1) / t->n_thread;

  void (*peek_fun) (void *arg, double frac) =
      t->peek_fun ? &quantize_codebook_peek : NULL;

  quantize_codebook_full (n1 - n0, t->nclust, t->d, t->k, t->codebook,
			  t->points + n0 * t->d, t->vw_weights,
			  t->vw + n0 * t->k, t->vwdis + n0 * t->k,
			  peek_fun, t);
}

/********** frontends */

void quantize_codebook_full_thread (int npt, int nclust, int d, int k,
                                    const float *codebook, const float *coords,
                                    const float *vw_weights,
                                    int *vw, float *vwdis2,
                                    int n_thread,
                                    void (*peek_fun) (void *arg,double frac),
                                    void *peek_arg) {
  if (npt < n_thread || n_thread == 1) {        /* too few pts */
    return quantize_codebook_full (npt, nclust, d, k, codebook, coords, vw_weights, vw, vwdis2, peek_fun, peek_arg);
  }

  quantize_codebook_input_t task = { 
    nclust, d, k, codebook, 
    npt, coords, vw_weights, vw, vwdis2,
    n_thread, peek_fun, peek_arg, PTHREAD_MUTEX_INITIALIZER
  };

  compute_tasks (n_thread, n_thread, &quantize_codebook_task, &task);

}


/***************** simplified calls */

float *quantize_codebook_multiple_thread (int npt, int nclust, int d, int k,
                                          const float *codebook, const float *coords,
                                          int *vw, 
                                          int n_thread,
                                          void (*peek_fun) (void *arg,
                                                                     double frac),
                                          void *peek_arg) {
  float *vwdis2=fvec_new(k*npt);
    
  quantize_codebook_full_thread (npt, nclust, d, k, codebook, coords, NULL, vw, vwdis2, n_thread, peek_fun, peek_arg);
  
  return vwdis2;
}



void quantize_codebook_thread (int npt, int nclust, int d,
                               const float *codebook, const float *coords,
                               int *vw, int n_thread,
                               void (*peek_fun) (void *arg, double frac),
                               void *peek_arg)
{
  float *vwdis2=fvec_new(npt);

  quantize_codebook_full_thread (npt, nclust, d, 1, codebook, coords, NULL, vw, vwdis2, n_thread, peek_fun, peek_arg);
  
  free(vwdis2);
}



/***********************************************************************
 *           Implementation of the threading part
 *
 * generic thread stuff */

typedef struct {
  pthread_mutex_t mutex;
  int i, n, tid;
  void (*task_fun) (void *arg, int tid, int i);
  void *task_arg;
} context_t;



static void *start_routine (void *cp)
{
  context_t *context = cp;
  int tid;
  pthread_mutex_lock (&context->mutex);
  tid = context->tid++;
  pthread_mutex_unlock (&context->mutex);

  for (;;) {
    int item;
    pthread_mutex_lock (&context->mutex);
    item = context->i++;
    pthread_mutex_unlock (&context->mutex);
    if (item >= context->n)
      break;
    else
      context->task_fun (context->task_arg, tid, item);
  }

  return NULL;
}

void compute_tasks (int n, int nthread,
                    void (*task_fun) (void *arg, int tid, int i),
                    void *task_arg)
{
  int i;
  context_t context;

  assert(nthread>0 || !"sombody has to do the job");

  context.i = 0;
  context.n = n;
  context.tid = 0;
  context.task_fun = task_fun;
  context.task_arg = task_arg;

  pthread_mutex_init (&context.mutex, NULL);

  if(nthread==1) 
    start_routine(&context);    
  else {
    pthread_t *threads = malloc (sizeof (pthread_t) * n);
      
    for (i = 0; i < nthread; i++) 
      pthread_create (&threads[i], NULL, &start_routine, &context);
      
    /* all computing */
      
    for (i = 0; i < nthread; i++)
      pthread_join (threads[i], NULL);
    
    free (threads);
  }

}


