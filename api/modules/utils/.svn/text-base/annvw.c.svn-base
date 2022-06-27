
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <string.h>

#include <yael/vector.h>
#include <yael/nn.h>
#include <yael/sorting.h>
#include <yael/clustering.h>
#include <yael/machinedeps.h>

#include "annvw.h"




#define NEWA(type,n) (type*)malloc((n)*sizeof(type))


/*
 * Computations are done by blocks (nice for cache access, distance matrix must fit in mem and allows beautiful progress bar)
 * blocks are BLOCK_N1 * BLOCK_N2
 */

#define BLOCK_N1 8192
#define BLOCK_N2 8192


#define MIN(a,b) ((a)<(b) ? (a) : (b))



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



/**********************************************************************************
 * searching the ANN structure 
 */


/* Load the ANN structure from two files: the coarse centroids and the edges 
   between the two layers of the structure. The parameter centroids contains 
   the vectors of the fine-grain centroids                                          */
ann_vw_t *load_ann_vw (const char *fclustcoarse_name, const char *fedges_name,
                       int nbvw, float *centroids)
{
  int i;
  ann_vw_t *ann = (ann_vw_t *) malloc (sizeof (ann_vw_t));


  FILE *fedges = fopen (fedges_name, "r");
  if (!fedges) {
    fprintf (stderr, "# Unable to open edges file %s\n", fedges_name);
    exit (4);
  }

  /* Load the coarse centroids */
  ann->nbvw_coarse =
      fvecs_new_read (fclustcoarse_name, &(ann->d), &(ann->nodes));

  ann->nbvw = nbvw;
  ann->leaves = centroids;

  /* Load the edges */
  ann->n_edge = NEWA (int, ann->nbvw_coarse);
  ann->edges = NEWA (int *, ann->nbvw_coarse);

  for (i = 0; i < ann->nbvw_coarse; i++) {
    int n;
    if (fread (&n, sizeof (int), 1, fedges) != 1) {
      fprintf (stderr, "load_ann_vw error 1");
      return NULL;
    }
    ann->n_edge[i] = n;
    ann->edges[i] = NEWA (int, n);
    if (fread (ann->edges[i], sizeof (int), n, fedges) != n) {
      fprintf (stderr, "load_ann_vw error 2");
      return NULL;
    }
  }

  fclose (fedges);

  ann->leave_blocks=NULL;

  return ann;
}


void ann_vw_make_blocks(ann_vw_t *ann) { 
  int i,j;

  ann->leave_blocks=NEWA(float*,ann->nbvw_coarse);
  
  for(i=0;i<ann->nbvw_coarse;i++) {
    int n=ann->n_edge[i];
    float *block=ann->leave_blocks[i]=fvec_new((ann->d+1)*n*sizeof(float));
    /* extra space used to store norms */
    float *norms=block+ann->d*n;

    for(j=0;j<n;j++) {
      memcpy(block+j*ann->d,
             ann->leaves+ann->d*ann->edges[i][j],
             sizeof(float)*ann->d);
      norms[j]=fvec_norm2sqr(block+j*ann->d,ann->d);
    }

  }

}


void free_ann_vw (ann_vw_t * ann)
{
  int i;
  free (ann->nodes);
  free (ann->n_edge);
  if(ann->edges) {
    for (i = 0; i < ann->nbvw_coarse; i++)
      free (ann->edges[i]);
    free (ann->edges);
  }

  if(ann->leave_blocks) {
    for (i = 0; i < ann->nbvw_coarse; i++)
      free (ann->leave_blocks[i]);
    free (ann->leave_blocks);
  }

  free (ann);
}

void ann_vw_describe(const ann_vw_t *ann) {
  int i;
  printf("ANN: nbvw=%d d=%d coarse_nbvw=%d\n",
         ann->nbvw,ann->d,ann->nbvw_coarse);
  long tot_edge=0;

  for(i=0;i<ann->nbvw_coarse;i++) 
    tot_edge+=ann->n_edge[i];

  printf("Average connections per node = %.2f\n",tot_edge/(double)ann->nbvw_coarse);
  printf("Connectivity = %ld / %ld = %.5f\n",
         tot_edge,ann->nbvw*(long)ann->nbvw_coarse,
         tot_edge/(double)(ann->nbvw*ann->nbvw_coarse));
  
}



typedef struct {
  void (*fun) (void *arg, double frac);
  void *arg;
  double a, b;
} sub_peek_fun_t;

static void sub_peek_fun (void *arg, double frac)
{
  sub_peek_fun_t *spf = arg;
  (*spf->fun) (spf->arg, spf->b + spf->a * frac);
}


/* Assign a cluster using the two-level structure of the ANN */
void quantize_codebook_annvw_single_full (const ann_vw_t * ann, int npt,
                                          const float *points, int *vw,float *distances,
                                          int *n_scalprod_out,
                                          void (*peek_fun) (void *arg, double frac),
                                          void *peek_arg)
{

  assert(ann->edges || !"call ann_vw_make_edges before using the ann_vw_t");

  int i, j;

  /* first use the top level to find the internal nodes */
  int *switch_assign = NEWA (int, npt);

  sub_peek_fun_t spf = { peek_fun, peek_arg, 0.3, 0 };

  quantize_codebook (npt, ann->nbvw_coarse, ann->d, ann->nodes, points,
                     switch_assign, peek_fun ? &sub_peek_fun : NULL, &spf);

  int n_scalprod=ann->nbvw_coarse*npt;
  
  /* Now select the best nodes in the second level */
  for (i = 0; i < npt; i++) {

    int node = switch_assign[i];
    int nb_leaves = ann->n_edge[node];

    /* Compute the distance from the point to the leaves reacheable from internal node */
    float min_dis = 1e9;
    int min_dis_index = -1;

    const float *v = points + i * ann->d;


    if(!ann->leave_blocks) {

      for (j = 0; j < nb_leaves; j++) {
        const float *leaf = ann->leaves + ann->edges[node][j] * ann->d;
        float dis=fvec_distance_L2sqr(leaf,v,ann->d);
        
        n_scalprod++;
        
        if (dis < min_dis) {
          min_dis = dis;
          min_dis_index = ann->edges[node][j];
        }
      }

    } else {
      const float *block = ann->leave_blocks[node];

      float *dists=fvec_new_copy(block+nb_leaves*ann->d,nb_leaves);
      
      { /* compute distances */
        float minus_2=-2,one=1;
        int inc=1;
        int d=ann->d;
        sgemv_("Transposed",&d,&nb_leaves,&minus_2,block,&d,v,&inc,&one,dists,&inc);
      }

      for (j = 0; j < nb_leaves; j++) {
        n_scalprod++;
      
        float dis=dists[j];
        
        if (dis < min_dis) {
          min_dis = dis;
          min_dis_index = ann->edges[node][j];
        }
      }
      min_dis+=fvec_norm2sqr(v,ann->d);
      free(dists);
    }

    /* assign the vector to the correct centroid */
    vw[i] = min_dis_index;
    distances[i]=min_dis;

    if (peek_fun && i % 128 == 0)
      (*peek_fun) (peek_arg, 0.3 + 0.7 * i / npt);

  }

  free (switch_assign);

  if(n_scalprod_out) *n_scalprod_out=n_scalprod;

  if (peek_fun)
    (*peek_fun) (peek_arg, 1);

}



void quantize_codebook_annvw_full (const ann_vw_t * ann, int npt,
                                   const float *points, int knn, int knn_coarse,
                                   int *vw, float *vwdis, 
                                   int *n_scalprod_out,
                                   void (*peek_fun) (void *arg, double frac),
                                   void *peek_arg)
{
  assert(ann->edges || !"call ann_vw_make_edges before using the ann_vw_t");

  if(knn==1 && knn_coarse==1) {
    quantize_codebook_annvw_single_full(ann,npt,points,vw,vwdis,n_scalprod_out,peek_fun,peek_arg);
    return;
  }

  int i, i2, j;

  /* first use the top level to find the internal nodes */
  int *switch_assign = ivec_new(npt*(long)knn_coarse);

  sub_peek_fun_t spf = { peek_fun, peek_arg, 0.3, 0 };

  float *sub_distances = 
    quantize_codebook_multiple (npt, ann->nbvw_coarse, ann->d, knn_coarse, ann->nodes, points,
                                switch_assign, peek_fun ? &sub_peek_fun : NULL, &spf);
  free(sub_distances);

  int n_scalprod=ann->nbvw_coarse*npt;

  /* Now select the best nodes in the second level */

  maxheap_t *mh = maxheap_new (knn);

  for (i = 0; i < npt; i++) {

    mh->i = 0;
    const float *v = points + i * (long)ann->d;

    for(i2=0;i2<knn_coarse;i2++) {

      int node = switch_assign[i*knn_coarse+i2];
      int nb_leaves = ann->n_edge[node];
      
      /* Compute the distance from the point to the leaves reacheable from internal node */
      
      for (j = 0; j < nb_leaves; j++) {
        const float *leaf = ann->leaves + ann->edges[node][j] * ann->d;
        float dis=fvec_distance_L2sqr(leaf,v,ann->d);

        n_scalprod++;
        maxheap_add (mh, ann->edges[node][j], dis);

      }
    }

    /*
       if(mh->i<knn) 
       fprintf(stderr,"quantize_codebook_annvw_multiple: Warning! not enough edges for multiple assignement. "
       "Some vw are set to -1.\n");
     */
    /* assign the vector to the correct centroid */
    for (j = 0; j < mh->i; j++) {
      vw[i * knn + j] = mh->elts[j].label;
      vwdis[i * knn + j] = mh->elts[j].val;
    }

    for (j = mh->i; j < knn; j++) {
      vw[i * knn + j] = -1;
      vwdis[i * knn + j] = 1e9;
    }

    if (peek_fun && i % 128 == 0)
      (*peek_fun) (peek_arg, 0.3 + 0.7 * i / npt);

  }

  maxheap_delete (mh);

  free (switch_assign);

  if(n_scalprod_out) *n_scalprod_out=n_scalprod;

  if (peek_fun)
    (*peek_fun) (peek_arg, 1);

}



/************* simplified calls */


void quantize_codebook_annvw (const ann_vw_t * ann, int npt,
                              const float * points, int * vw,
                              void (*peek_fun) (void *arg, double frac),
                              void *peek_arg) {

  float *dists=fvec_new(npt);
  quantize_codebook_annvw_single_full (ann,npt,points,vw,dists,NULL,peek_fun,peek_arg);
  free(dists);
}


float* quantize_codebook_annvw_multiple (const ann_vw_t * ann, int npt, float * points,
                                         int k, int * vw) {
  float *dists=fvec_new(npt*k);
  quantize_codebook_annvw_full (ann, npt, points, k, 1, vw, dists, NULL, NULL, NULL);
  return dists;
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

  /* for ANN */
  int k_coarse;
  const ann_vw_t *ann;
  int *n_scalprod;

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

  if (t->ann)
    quantize_codebook_annvw_full (t->ann, n1 - n0, t->points + n0 * t->ann->d,
                                  t->k,t->k_coarse,
                                  t->vw + n0 * t->k, t->vwdis + n0 * t->k,
                                  t->n_scalprod ? t->n_scalprod+i : NULL,
                                  peek_fun, t);
  else 
    quantize_codebook_full (n1 - n0, t->nclust, t->d, t->k, t->codebook,
                            t->points + n0 * t->d, t->vw_weights,
                            t->vw + n0 * t->k, t->vwdis + n0 * t->k,
                            peek_fun, t);
}

void quantize_codebook_annvw_full_thread (const ann_vw_t * ann, int npt,
                                          const float *points, int knn, int knn_coarse, 
                                          int *vw, float *vwdis, 
                                          int *n_scalprod_out,
                                          int n_thread,
                                          void (*peek_fun) (void *arg, double frac),
                                          void *peek_arg)
{
  if (npt < n_thread || n_thread == 1) {        /* too few pts */
    quantize_codebook_annvw_full (ann, npt, points, knn, knn_coarse,vw,vwdis,n_scalprod_out,peek_fun, peek_arg);
    return;
  }
  
  int *n_scalprod_tab=n_scalprod_out ? ivec_new(n_thread) : NULL;
  
  /* divide assignment into slices */
  quantize_codebook_input_t task = { 
    0, 0, knn, NULL, 
    knn_coarse,ann,n_scalprod_tab,
    npt, points, NULL, vw, vwdis,
    n_thread, peek_fun, peek_arg, PTHREAD_MUTEX_INITIALIZER
  };

  compute_tasks (n_thread, n_thread, &quantize_codebook_task, &task);
  
  if(n_scalprod_tab) {
    *n_scalprod_out=ivec_sum(n_scalprod_tab,n_thread);
    free(n_scalprod_tab);
  }
  
}


void quantize_codebook_annvw_thread (const ann_vw_t * ann, int npt,
                                     const float * points, int * vw,
                                     int n_thread,
                                     void (*peek_fun) (void *arg, double frac),
                                     void *peek_arg) {
  float *vwdis=fvec_new(npt);
  quantize_codebook_annvw_full_thread (ann,npt,points,1,1,vw,vwdis,NULL,n_thread,peek_fun,peek_arg);
  free(vwdis);
}


float* quantize_codebook_annvw_multiple_thread (const ann_vw_t * ann, int npt, float * points,
                                                int k, int * vw, int n_thread) {
  float *dists=fvec_new(npt*k);
  quantize_codebook_annvw_full_thread (ann, npt, points, k, 1, vw, dists, NULL, n_thread, NULL, NULL);
  return dists;
}



/**********************************************************************************
 * Building the ANN structure 
 */

ann_vw_t * ann_vw_new(int d,int nbvw, const float * centroids,int nbvw_coarse,int max_iter,
                      int **clust_assign_out) {
  ann_vw_t *ann = (ann_vw_t *) malloc (sizeof (ann_vw_t));

  ann->d=d;
  ann->nbvw=nbvw;
  ann->nbvw_coarse=nbvw_coarse;
  int *clust_assign;
  ann->nodes=clustering_kmeans_assign (nbvw, d, centroids, nbvw_coarse, max_iter, 0, &clust_assign);
  if(clust_assign_out) 
    *clust_assign_out=clust_assign;
  else
    free(clust_assign);
  
  ann->leaves=centroids;

  ann->n_edge=ivec_new(nbvw_coarse);
  ann->edges=NULL;
  ann->leave_blocks=NULL;
  return ann;
}


/* perform assignment on level 0 (coarse) or 1 (leaves) */
void ann_vw_learn_assign(const ann_vw_t *ann,int level,float *points,int n,int *clust_assign) {

  if(level==1) 
    quantize_codebook_thread (n,ann->nbvw,ann->d,ann->leaves,points,clust_assign,count_cpu(),NULL,NULL);
  else if(level==2)
    quantize_codebook_thread (n,ann->nbvw_coarse,ann->d,ann->nodes,points,clust_assign,count_cpu(),NULL,NULL);
  else assert(0);

}



/* fill in edges table from co-occurrences */
void ann_vw_make_edges(ann_vw_t *ann,const int *clust_assign1,const int *clust_assign2,int n,int thresh) {
  int i,j;
  long nbvw2=ann->nbvw_coarse;
  long nbvw1=ann->nbvw; 

  if(ann->edges) {
    for(i=0;i<nbvw2;i++) free(ann->edges[i]);
    free(ann->edges);
  }
  ann->edges=NEWA(int*,nbvw2);

  int *accu=ivec_new_0(nbvw1*nbvw2);
  /* TODO transpose */
#define M(i,j) (accu[(i)*nbvw2+(j)])

  for (i = 0 ; i < n ; i++) {
    M (clust_assign1[i], clust_assign2[i])++;
  }    

  int *edges=ivec_new(nbvw1);
  for (i = 0 ; i < nbvw2 ; i++) {
    int ne=0;
    
    for (j = 0 ; j < nbvw1 ; j++)
      if (M(j,i) >= thresh)
        edges[ne++]=j;

    ann->n_edge[i]=ne;
    ann->edges[i]=ivec_new_copy(edges,ne); /* don't trust realloc */
  }
 
  free(edges);
  free(accu);
#undef M

  
}


/* write centroids and edges */
void ann_vw_write(const ann_vw_t *ann,const char * fclustcoarse_name, const char * fedges_name) {
  int i;

  int nbvw2=ann->nbvw_coarse;
  
  fvecs_write (fclustcoarse_name, ann->d, nbvw2, ann->nodes);
  
  FILE * fedges = fopen (fedges_name, "w");

  if(!fedges) {
    fprintf(stderr,"could not open %s",fedges_name);
    perror("");
    return;
  }

  for(i=0;i<nbvw2;i++) {
    int ne=ann->n_edge[i];
    fwrite(&ne,sizeof(int),1,fedges);
    fwrite(ann->edges[i],sizeof(int),ne,fedges);
  }
  
  fclose(fedges);
}
