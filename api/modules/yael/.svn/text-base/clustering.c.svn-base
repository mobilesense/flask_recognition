#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <sys/time.h>

#include "clustering.h"
#include "vector.h"
#include "nn.h"
#include "machinedeps.h"
#include "sorting.h"


#define NEWA(type,n) (type*)malloc(sizeof(type)*(n))
#define NEWAC(type,n) (type*)calloc(sizeof(type),(n))
#define NEW(type) NEWA(type,1)

static void floats_mul (float *a, float b, int n)
{
  while (n--)
    a[n] *= b;
}

static double sqr (double a)
{
  return a * a;
}


/*
 * static double floats_assign_dist2 (float * centroids, float * points, int * clust_assign,
 * 	int n, int d)
 * {
 *     double accu=0;
 *     while (n--) 
 * 	accu += floats_dist2(points+(n*d),centroids+(clust_assign[n]*d),d);
 *     return accu;
 * }
 */


static void floats_normalize2 (float *a, int n, int d, int newnorm)
{
  int i, j;
  double vecnorm, f;

  for (i = 0; i < n; i++) {
    float *al = a + d * i;

    /* compute the norm */
    vecnorm = 0;
    for (j = 0; j < d; j++)
      vecnorm += al[j] * al[j];

    f = newnorm / sqrt (vecnorm);

    /* Normalize the corresponding vector */
    for (j = 0; j < d; j++)
      al[j] *= f;
  }
}


static void progress_bar_peek_fun (void *arg, double frac)
{
  int n = (int) (frac * 50), i;
  for (i = 0; i < 50 ; i++)
    putc (i < n ? '+' : '-',stderr);
  putc ('\r',stderr);
  fflush (stderr);
}



float *clustering_kmeans_assign_with_score (int n, int di,
                                            const float *points, int k, int nb_iter_max,                                            
                                            double normalize, 
                                            int n_thread,
                                            double *score, int **clust_assign_out)
{

  assert (n > k);

  int i, c, iter = 0;
  long d=di; /* to force 64-bit address computations */ 

  float *centroids = fvec_new (k * d);
  float *centroids_new = fvec_new (k * d);

  /* First select a subset of vectors */
  int *randperm = ivec_new_random_perm (n);
  for (i = 0; i < k; i++)
    memcpy (centroids + i * d, points + randperm[i] * d, sizeof (float) * d);
  free (randperm);

  /* vector counting the number of vectors in each centroids */
  int *cluster_size = (int *) malloc (sizeof (int) * n);

  /* to which cluster a vector is assigned ? */
  int *clust_assign = (int *) malloc (sizeof (int) * n);
  
  /* the distance between a vector and its cluster centroid */
  float *vwdist = (float *) malloc (sizeof (float) * n);

  /* a value to check if the centroids have been updated */
  double modiftestvalue, modiftestvalue_new;
  if (score) {
    memset (clust_assign, 0, sizeof (int) * n);
    modiftestvalue_new=fvec_norm2sqr (points, n * d);
  } else {
    modiftestvalue_new=fvec_sum (centroids, k * d);
  }
  assert(finite(modiftestvalue_new));

  void (*peek_fun) (void *arg, double frac);

  int verbose=!(n_thread & CLUSTERING_KMEANS_QUIET);
  n_thread &= ~ CLUSTERING_KMEANS_QUIET; 
  peek_fun=verbose ? &progress_bar_peek_fun : NULL; 

  /*  printf("nb_iter_max=%d\n",nb_iter_max); */

  /* launch the construction iterating process */
  do {
    if (normalize != 0)
      floats_normalize2 (centroids, k, d, normalize);

    iter++;

    if(verbose)
      fprintf (stderr, "                                                   [%3d/%-3d]\r", 
               iter, nb_iter_max);
    modiftestvalue = modiftestvalue_new;

    memset (centroids_new, 0, sizeof (float) * k * d);
    memset (cluster_size, 0, sizeof (int) * k);

    if (score) {
      quantize_codebook_full_thread(n,k,d,1,centroids, points,NULL,
                                    clust_assign, vwdist, n_thread, peek_fun, NULL); 
    } else {
      quantize_codebook_thread (n, k, d, centroids, points, clust_assign,
                                n_thread, peek_fun, NULL);
    }

    for (i = 0; i < n; i++) {
      assert((clust_assign[i]>=0 && clust_assign[i]<k) || !"probably a NaN in the input data"); 
      fvec_add (centroids_new + clust_assign[i] * d, points + i * d, d);
      cluster_size[clust_assign[i]]++;
    }

    /* check here because we want clust_assign to be in sync with centroids */
    if (iter == nb_iter_max)
      break;

    /*
       for (i = 0 ; i < k ; i++) 
       fprintf (stderr, "%d ", cluster_size[i]); 
     */
    //    fprintf (stderr, "\n");

    for (i = 0; i < k; i++) {
      if (cluster_size[i] > 0) {
        floats_mul (centroids_new + (i) * d, 1.0 / cluster_size[i], d);
      }
      /* Do no allow a cluster to be void. Split another cluster into two */
      else {
        /* Find a non-void cluster */
        do
          c = floor (drand48 () * k);
        while (cluster_size[c] < 4);    /* minimum size for splitting a cluster */

        /* Add a very small gaussian noise to the previous cluster */

        float noise_amount = 0.00001;
        {
          float *ci = centroids_new + i * d;
          float *cc = centroids_new + c * d;
          int j;
          for (j = 0; j < d; j++)
            ci[j] = cc[j] + (drand48 () * 2 - 1) * noise_amount;
        }

        fprintf (stderr, "r");
      }
    }

    /* make the sum to check if the clusters have changed or not */
    if (score) {
      modiftestvalue_new = fvec_sum (vwdist, n);
    } else {
      modiftestvalue_new = fvec_sum (centroids_new, k * d);
    }
    assert(finite(modiftestvalue_new));

    {
      float *mat_tmp = centroids_new;     
      centroids_new = centroids;
      centroids = mat_tmp;
    }
/*    printf("%g ?= %g\n",modiftestvalue, modiftestvalue_new);  */
  } while (modiftestvalue != modiftestvalue_new);


  {
    int i;
    double sum=0,sum2=0;

    for(i=0;i<k;i++) {
      sum+=cluster_size[i];
      sum2+=cluster_size[i]*(double)cluster_size[i];
    }
    //    fprintf(stderr,"imbalance factor=%.3f\n",sum2*k/(sum*sum));

    if (score) {
	*score = modiftestvalue;
	fprintf(stderr,"sum of sq distances=%g\n",modiftestvalue);
    }
    
  }

  if (normalize)
    floats_normalize2 (centroids, k, d, normalize);

  free (centroids_new);

  if (clust_assign_out)
    *clust_assign_out = clust_assign;
  else
    free (clust_assign);

  free (cluster_size);
  fprintf (stderr, "\n");

  return centroids;
}

float *clustering_kmeans_assign (int n, int d,
                                 const float *points, int k, int nb_iter_max,
                                 double normalize, int **clust_assign_out)
{
   return clustering_kmeans_assign_with_score (n, d, points, k, nb_iter_max,
                                               normalize, count_cpu(), NULL, clust_assign_out);
}

float *clustering_kmeans (int n, int d,
                          const float *points, int k, int nb_iter_max,
                          double normalize)
{

  int *clust_assign;

  float *centroids = clustering_kmeans_assign_with_score (n, d, points, k, nb_iter_max,
                                               normalize, count_cpu(), NULL, &clust_assign);
  free (clust_assign);

  return centroids;
}




static void clustering_cdm_factors(void (*find_knn)(void *arg,float *neigh_dists),
                                   void *find_arg,
                                   int n, 
                                   double alpha, int n_iter,
                                   int n_neigh,
                                   float *cdm_factors) {
  
  int i,j,k;

  /* reverse assignement: we want neighbours of the cluster centroids */ 
  int *centroid_neighs=ivec_new(n*n_neigh);

  /* average distance to neighbourhood */
  float *r=fvec_new(n),*f2=fvec_new(n),*neigh_dists=fvec_new(n*n_neigh);

  for(i=0;i<n;i++) 
    cdm_factors[i]=1.0;

  for(i=0;i<n_iter;i++) {
    
    /**** compute nearest neighbours using distance correction */
   
    /* square correction term because all distances are squared */

    /* 

    for(j=0;j<n;j++) f2[j]=sqr(cdm_factors[j]);    

    quantize_codebook_full_thread (n,n,d,n_neigh,
                                   centroids,centroids,f2,
                                   centroid_neighs,neigh_dists,
                                   count_cpu(),NULL,NULL);

    */

    (*find_knn)(find_arg,neigh_dists);
    
    /* new distances */

    for(j=0;j<n;j++) {
      for(k=0;k<n_neigh;k++) {
        assert(finite(neigh_dists[j*n_neigh+k]));
        float d=neigh_dists[j*n_neigh+k];
        if(d<0) d=0;
        neigh_dists[j*n_neigh+k]=sqrt(d)*cdm_factors[j];
        assert(finite(neigh_dists[j*n_neigh+k]));
      }
    }

    /**** compute average distances to neighbourhood */

    for(j=0;j<n;j++) {
      r[j]=fvec_sum(neigh_dists+j*n_neigh,n_neigh)/n_neigh;
      assert(finite(r[j]));
    }
    /**** compute overall geometrical mean average */

    double sum=0;

    for(j=0;j<n;j++) sum+=log(r[j]);

    double r_bar=exp(sum/n);

    /**** compute distance correction term */ 
    
    for(j=0;j<n;j++) 
      cdm_factors[j]*=pow(r_bar/r[j],alpha);

  }

  free(r);
  free(f2);
  free(neigh_dists);  
  free(centroid_neighs);

}


typedef struct {
  int n,d,n_neigh;
  const float *centroids;
  float *cdm_factors;
} cmd_l2_args_t;


static void cdm_l2_find_knn(void *arg,float *dists) {
  cmd_l2_args_t *a=arg;
  int j;
  int n=a->n,k=a->n_neigh;
  const float *cdm_factors=a->cdm_factors;

  float *f2=fvec_new(n);
  
  for(j=0;j<n;j++) f2[j]=sqr(cdm_factors[j]);    

  int *centroid_neighs=ivec_new(k*n);

  quantize_codebook_full_thread (n,n,a->d,k,
                                 a->centroids,a->centroids,f2,
                                 centroid_neighs,dists,
                                 count_cpu(),NULL,NULL);
  for(j=0;j<n*k;j++) {
    float d=dists[j];
    dists[j]=d<0 ? 0 : sqrt(d);
  }
  
  free(f2);
  free(centroid_neighs); 

}

void clustering_cdm_factors_l2(int n, int d,  
                               const float *centroids, 
                               double alpha, int n_iter,
                               int n_neigh,
                               float *cdm_factors) {
  cmd_l2_args_t args={n,d,n_neigh,centroids,cdm_factors};

  clustering_cdm_factors(&cdm_l2_find_knn,&args,n,
                         alpha, n_iter, n_neigh,
                         cdm_factors);
}

typedef struct {
  int n,n_neigh;
  const float *distances;
  float *cdm_factors;
} cmd_dists_args_t;


static void cdm_dists_find_knn(void *arg,float *dists) {
  cmd_dists_args_t *a=arg;
  int n=a->n,k=a->n_neigh;
  const float *cdm_factors=a->cdm_factors;
  int i,j;
  
  float *wd=fvec_new(n);

  for(i=0;i<n;i++) {

    const float *dl=a->distances+i*n;

    for(j=0;j<n;j++) 
      wd[j]=cdm_factors[j]*dl[j];

    fvec_quantile(wd,n,k);
    
    memcpy(dists+i*k,wd,k*sizeof(float));

  }
  
  free(wd);

}

void clustering_cdm_factors_dists(int n, 
                                  const float *distances, 
                                  double alpha, int n_iter,
                                  int n_neigh,
                                  float *cdm_factors) {

  cmd_dists_args_t args={n,n_neigh,distances,cdm_factors};

  clustering_cdm_factors(&cdm_dists_find_knn,&args,n,
                         alpha, n_iter, n_neigh,
                         cdm_factors);
  
}






static double getmillisecs ()
{
  struct timeval tv;
  gettimeofday (&tv, NULL);
  return tv.tv_sec * 1e3 + tv.tv_usec * 1e-3;
}


double clustering_kmedoids_from_dists (int n, int n_med, int k,
                                       int nb_iter_max, const float *alldists,
                                       int **clust_assign_out,
                                       int **med_subset_out)
{
  int i, j;

  double t0 = getmillisecs ();

#define ALLDISTS(clustno,ptno) alldists[(ptno)+(clustno)*n]

  /* begin with a random subset of vectors */
  int *med_subset = ivec_new_random_perm (n_med);

  /* compute assignement and sse */
  double sse = 0;
  int *clust_assign = NEWA (int, n);

  for (i = 0; i < n; i++) {
    int clus = -1;
    double se = 1e50;
    for (j = 0; j < k; j++) {
      double new_se = ALLDISTS (med_subset[j], i);
      if (new_se < se) {
        se = new_se;
        clus = j;
      }
    }
    sse += se;
    clust_assign[i] = clus;
  }


  /* iterations */
  int *clust_assign_new = NEWA (int, n);
  int iter;
  int n_change = 0;

  for (iter = 0; iter < nb_iter_max; iter++) {

    /* select random swap */
    int to_p = random () % k;
    int from;
    do {                        /* steal from another cluster */
      from = random () % n_med;
    } while (clust_assign[from] == to_p);

    /* compute new sse and cluster assignement */
    double new_sse = 0;
    for (i = 0; i < n; i++) {

      int clust;
      double se;

      if (clust_assign[i] == to_p) {
        /* re-assign completely */
        clust = -1;
        se = 1e50;
        for (j = 0; j < k; j++) {
          int l = j == to_p ? from : med_subset[j];
          double d2 = ALLDISTS (l, i);
          if (d2 < se) {
            se = d2;
            clust = j;
          }
        }
      } else {
        clust = clust_assign[i];
        se = ALLDISTS (med_subset[clust], i);
        /* check if point would be assigned to new cluster */
        double d2 = ALLDISTS (from, i);
        if (d2 < se) {
          se = d2;
          clust = to_p;
        }
      }
      clust_assign_new[i] = clust;
      new_sse += se;
    }
    /*
       printf("replace centroid %d with %d: new sse=%g/%g, %s \n",
       perm[to_p],from,new_sse,sse,new_sse<sse ? "keep" : "reject");
     */

    /* check if improves sse */
    if (new_sse < sse) {
      n_change++;
      sse = new_sse;

      med_subset[to_p] = from;

      {                         /* swap current & new */
        int *tmp = clust_assign_new;
        clust_assign_new = clust_assign;
        clust_assign = tmp;
      }
    }

    if (iter % 1000 == 0) {     /* verify every 1000 iterations */

      printf ("iter %d n_change %d sse %g t=%.3f s\n",
              iter, n_change, sse, (getmillisecs () - t0) / 1000.0);

    }

  }


  free (clust_assign_new);

  if (clust_assign_out)
    *clust_assign_out = clust_assign;
  else
    free (clust_assign);

  if (med_subset_out)
    *med_subset_out = med_subset;
  else
    free (med_subset);

#undef ALLDISTS

  return sse;
}

float *clustering_kmedoids (int n, int d,
                            const float *points, int k, int nb_iter_max,
                            int **clust_assign_out)
{
  int i;

  /* arbitrarily chose medoids among the n_subset first points */
  int n_subset = 20 * k;
  if (n_subset > n)
    n_subset = n;
  n_subset = n;

  float *alldists = fvec_new (n_subset * n);

  compute_cross_distances (d, n, n_subset, points, points, alldists);

  int *clust_assign;
  int *med_subset;

  double sse =
      clustering_kmedoids_from_dists (n, n_subset, k, nb_iter_max, alldists,
                                      &clust_assign, &med_subset);

  free (alldists);

  float *centroids = fvec_new (k * d);

  for (i = 0; i < k; i++)
    memcpy (centroids + i * d, points + med_subset[i] * d,
            d * sizeof (float));

  {
    double verif_sse = 0;

    for (i = 0; i < n; i++)
      verif_sse +=
          fvec_distance_L2sqr (points + med_subset[clust_assign[i]] * d,
                        points + i * d, d);
    printf ("verif_sse=%g ?= %g\n", verif_sse, sse);
  }

  free (med_subset);

  if (clust_assign_out) {
    *clust_assign_out = clust_assign;
  } else
    free (clust_assign);

  return centroids;
}


/* Algorithm from
@CONFERENCE{trikmeans03,
  author = {Elkan, Charles},
  title = {Using the Triangle Inequality to Accelerate K-Means},
  booktitle = {International Conference on Machine Learning},
  year = {2003},
}
*/

typedef struct {
  int n,d,k;
  const float *points;
  float *centroids;
  int *clust_assign;
  int n_eval;

  /* internal */

  float *c2_dists,*s_dist;
  float *upper,*lower;
  
} trig_t;


static void trig_init(trig_t *trig);

static void trig_update(trig_t *trig,float *new_centroids);

static void trig_compute(trig_t *trig);

static void trig_dealloc(trig_t *trig);





float* clustering_kmeans_assign_trig (int n, int d,
                                      const float *points, int k, int nb_iter_max, 
                                      double normalize, 
                                      int ** clust_assign_out) {
  assert (n > k);

  int i, c, iter = 0;

  float *centroids = fvec_new (k * d);
  float *centroids_new = fvec_new (k * d);

  /* First select a subset of vectors */
  int *randperm = ivec_new_random_perm (n);
  for (i = 0; i < k; i++)
    memcpy (centroids + i * d, points + randperm[i] * d, sizeof (float) * d);
  free (randperm);

  /* vector counting the number of vectors in each centroids */
  int *cluster_size = (int *) malloc (sizeof (int) * n);

  /* to which cluster a vector is assigned ? */
  int *clust_assign = (int *) malloc (sizeof (int) * n);


  /* init public fields */
  trig_t trig={n,d,k,points,centroids,clust_assign,0};

  trig_init(&trig);


  /* a value to check if the centroids have been updated */
  double modiftestvalue, modiftestvalue_new = fvec_sum (centroids, k * d);


  void (*peek_fun) (void *arg, double frac);

  peek_fun = &progress_bar_peek_fun;


  /* launch the construction iterating process */
  do {
    if (normalize != 0)
      floats_normalize2 (centroids, k, d, normalize);

    iter++;
    fprintf (stderr, "Iteration %-3d\n", iter);
    modiftestvalue = modiftestvalue_new;

    memset (centroids_new, 0, sizeof (float) * k * d);
    memset (cluster_size, 0, sizeof (int) * k);
   
    trig_compute(&trig);
    
    for (i = 0; i < n; i++) {
      fvec_add (centroids_new + clust_assign[i] * d, points + i * d, d);
      cluster_size[clust_assign[i]]++;
    }

    /* check here because we want clust_assign to be in sync with centroids */
    if (iter == nb_iter_max)
      break;

    /*
       for (i = 0 ; i < k ; i++) 
       fprintf (stderr, "%d ", cluster_size[i]); 
     */
    fprintf (stderr, "\n");

    for (i = 0; i < k; i++) {
      if (cluster_size[i] > 0) {
        floats_mul (centroids_new + (i) * d, 1.0 / cluster_size[i], d);
      }
      /* Do no allow a cluster to be void. Split another cluster into two */
      else {
        /* Find a non-void cluster */
        do
          c = floor (drand48 () * k);
        while (cluster_size[c] < 4);    /* minimum size for splitting a cluster */

        /* Add a very small gaussian noise to the previous cluster */

        float noise_amount = 0.00001;
        {
          float *ci = centroids_new + i * d;
          float *cc = centroids_new + c * d;
          int j;
          for (j = 0; j < d; j++)
            ci[j] = cc[j] + (drand48 () * 2 - 1) * noise_amount;
        }

        fprintf (stderr, "r");
      }
    }

    /* make the sum to check if the clusters have changed or not */
    modiftestvalue_new = fvec_sum (centroids_new, k * d);

    trig_update(&trig,centroids_new);    
    {
      float *mat_tmp = centroids_new;
      centroids_new = centroids;
      centroids = mat_tmp;
    }
    

  } while (modiftestvalue != modiftestvalue_new);

  fprintf (stderr, "\n");

  {
    int i;
    double sum=0,sum2=0;

    for(i=0;i<k;i++) {
      sum+=cluster_size[i];
      sum2+=cluster_size[i]*(double)cluster_size[i];
    }
    printf("imbalance factor=%.3f,n_eval=%d\n",sum2*k/(sum*sum),trig.n_eval);
    
  }

  trig_dealloc(&trig);
  

  if (normalize)
    floats_normalize2 (centroids, k, d, normalize);

  free (centroids_new);

  if (clust_assign_out)
    *clust_assign_out = clust_assign;
  else
    free (clust_assign);

  free (cluster_size);



  return centroids;
  
  
}

static float Inf=1e30;

static void trig_c2_dists(trig_t *trig) {
  int k=trig->k;
  int i,j;
  compute_cross_distances(trig->d,k,k,trig->centroids,trig->centroids,trig->c2_dists);
  trig->n_eval+=k*(k-1)/2;

  for(i=0;i<k;i++) {
    float closest=Inf;
    for(j=0;j<k;j++) 
      if(trig->c2_dists[j+k*i]<closest && j!=i) 
        closest=trig->c2_dists[j+k*i];
    trig->s_dist[i]=(1/4.0)*closest;
  }
  
}

static void trig_init(trig_t *trig) {
  int k=trig->k;
  int n=trig->n;
  trig->c2_dists=fvec_new(k*k);
  trig->s_dist=fvec_new(k);
  trig->upper=fvec_new(n);
  trig->lower=fvec_new(n*k);

  int i,j;

  /* assume initial assignement of all to centroid 0 */
  for(i=0;i<n;i++) {
    trig->clust_assign[i]=0;
    trig->upper[i]=Inf;
    for(j=0;j<k;j++) 
      trig->lower[i*k+j]=0;
  }  
  
}

static void trig_compute(trig_t *trig) {
  int k=trig->k;
  int n=trig->n;
  int d=trig->d;
  int *c=trig->clust_assign;
  float *u=trig->upper;
  float *s=trig->s_dist;
  int i,j;
  trig_c2_dists(trig);
  
#ifdef TRIG_DB
  float *gt_d2=fvec_new(k*n);
  compute_cross_distances(d,k,n,trig->centroids,trig->points,gt_d2);
#endif

  for(i=0;i<n;i++) {

#ifdef TRIG_DB
    {
      float smallest=Inf;
      for(j=0;j<k;j++) {
        if(gt_d2[j+i*k]<smallest) smallest=gt_d2[j+i*k];
        assert(trig->lower[i*k+j]<=gt_d2[j+i*k]+1e-5);
      }
      assert(u[i]+1e-5>=smallest);
    }
#endif
    
    if(u[i]<=s[c[i]])  /* no chance that there's a better centroid */
      continue;

    float *l=trig->lower+i*k;
    
    /* boolean: u[x] != d(x,c[x]) */
    float d_x_cx_u=u[i]; /* upper bound on distance with current centroid */
    int r_x=1;           /* means d_x_cx_u is != d_x_cx_u */


    for(j=0;j<k;j++) {


      if(j==c[i])     /* our current choice already */
        continue; 

      if(l[j]>=d_x_cx_u)  /* centroid j cannot be better  */
        continue; 

      /* distance of candiate centroid to the current choice */
      float d_cx_c=trig->c2_dists[j+c[i]*k];
     
      if(d_x_cx_u<=(1/4.0)*d_cx_c) /* will not get closer */  
        continue;

      float d_x_cx;

      if(r_x) {
        d_x_cx=d_x_cx_u=fvec_distance_L2sqr(trig->centroids+c[i]*d,trig->points+i*d,d);
        trig->n_eval++;
        r_x=0;
        l[c[i]]=d_x_cx;

        /* we have a tighter bound: check again */

        if(l[j]>=d_x_cx)
          continue; 
      
        if(d_x_cx<=(1/4.0)*d_cx_c)
          continue;

      } else 
        d_x_cx=d_x_cx_u;

      /* now comes the real check.... */
      
      float d_x_c=fvec_distance_L2sqr(trig->centroids+j*d,trig->points+i*d,d);
      trig->n_eval++;
      l[j]=d_x_c;

      if(d_x_c<d_x_cx) {
        c[i]=j;          
        d_x_cx_u=d_x_c;        
      }
      
    }       
    
    u[i]=d_x_cx_u;

#ifdef TRIG_DB
    {
      float smallest=Inf;
      int si=-1;
      for(j=0;j<k;j++) {
        if(gt_d2[j+i*k]<smallest) {smallest=gt_d2[j+i*k]; si=j; }
        assert(l[j]<=gt_d2[j+i*k]+1e-5);
      }
      assert(u[i]+1e-5>=smallest);
      assert(si==c[i]);
    }
#endif
  }

 /*
  for(i=0;i<n;i++) 
    printf("%d ",c[i]);

  printf("\n");
*/

}

static void trig_update(trig_t *trig,float *new_centroids) {
  int i,j;
  int k=trig->k;
  int n=trig->n;
  int d=trig->d;
  float *u=trig->upper;
  float *ds=fvec_new(k);
  int *c=trig->clust_assign;

#ifdef TRIG_DB
  float *gt_d2=fvec_new(k*n);
  compute_cross_distances(d,k,n,new_centroids,trig->points,gt_d2);
#endif

  trig->n_eval+=k;
  for(i=0;i<k;i++ )
    ds[i]=fvec_distance_L2(trig->centroids+i*d,
                           new_centroids+i*d,d);

  for(i=0;i<n;i++) {    
    float d=ds[c[i]];
    u[i]=u[i]+d*d+2*d*sqrt(u[i]);
    
    float *l=trig->lower+i*k;    
    for(j=0;j<k;j++) {
      
      if(l[j]>ds[j]*ds[j]) {
        float new_l=sqrt(l[j])-ds[j];
        l[j]=new_l*new_l;        
      } else 
        l[j]=0;

#ifdef TRIG_DB
      assert(trig->lower[i*k+j]<=gt_d2[j+i*k]+1e-5);      
#endif
    }

#ifdef TRIG_DB
    {
      float smallest=Inf;
      for(j=0;j<k;j++) {
        if(gt_d2[j+i*k]<smallest) smallest=gt_d2[j+i*k];
        assert(trig->lower[i*k+j]<=gt_d2[j+i*k]+1e-5);
      }
      assert(u[i]+1e-5>=smallest);
    }
#endif    

  }
  
  trig->centroids=new_centroids;  
  free(ds);
}

static void trig_dealloc(trig_t *trig) {
  free(trig->lower);
  free(trig->upper);
  free(trig->s_dist);
  free(trig->c2_dists);
}


