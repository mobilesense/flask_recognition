#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <float.h>
#include <sys/time.h>

#include <yael/machinedeps.h>
#include <yael/vector.h>
#include <yael/clustering.h>
#include <yael/nn.h>
#include <utils/hkm.h>

#include <siftgeo/siftgeo.h>
#include <siftgeo/siftgeo_and_nn.h>


#define INFMT_FVEC             0
#define INFMT_SIFTGEO          1

#define CLUST_TYPE_KMEANS      0
#define CLUST_TYPE_HKM         1
#define CLUST_TYPE_KMEDOIDS    2

#include "siftgeo_binarize.h"


void display_help (const char * progname)
{
    fprintf (stderr, "Usage: %s  -k # -i filename [-o clustfile] [-nrlz #] [-nredo #] [-seed seed]\n"
	     "                  [-infvecs] [-he paramsfile -nproj # [-he_pca]] [-maxn maxn]\n"
             "                  [-uf] [-hkm -bf # -nlevel]\n\n", 
	     progname);
    
    fprintf (stderr, 
             "Input\n"
             "    -i filename     input file containing the set of vectors to cluster\n"
             "    -infvecs        use fvecs instead of siftgeo (default)\n"
             "    -maxn #         maximum number of input points to use\n"
             "\nClustering parameters\n"                  
             "    -k #            number of centroids to produce\n"
             "    -nredo #        number of runs to perform\n"
             "    -nrlz #         perform normalization of vectors while iterating\n"
             "    -he_pca         use the PCA version of Hamming Embedding\n"
             "    -he_pca_nowhite idem without whitening\n"
             "    -nproj #        nb projections used in the HE (default=64)\n"
             "    -hkm            use a hierachical k-means (HKM)\n"
             "    -bf             branching factor of the HKM\n"
             "    -nlevel         number of tree levels of the HKM\n"
             "\nOuput\n"                          
             "    -uf             display the imbalance factor\n"
             "    -o clustfile    output cluster file (fvec format)\n"
             "    -he paramsfile  output parameter file for Hamming Embedding\n");
    exit (1);
}



/* compute the unbalanced factor associated with a given visual vocabulary. 
   Input: histogram of occurences of visual words (k distinct visual words) */
double compute_uf (int * h, int k)
{
  int i, n = 0;     /* n gather the total number of occurences */
  double uf = 0;

  for (i = 0 ; i < k ; i++) 
    n += h[i];


  for (i = 0 ; i < k ; i++) 
    uf += ((double) h[i] / n) * ((double) h[i] / n);

  return uf * k;
}



int main (int argc, char **argv)
{
  int i;

  char * filesource_name = NULL;
  /* The file where the clusters are written */
  char * fileout_name = NULL;
  int k = 1000;
  int nredo=1;
  int max_iter = 30;
  int maxn=10*1000*1000;
  int nrlz = 0;
  int infmt=INFMT_SIFTGEO;
  int clust_type = CLUST_TYPE_KMEANS;
  int seed_arg=-1;
  int report_uf=0;

  /* for HE */
  char *he_filename=NULL; /* none by default */
  int nproj = 64;
  int he_pca=0;

  /* for HKM */
  int nlevel = -1;
  int bf = -1;

  if (argc == 1)
    display_help (argv[0]);

  for(i=1;i<argc;i++) {
    char *a=argv[i];

    if(!strcmp(a,"-h") || !strcmp(a,"--help")) {
      display_help (argv[0]);
    } else if(!strcmp(a,"-i") && i+1<argc) {
      filesource_name=argv[++i];
    } else if(!strcmp(a,"-o") && i+1<argc) {
      fileout_name=argv[++i];
    } else if(!strcmp(a,"-k") && i+1<argc && sscanf(argv[i+1],"%d",&k)==1) {
      i++;
    } else if(!strcmp(a,"-nredo") && i+1<argc && sscanf(argv[i+1],"%d",&nredo)==1) {
      i++;
    } else if(!strcmp(a,"-max_iter") && i+1<argc && sscanf(argv[i+1],"%d",&max_iter)==1) {
      i++;
    } else if(!strcmp(a,"-maxn") && i+1<argc && sscanf(argv[i+1],"%d",&maxn)==1) {
      i++;
    } else if(!strcmp(a,"-nrlz")) {
      nrlz = 1;
    } else if(!strcmp(a,"-uf")) {
      report_uf=1;
    } else if(!strcmp(a,"-medoids")) {
      clust_type=CLUST_TYPE_KMEDOIDS;
    } else if(!strcmp(a,"-infvecs")) {
      infmt=INFMT_FVEC;
    } else if(!strcmp(a,"-hkm")) {
      clust_type = CLUST_TYPE_HKM;
    } else if(!strcmp(a,"-seed") && i+1<argc && sscanf(argv[i+1],"%d",&seed_arg)==1) {
      i++;
    } else if(!strcmp(a,"-he") && i+1<argc) {
      he_filename=argv[++i];
    } else if(!strcmp(a,"-nproj") && i+1<argc && sscanf(argv[i+1],"%d",&nproj)==1) {
      i++;
    } else if(!strcmp(a,"-he_pca")) {
      he_pca=1;
    } else if(!strcmp(a,"-he_pca_nowhite")) {
      he_pca=2;
    } else if(!strcmp(a,"-nlevel") && i+1<argc && sscanf(argv[i+1],"%d",&nlevel)==1) {
      i++;
    } else if(!strcmp(a,"-bf") && i+1<argc && sscanf(argv[i+1],"%d",&bf)==1) {
      i++;
    } else {
      fprintf(stderr,"could not parse argument %s\n",a);
      display_help (argv[0]);
    }
  }

  

  /* optionally, initialize the random generator */
  if (seed_arg>=0) {
    srandom (seed_arg);
  }


  int n, d;

  float *points;

  double score_new, score_min = DBL_MAX;


  switch (infmt) {
  case INFMT_FVEC:
    {
      n = fvecs_new_read (filesource_name, &d, &points);
      if (n < 0) {
	fprintf (stderr, "# Unable to load vectors from file %s\n", filesource_name);
	exit (2);
      }
      else
	fprintf (stderr, "* Found %d distinct vectors\n", n);
    }
    break;

  case INFMT_SIFTGEO:
    { 
      point_t *corners;
      n=read_points (filesource_name,&corners,0);
      if(n<0) return 1;
      points=siftgeo_to_fvecs(corners,n,&d);
      delete_points(corners,n);
    }
    break;

  default:
    fprintf (stderr, "# Invalid input format\n");
    exit (3);
  }
  
  fprintf (stdout, "[k-means] found %d vectors, dim %d -> ", n,d);
  if (n > maxn)
    n = maxn;
  fprintf (stdout, "use %d\n", n);

  while ( 0 <nredo--) {
      fprintf(stdout, "%d remaining runs\n",nredo);
      fflush(stdout);

      int * clust_assign = NULL;
      int * clust_size = NULL;

      /* Standard k-means and Hamming Embedding */
      float * centroids = NULL;
      hkm_t * hkm = NULL;

      if (clust_type == CLUST_TYPE_KMEANS) {

        printf("max_iter=%d\n",max_iter);

        float * centroids = clustering_kmeans_assign_with_score (n, d, points, k, 
                                                                 max_iter, nrlz, count_cpu(), &score_new, 
                                                                 &clust_assign);
        
        if (score_new<score_min) {
          fprintf(stdout,"best clustering (togo=%d,new=%g,previous=%g), saving...\n",nredo,score_new,score_min);
          
          score_min = score_new;
          
          clust_size=ivec_new_histogram(k,clust_assign,n);
          
          fvecs_write (fileout_name, d, k, centroids);
        }

      } 
      else if (clust_type == CLUST_TYPE_KMEDOIDS) {

	  float * centroids = clustering_kmedoids (n, d, points, k, 
		  max_iter,  
		  &clust_assign);

	  clust_size=ivec_new_histogram(k,clust_assign,n);

	  fvecs_write (fileout_name, d, k, centroids);

      } 
      else if (clust_type == CLUST_TYPE_HKM) {

	  /* Required information was not entered */
	  if (bf == -1 || nlevel == -1)
	      display_help (argv[0]);

	  fprintf (stdout, "* hierarchical clustering: bf=%d, lev=%d, n=%d, d=%d\n", 
		   bf, nlevel, n, d);

	  hkm = hkm_learn (n, d, nlevel, bf, points, max_iter, 0, &clust_assign); 
	  fprintf (stdout, "* hierarchical clustering: k=%d, bf=%d, lev=%d\n", 
		  hkm->k, hkm->bf, hkm->nlevel);

	  k = hkm->k;
	  clust_size = ivec_new_histogram (hkm->k, clust_assign, n);
	  //	  ivec_print (clust_size, k);

	  hkm_write (fileout_name, hkm);
      }



      /*-------------------------------------------------------------*/
      /* optionally, create a set of parameter for Hamming Embedding */
      if (he_filename) {
	  fprintf (stdout, "* Compute binary signature\n");

	  siftgeo_binarize_t * sb;

	  /* standard matrix projection or use principal subspace ? */
	  if (he_pca == 0)
	      sb = siftgeo_binarize_new (n, d, k, nproj, points, 
                                         clust_assign);
	  else
            sb = siftgeo_binarize_pca_new (n, d, k, nproj, he_pca==1,
                                           points, clust_assign, centroids);

	  /* Write the binary structure to disk */
	  siftgeo_binarize_write (he_filename, sb);
      }


      /* Optionally, display the unbalance factor */
      //  ivec_print (clust_size, k);

      if (report_uf)
	  fprintf (stdout, "uf=%.3f\n", compute_uf (clust_size, k));

      free (centroids);
      free (clust_assign);
      free (clust_size);

  }

  free (points);
  return 0;
}

