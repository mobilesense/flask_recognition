#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include <yael/clustering.h>
#include <yael/vector.h>
#include <yael/nn.h>
#include <yael/machinedeps.h>
#include <yael/sorting.h>

#include <utils/annvw.h>

#include <siftgeo/siftgeo.h>
#include <siftgeo/siftgeo_and_nn.h>




#define SLICE_SIZE (8192*8*8) 
/* #define SLICE_SIZE (512) */


#include "siftgeo_binarize.h"


void display_help (const char * progname)
{
  fprintf (stderr, "Usage: %s options\n\n", progname);
  fprintf (stderr, 
           "Learns an approximate cluster assignement\n"
           "\n"
           "Input\n"
           "-clust clustfile         the cluster file for which we want ANN search\n"
           "-i learn.siftgeo         the siftgeo input file used to learn the connectivity between layers\n"
           "[-inassign ca.ivec]      cluster assignement of siftgeo file (if known)\n"
           "\n"
           "Parameters\n"
           "[-ann_k nbvw_coarse]     number of centroids on the intermediate level\n"
           "[-ann_thres thres]       keep only edges that are followed by at least this many learning points\n"
           "[-nodefaultedges]        don't include edges from fine clusters to the coarse centroids they belong to\n"
           "\n"
           "Output\n"
           "-o clustfile2            the file of coarse clusters\n"
           "-edges edgesfile         store the 2stage approximate assignment connectivity in edgesfile\n"
           "[-outassign ca.ivec]     cluster assignement of siftgeo file\n");
  exit (1);
}


void imf(int *cluster_size,int k) {
  double sum=0,sum2=0;
  int i;
  for(i=0;i<k;i++) {
    sum+=cluster_size[i];
    sum2+=cluster_size[i]*(double)cluster_size[i];
  }
  fprintf(stderr,"imbalance factor=%.3f\n",sum2*k/(sum*sum));
}

int main (int argc, char **argv)
{
  int i, j, d;

  char * fclust1_name = NULL;
  char * fedges_name = NULL;
  
  char *inassign=NULL,*outassign=NULL;
  
  int max_iter = 20;
  long nbvw2 = -1;

  char * fsiftgeo_name = "/dev/stdin";

  /* The file where the clusters are written. */
  char * fclustfile_name = NULL;

  int thres = 1;

  int defaultedges=1;
 
  for(i=1;i<argc;i++) {
    char *a=argv[i];

    if(!strcmp(a,"-h") || !strcmp(a,"--help")) {
      display_help (argv[0]);
    } else if(!strcmp(a,"-clust") && i+1<argc) {
      fclust1_name=argv[++i];
    } else if(!strcmp(a,"-edges") && i+1<argc) {
      fedges_name=argv[++i];
    } else if(!strcmp(a,"-inassign") && i+1<argc) {
      inassign=argv[++i];
    } else if(!strcmp(a,"-outassign") && i+1<argc) {
      outassign=argv[++i];
    } else if(!strcmp(a,"-max_iter") && i+1<argc && sscanf(argv[i+1],"%d",&max_iter)==1) {
      i++;
    } else if(!strcmp(a,"-ann_k") && i+1<argc && sscanf(argv[i+1],"%ld",&nbvw2)==1) {
      i++;
    } else if(!strcmp(a,"-i") && i+1<argc) {
      fsiftgeo_name=argv[++i];
    } else if(!strcmp(a,"-o") && i+1<argc) {
      fclustfile_name=argv[++i];
    } else if(!strcmp(a,"-ann_thres") && i+1<argc && sscanf(argv[i+1],"%d",&thres)==1) {
      i++;
    } else if(!strcmp(a,"-nodefaultedges")) {
      defaultedges=0;
    } else {
      fprintf(stderr,"could not parse argument %s\n",a);
      display_help (argv[0]);
    }

  }

  assert (thres >= 1);
  fprintf (stderr, "* Threshold on connections: %d\n", thres);

  if(!fclust1_name || !fedges_name || nbvw2<0 || !fclustfile_name) {
    fprintf(stderr,"Mandatory argument missing\n");    
    display_help (argv[0]);
  }

  /* Load the cluster file */

  float * centroids1;
  long nbvw1 = fvecs_new_read (fclust1_name, &d, &centroids1);
  assert (nbvw1 > 0);
  fprintf (stderr, "* Found %ld clusters in loaded cluster file\n", nbvw1);

  assert (nbvw1 > nbvw2);

  /* Read the learning set from disk */
  fprintf (stderr, "* Read the siftgeo learning set for ANN\n");
  int npt = count_points (fsiftgeo_name, 0);
  if (npt < 0) {
    fprintf (stderr, "Error while reading the learning siftgeo file %s\n", fsiftgeo_name);
    exit (1);
  }
  fprintf (stderr, "* Found %d points in the learning siftgeo set\n", npt);

 
  /* Compute the coarse clusters using the second level as input and write it down to disk */
  fprintf (stderr, "* Produce the coarse clustering\n");
  
  srandom(0);

  int *clust_assign0;
  ann_vw_t *ann=ann_vw_new(d,nbvw1,centroids1,nbvw2,max_iter,&clust_assign0);
    
  /* Assign the points to the two levels of clusters */
  int *clust_assign1 = NULL;
  int *clust_assign2 = malloc(sizeof(int)*npt);

  if(inassign) {
    printf("loading initial assignement from %s\n",inassign);
    int npt2;
    clust_assign1=ivec_new_read(inassign,&npt2);
    assert(npt2==npt);
  } else 
    clust_assign1 = malloc(sizeof(int)*npt);

  FILE *fsiftgeo=fopen(fsiftgeo_name,"r");
  
  for (i = 0 ; i < npt ; i += SLICE_SIZE) {
    
    int slice_npt=i+SLICE_SIZE<npt ? SLICE_SIZE : npt-i;

    fprintf (stderr, "slice %d+%d / %d\n", i,slice_npt,npt);
    
    pointset_t *ps=pointset_alloc(slice_npt);
    
    for(j=0;j<slice_npt;j++) {
      int ret=read_point_t(fsiftgeo,ps->pts+j,0);
      assert(ret==1 && ps->pts[j].dim==d);
    }
    
    float * points = siftgeo_to_fvecs (ps->pts, slice_npt, &d);
 
    if(!inassign) {
      fprintf (stdout, "* Quantization of the learning set -> fine grain quantizer\n");
      ann_vw_learn_assign(ann,1,points,slice_npt,clust_assign1 + i);
    }

    fprintf (stdout, "* Quantization of the learning set -> coarse quantizer\n");
    ann_vw_learn_assign(ann,2,points,slice_npt,clust_assign2 + i);

    free (points);

    pointset_delete(ps);
  }
  fclose(fsiftgeo);
  
/*
  {
    imf(ivec_new_histogram(nbvw1,clust_assign1,npt),nbvw1);
    imf(ivec_new_histogram(nbvw2,clust_assign2,npt),nbvw2);    
  }
*/
  
  if(outassign) {
    printf("storing fine assignement in %s\n",outassign);
    ivecs_write (outassign, npt, 1, clust_assign1);    
  }


  printf("making edges\n");  
  ann_vw_make_edges(ann,clust_assign1,clust_assign2,npt,thres);

  if(defaultedges) {
    printf("including default edges\n");
    int n_added=0;
    for(i=0;i<nbvw1;i++) {
      j=clust_assign0[i];      
      assert(j<nbvw2);
      int n=ann->n_edge[j];
      int idx=ivec_sorted_find(ann->edges[j],n,i);
      if(idx>=0 && ann->edges[j][idx]==i) continue; /* already included */
      
      /* shift array right */
      ann->n_edge[j]++;
      int *e=ann->edges[j]=realloc(ann->edges[j],sizeof(int)*(n+1));
      memmove(e+idx+2,e+idx+1,sizeof(int)*(n-idx-1));
      e[idx+1]=i;

      n_added++;
    }    
    printf("   added %d edges\n",n_added);
  }
  

/*
  for(i=0;i<nbvw1;i++) {
    j=clust_assign2[i];
    assert(ivec_sorted_count_occurrences(ann->edges[j],ann->n_edge[j],clust_assign1[i])==1);
  }
*/
  printf("stats:\n");
  ann_vw_describe(ann);  

  printf("storing result to %s %s\n",fclustfile_name,fedges_name);  
  ann_vw_write(ann,fclustfile_name,fedges_name);

  return 0;
}


