
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include <mpi.h>

#include <yael/nn.h>
#include <yael/vector.h>
#include <siftgeo/siftgeo.h>
#include <siftgeo/siftgeo_binarize.h>

/**********************************************************************
 * Misc
 **********************************************************************/

#define NEWA(type,n) (type*)malloc((n)*sizeof(type))
#define NEW(type) NEWA(type,1)

#define __USE_GNU
#include <sched.h>

static int count_cpu() {
  cpu_set_t set;
  sched_getaffinity(0,sizeof(cpu_set_t),&set);
  int i,count=0;
  for(i=0;i<CPU_SETSIZE;i++)
    if(CPU_ISSET(i,&set)) count++;
  return count;
}



static double floats_sum(float *a,int n) {
  double accu=0;
  while(n--) accu+=a[n];
  return accu;
}

static void floats_add(float *a,float *b,int n) {
  while(n--) a[n]+=b[n];
}

static void floats_mul(float *a,float b,int n) {
  while(n--) a[n]*=b;
}

static void floats_normalize2 (float *a,int n,int d, int newnorm) {
  int i, j;
  double vecnorm,f;

  for (i = 0; i < n; i++) {
    float *al=a+d*i;

    /* compute the norm */
    vecnorm = 0;
    for (j = 0; j < d; j++)
      vecnorm += al[j]*al[j];
    
    f=newnorm/sqrt (vecnorm);
    
    /* Normalize the corresponding vector */
    for (j = 0; j < d; j++)
      al[j]*=f;
  }
}


/**********************************************************************
 * Descriptor file random access 
 **********************************************************************/


#define FVEC_SIZE(d) (sizeof(int)+d*sizeof(float))
#define SIFTGEO_SIZE(d) (sizeof (geom_t) + sizeof (int) + (d) * sizeof (char))


FILE *open_ptsfile(const char*fname,int input_is_fvecs,int *dim_out,int*n_out) {
  FILE *f=fopen(fname,"r");

  if(!f) {
    perror("could not open input");
    exit(1);
  }

  struct stat sb;
  fstat (fileno(f), &sb);
  
  size_t fsize=sb.st_size;

  if(fsize==0) {
    fprintf(stderr,"empty point file\n");
    exit(0);
  }
  int d;
  
  int entry_size;
  if(input_is_fvecs) {
    fread(&d,sizeof(int),1,f);
    entry_size=FVEC_SIZE(d);
  } else {
    fseek (f, sizeof (geom_t), SEEK_SET);
    fread (&d, sizeof (int), 1, f);
    entry_size=SIFTGEO_SIZE(d);
  }

  *dim_out=d;

  if(fsize%entry_size!=0) {
    fprintf(stderr,"weird fvecs input size\n");
    exit(0);
  } 
 
  *n_out=fsize/entry_size;

  return f;
}

void load_descs(FILE *infile,int input_is_fvecs,long ofs,
                long n,long d,float *out) {
  int i,j;
  if(input_is_fvecs) {  
    fseek(infile,FVEC_SIZE(d)*ofs,SEEK_SET);
    for(i=0;i<n;i++) {
      int d2;    
      fread(&d2,sizeof(int),1,infile);
      assert(d2==d);    
      fread(out+i*d,sizeof(float),d,infile);    
    }
  } else {
    unsigned char *buf=malloc(d);
    geom_t geom;
    fseek(infile,SIFTGEO_SIZE(d)*ofs,SEEK_SET);
    for (i = 0 ; i < n ; i++) {
      fread (&geom, sizeof (geom_t), 1, infile);
      int d2;
      fread (&d2, sizeof (int), 1, infile);
      assert(d2==d);    
      fread(buf,1,d,infile);
      for(j=0;j<d;j++) 
        out[i*d+j]=buf[j];
    }
  }
  

}


/**********************************************************************
 * MPI root implementation
 **********************************************************************/



#define TAG_SLICE 12012
#define TAG_ASSIGN 12013
#define TAG_PROJ 12014
#define TAG_PROJ2 12015
#define ROOT 0


static void progress_bar(double frac) {
  int n=(int)(frac*70),i;
  char line[72];
  for(i=0;i<n;i++) line[i]='+';
  for(i=n;i<70;i++) line[i]='-';
  line[70]='\n';
  line[71]=0;
  printf(line);
  fflush(stdout);
}

static int verbose=0;


/* assign slices of 0..n-1 to nslave slaves. The slaves return the
   first assignement, which ensures that the load is balanced wrt the
   processing speed. */
static void root_assign_slices(int nslave,int n,int k,int *clust_assign) {
  int i;

  int nslice=16*nslave; /* 16 slices per machine */
  
  int slices[nslice+1];  /* slice boundaries */
  for(i=0;i<=nslice;i++) slices[i]=i*(long)n/nslice;
  int stop[]={-1,-1}; /* send when slave must stop */
  
  /* current requests */
  MPI_Request requests[nslave];
  int n_ret[nslave][2];
  
  /* send request to all slaves */
  int next_slice=0;
  for(i=0;i<nslave;i++) {
    MPI_Send(slices+next_slice++,2,MPI_INT,i+1,TAG_SLICE,MPI_COMM_WORLD);
    MPI_Irecv(n_ret[i],2,MPI_INT,i+1,TAG_SLICE,MPI_COMM_WORLD,&requests[i]);
  }
        
  int n_active=nslave;
  
  while(n_active>0) {
    MPI_Status status[nslave];
    
    progress_bar(next_slice/(float)nslice);
    
    MPI_Waitany(nslave,requests,&i,status);
    
    /* here, i = a slave that returned something */

    if(verbose>1)
      printf("return waitany i=%d %d-%d n_active=%d\n",i,n_ret[i][0],n_ret[i][1],n_active);
    
    /* fill in the corresponding segment in clust_assign */ 
    MPI_Recv(clust_assign+n_ret[i][0],n_ret[i][1]-n_ret[i][0],
             MPI_INT,i+1,TAG_ASSIGN,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

    /* sanity */
    assert(clust_assign[n_ret[i][0]]<k);
    
    if(next_slice<nslice) { /* there are other slices to process */
      if(verbose>1) printf("assign to i=%d slice %d-%d\n",i,slices[next_slice],slices[next_slice+1]);      
      MPI_Send(slices+next_slice++,2,MPI_INT,i+1,TAG_SLICE,MPI_COMM_WORLD);
      MPI_Irecv(n_ret[i],2,MPI_INT,i+1,TAG_SLICE,MPI_COMM_WORLD,&requests[i]);
    } else if(n_active>0) { /* send stop order */
      if(verbose>1) printf("stop i=%d\n",i);     
      MPI_Send(stop,2,MPI_INT,i+1,TAG_SLICE,MPI_COMM_WORLD);            
      n_active--;
    } 

  }

}

/* For the next iterations, the slaves report what slices they take
   care of, followed by the assignment. */
static void root_retrieve_assignement(int nslave,int n,int k,int *clust_assign) {
  int i;

  
  int n_active=nslave;
  
  int n_ret[nslave][2];
  MPI_Request requests[nslave];
  
  for(i=0;i<nslave;i++) 
    MPI_Irecv(n_ret[i],2,MPI_INT,i+1,TAG_SLICE,MPI_COMM_WORLD,&requests[i]);
    
  int ndone=0; /* nb of cells in clust_assign that are filled in */
  while(n_active>0) {
    
    MPI_Status status[nslave];
    if(verbose>1) printf("waitany n_active=%d/%d ndone=%d/%d\n",n_active,nslave,ndone,n);
    MPI_Waitany(nslave,requests,&i,status);

    /* slave i reports... */
    
    if(verbose>1) printf("from waitany i=%d ret=%d %d\n",i,n_ret[i][0],n_ret[i][1]);
    
    if(n_ret[i][0]<0) { /* stop */
      n_active--;            
    } else { /* fill in clust_assign */
      MPI_Recv(clust_assign+n_ret[i][0],n_ret[i][1]-n_ret[i][0],
               MPI_INT,i+1,TAG_ASSIGN,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      ndone+=n_ret[i][1]-n_ret[i][0];
      MPI_Irecv(n_ret[i],2,MPI_INT,i+1,TAG_SLICE,MPI_COMM_WORLD,&requests[i]);      
    }
    
    progress_bar(ndone/(float)n);
  }
  assert(ndone==n || !"not all slices assigned");
      
  if(verbose>1)
    printf("exit root_retrieve_assignement loop\n");
  
}

static float Inf=1.0/0.0;
static float NaN=0.0/0.0;

static int feq(float a,float b) {
  return (fabs(a-b)/(fabs(a)+fabs(b)))<1e-5;
}


/* Principle: here again we get as little information as possible from
   the slaves. The medians are found by bissection. The root only
   manages thresholds.
*/

static siftgeo_binarize_t *  root_compute_he(int n,int k,long d,long nproj,
                                             int he_pca,int *clust_assign,
                                             int *cluster_size) {

  long i,j;

  printf ("* Compute binary signature\n");

  /* Prepare random projection matrix */ 

  siftgeo_binarize_t * sb;

  /* do not compute medians */
  if (he_pca == 0)
    sb = siftgeo_binarize_new (n, d, k, nproj, NULL, NULL);
  else
    assert(!"not implemented");
  
  /* send projection matrix to slaves */
  if(verbose>1) printf("root broadcast projection matrix...\n");
  MPI_Bcast(sb->p,d*nproj,MPI_FLOAT,ROOT,MPI_COMM_WORLD);
  if(verbose>1) printf("...done\n");
  
  /* Retrieve minimas and maximas for each visual word and each
     dimension. Maxima are negated so that we can reduce the whole
     array with a single minimum  */
  
  float *minimaxes=fvec_new(nproj*k*2);
  for(i=0;i<nproj*k*2;i++) minimaxes[i]=Inf;
  if(verbose>1) printf("root retrieve minimaxes...\n");
  MPI_Reduce(MPI_IN_PLACE,minimaxes,nproj*k*2,MPI_FLOAT,MPI_MIN,ROOT,MPI_COMM_WORLD);
  if(verbose>1) printf("... done\n");

  /* Parallel bisection. We initialize the bounds with the minima & maxima */

  long kp_remain=nproj*k;

  typedef struct {
    float t0,t1;
  } bound_t;

  bound_t *bounds=NEWA(bound_t,kp_remain);

  for(i=0;i<kp_remain;i++) {
    bound_t *bound=bounds+i;
    bound->t0=minimaxes[2*i]-0.1;
    bound->t1=-minimaxes[2*i+1]+0.1;
  }

  int *n_below=ivec_new(kp_remain); 
  
  /* map[i] for 0<=i<kp_remain is maps from indices in bounds array
     into the medians array. We remove elements from map when the
     bisection succeeds. Eventually kp_remain drops to 0 */

  int *map=ivec_new(kp_remain);
  for(i=0;i<kp_remain;i++) map[i]=i;    
    
  float *thresholds=fvec_new(nproj*k);
  
  for(;;) {
    
    if(verbose>1) printf("bounds[0]={vals=%g~%g}\n",bounds[0].t0,bounds[0].t1);
    
    for(i=0;i<kp_remain;i++) {
      bound_t *bound=bounds+i;
      thresholds[i]=0.5*(bound->t0+bound->t1);                
    }

    if(verbose>1) printf("root broadcast %ld thresholds...\n",kp_remain);
    MPI_Bcast(thresholds,kp_remain,MPI_FLOAT,ROOT,MPI_COMM_WORLD);
    if(verbose>1) printf("... done\n");

    /* Compact arrays. The slaves will do the same at the other side,
       so the numbering remains in sync */
    int kp_new=0;
    for(i=0;i<kp_remain;i++) if(map[i]>=0) {
      map[kp_new]=map[i];
      bounds[kp_new]=bounds[i];
      thresholds[kp_new]=thresholds[i];
      kp_new++;
    }
    kp_remain=kp_new;
    if(verbose>1) printf("new kp_remain=%ld\n",kp_remain);
    
    if(kp_remain==0) break;

    for(i=0;i<kp_remain;i++) n_below[i]=0;
    if(verbose>1) printf("root getting n_below\n");   
    MPI_Reduce(MPI_IN_PLACE,n_below,kp_remain,MPI_INT,MPI_SUM,ROOT,MPI_COMM_WORLD);    
    if(verbose>1) printf("root synth...\n");
    

    /* Bissection update bounds */
    for(i=0;i<kp_remain;i++) {
      bound_t *bound=bounds+i;

      long p=map[i];
      int cs=cluster_size[p / nproj];

      assert(n_below[i]>=0 && n_below[i]<=cs);

      /*  t0==t1 test useful because there may be multiple occurences
          of values */
      if(n_below[i]==cs/2 || feq(bound->t0,bound->t1)) { 
        sb->medians[p]=thresholds[i];  /* temporary */
        bound->t0=NaN;
        map[i]=-1;
      } else if(n_below[i]<cs/2) {
        bound->t0=thresholds[i];
      } else {
        bound->t1=thresholds[i];
      } 

    }

    if(verbose>1) printf("root ...done\n");

  }

  /* Ask slaves fot the maximum value below the threshold and the
     minimum value above */

  if(verbose>1) printf("root broadcast medians...\n");
  MPI_Bcast(sb->medians,k*nproj,MPI_FLOAT,ROOT,MPI_COMM_WORLD);
  if(verbose>1) printf("...done\n");
  
  for(i=0;i<nproj*k*2;i++) minimaxes[i]=Inf;
  if(verbose>1) printf("root retrieve minimaxes...\n");
  MPI_Reduce(MPI_IN_PLACE,minimaxes,nproj*k*2,MPI_FLOAT,MPI_MIN,ROOT,MPI_COMM_WORLD);
  if(verbose>1) printf("... done\n");

  /* update medians using these */
  for(i=0;i<k;i++) {
    int cs=cluster_size[i];
    
    for(j=0;j<nproj;j++) {
      long p=i*nproj+j;
    
      if(cs==0) {
        sb->medians[p]=0;
      } else if(cs%2==0) {
        sb->medians[p]=(minimaxes[2*p]-minimaxes[2*p+1])*0.5;
      } else 
        sb->medians[p]=minimaxes[2*p];
    }          

  }

  return sb;
}



void root_process(int size,
                  int d,int k,int n,
                  FILE *infile,int input_is_fvecs,
                  int nb_iter_max,double normalize,
                  const char *centroids_fname,
                  const char*clust_assign_fname,
                  const char *he_filename,
                  int nproj,
                  int he_pca) {

  float *centroids=fvec_new(k*d);
  float *centroids_new=fvec_new(k*d);

  int * clust_assign = ivec_new(n);
  int * cluster_size = ivec_new(n);

  int i,iter;

  /* First select a subset of vectors */
  int *randperm=ivec_new_random_perm(n);

  for (i = 0 ; i < k ; i++) 
    load_descs(infile,input_is_fvecs,randperm[i],1,d,centroids+i*d); 

  free(randperm);

  const int nslave=size-1;
  
  
  for(iter=0;iter<nb_iter_max;iter++) { 

    if (normalize != 0)
      floats_normalize2(centroids, k,d,normalize);
       
    /* 1. send centroids to all */
    
    if(verbose>1) printf("sending %d centroids\n",k);
    MPI_Bcast(centroids,k*d,MPI_FLOAT,ROOT,MPI_COMM_WORLD);
    if(verbose>1) printf("sent\n");
        
    fprintf (stderr, "Iteration %d\n", iter);

           
    if(iter==0) { /* assign slices when connections come in */
      
      root_assign_slices(nslave,n,k,clust_assign);

    } else { /* only retrieve slices */

      root_retrieve_assignement(nslave,n,k,clust_assign);
      
    }    

    memset(cluster_size, 0, sizeof (int) * k);
    for (i = 0 ; i < n ; i++) {
      cluster_size[clust_assign[i]]++;
    }    

    /* check here because we want clust_assign to be in sync with centroids */ 
    if(iter==nb_iter_max-1) break;

    /* make new centroids table */
    memset(centroids_new,0,sizeof(float)*k*d);
    if(verbose>1) printf("root reduce...\n");
    MPI_Reduce(MPI_IN_PLACE,centroids_new,k*d,MPI_FLOAT,MPI_SUM,ROOT,MPI_COMM_WORLD);
    if(verbose>1) printf("...done\n");

    if(verbose)
      for (i = 0 ; i < k ; i++)
        fprintf (stderr, "%d ", cluster_size[i]);
    fprintf (stderr, "\n");
    
    /* recompute centroids */
    for (i = 0 ; i < k ; i++) {
      if (cluster_size[i] > 0) {
        floats_mul (centroids_new+i*d, 1.0/cluster_size[i], d);
      } else { /* Do no allow a cluster to be void. Split a cluster into two */
        int c;
	/* Find a non-void cluster */
	do
	  c = floor (drand48 () * k);
	while (cluster_size[c] < 4); /* minimum size for splitting a cluster */

	/* Add a very small gaussian noise to the previous cluster */
		
        float noise_amount=0.00001;
        {
          float *ci=centroids_new + i * d;
          float *cc=centroids_new + c * d;
          int j;
          for(j=0;j<d;j++) 
            ci[j]=cc[j]+(drand48()*2-1)*noise_amount;
        }     	
	
	fprintf (stderr, "r");
      }
    }

   
    /* swap centroid tables */

    {
      float* mat_tmp = centroids_new;
      centroids_new = centroids;
      centroids = mat_tmp;
    }

  }  

  fvecs_write (centroids_fname,d,k,centroids);

  if(clust_assign_fname) {
    ivecs_write(clust_assign_fname, n, 1, clust_assign);
  }

  if(he_filename) {
    
    siftgeo_binarize_t *sb=root_compute_he(n,k,d,nproj,he_pca,clust_assign,cluster_size);

    siftgeo_binarize_write (he_filename, sb);
    
  }


  free(centroids_new);
  free(cluster_size);
  free(clust_assign);
    

}


/**********************************************************************
 * MPI slave implementation
 **********************************************************************/

typedef struct slice_t { 
  int n[2];                    /* begin, end */
  float *pts;                  /* point coordinates */
  float *pts_projected;        /* by HE projection matrix */
  int *assignement;            /* closest centroid */
  struct slice_t *next;        /* linked list */
} slice_t; 


static slice_t *slave_get_new_slice(int rank,int k,int d,FILE *infile,int input_is_fvecs) {

  int ns[2];
  
  if(verbose>2) printf("rank %d recv\n",rank);
  MPI_Recv(ns,2,MPI_INT,ROOT,TAG_SLICE,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
  if(verbose>2) printf("rank %d got %d-%d\n",rank,ns[0],ns[1]);
  
  if(ns[0]<0) return NULL;
        
  int nsl=ns[1]-ns[0];
  
  slice_t *slice=NEW(slice_t);
  slice->n[0]=ns[0];
  slice->n[1]=ns[1];
  slice->assignement=NEWA(int,nsl);
  slice->pts=fvec_new(nsl*d);
  
  if(verbose>2) printf("rank %d loading...\n",rank);
  load_descs(infile,input_is_fvecs,slice->n[0],nsl,d,slice->pts);
  if(verbose>2) printf("rank %d ...done\n",rank);

  return slice;
}


static void slave_slice_assign_and_send(int rank,int k,long d,
                                        slice_t *slice,
                                        const float *centroids,
                                        float *centroids_new) {
  int i;

  int nsl=slice->n[1]-slice->n[0];
  
  if(verbose>2) printf("rank %d quantize %d %d %ld \n",rank,nsl,k,d);
  quantize_codebook_thread(nsl,k,d,centroids,slice->pts,slice->assignement,count_cpu(),
                           NULL,NULL);    
  if(verbose>2) printf("rank %d quantize done\n",rank);

  /* accumulate */
  for(i=0;i<nsl;i++) 
    floats_add(centroids_new+slice->assignement[i]*d,slice->pts+i*d,d);

  if(verbose>2)
    printf("rank %d send %d-%d\n",rank,slice->n[0],slice->n[1]);
  
  MPI_Send(slice->n,2,MPI_INT,ROOT,TAG_SLICE,MPI_COMM_WORLD);
  assert(slice->assignement[0]<k);
  MPI_Send(slice->assignement,nsl,MPI_INT,ROOT,TAG_ASSIGN,MPI_COMM_WORLD);

}


static void slave_bounds_for_he(int rank,slice_t *slices,long d,long k,long nproj) {
  long i,j;
  
  /* Retrieve projection matrix  */
  float *projmat=fvec_new(nproj*d);
  if(verbose>2) printf("rank %d receive proj matrix...\n",rank);
  MPI_Bcast(projmat,nproj*d,MPI_FLOAT,ROOT,MPI_COMM_WORLD);
  if(verbose>2) printf("rank %d ... done\n",rank);
  
  /* Project slices and compute bounding boxes of each point in
     projected space. The maxima are coded as -the value to reduce
     with a single MPI_MIN */
  
  float *minimaxes=fvec_new(k*nproj*2);
  for(i=0;i<k*nproj*2;i++) minimaxes[i]=Inf;

  slice_t *slice;
  for(slice=slices;slice;slice=slice->next) {
    int nsl=slice->n[1]-slice->n[0];

    if(verbose>2) 
      printf("rank %d projecting slice %d-%d \n",rank,slice->n[0],slice->n[1]);

    siftgeo_binarize_t sb={ /* fake */
      d,nproj,k,projmat,NULL
    };
    
    slice->pts_projected=he_transform_points(&sb, slice->pts,nsl);
    
    for(i=0;i<nsl;i++) {
      int vw=slice->assignement[i];
      float *vals=slice->pts_projected+nproj*i;
      float *mm=minimaxes+nproj*vw*2;
      for(j=0;j<nproj;j++) {
        if(vals[j]<mm[2*j]) mm[2*j]=vals[j];
        if(-vals[j]<mm[2*j+1]) mm[2*j+1]=-vals[j];
      }      
    }
    
  }

  if(verbose>2) printf("rank %d sending minimaxes...\n",rank);  
  float *minimaxes_aux=fvec_new(k*nproj*2);
  MPI_Reduce(minimaxes,minimaxes_aux,2*nproj*k,MPI_FLOAT,MPI_MIN,ROOT,MPI_COMM_WORLD);
  if(verbose>2) printf("rank %d ... done\n",rank);
  

  /* now repeately get thresholds for a decreasing subset of the
     points and count the number of values below the threshold. */

  long kp_remain=k*nproj;
  float *thresholds=fvec_new(kp_remain);
  int *map=ivec_new(kp_remain); /* maps initial indices to subset indices */
  for(i=0;i<kp_remain;i++) map[i]=i;

  int *n_below=ivec_new(kp_remain),*n_below_aux=ivec_new(kp_remain);
  for(;;) {


    if(verbose>2) printf("rank %d getting thresholds, kp_remain=%ld...\n",rank,kp_remain);    
    MPI_Bcast(thresholds,kp_remain,MPI_FLOAT,ROOT,MPI_COMM_WORLD);
    if(verbose>2) printf("rank %d ... done\n",rank);


    /* switch to new subset */
    long ii=0,old_kp_remain=kp_remain;
    kp_remain=0;
    for(i=0;i<k*nproj;i++) if(map[i]>=0) {
      if(!isnan(thresholds[ii])) {
        map[i]=kp_remain;
        thresholds[kp_remain]=thresholds[ii];
        kp_remain++;
      } else {
        map[i]=-1;
      }
      ii++;
    }        
    assert(ii==old_kp_remain);
    if(verbose>2) printf("rank %d new kp_remain=%ld\n",rank,kp_remain);
   

    /* stop if subset becomes empty */
    if(kp_remain==0) break;

    
    /* count values below the threshold in each slice */
    for(i=0;i<kp_remain;i++) n_below[i]=0;  
    for(slice=slices;slice;slice=slice->next) {
      int nsl=slice->n[1]-slice->n[0];
       
      for(i=0;i<nsl;i++) {
        int vw=slice->assignement[i];
        float *vals=slice->pts_projected+nproj*i;
        
        for(j=0;j<nproj;j++) {
          int p=map[vw*nproj+j];
          if(p<0) continue;
                    
          if(vals[j]<thresholds[p]) n_below[p]++;
          
        }
      }
    }

    if(verbose>2) printf("rank %d sending n_below\n",rank);
    MPI_Reduce(n_below,n_below_aux,kp_remain,MPI_INT,MPI_SUM,ROOT,MPI_COMM_WORLD);    
    if(verbose>2) printf("rank %d loop\n",rank);
    
  }


  /* Last stage: root needs to know real values of the points. For
     each threshold, we return the smallest valeue above and the
     largest value below (negated for reduce) */

  if(verbose>2) printf("rank %d getting thresholds\n",rank);  
  MPI_Bcast(thresholds,k*nproj,MPI_FLOAT,ROOT,MPI_COMM_WORLD);
  if(verbose>2) printf("rank %d ... done\n",rank);
    
  for(i=0;i<2*k*nproj;i++) minimaxes[i]=Inf;
  for(slice=slices;slice;slice=slice->next) {
    int nsl=slice->n[1]-slice->n[0];
       
    for(i=0;i<nsl;i++) {
      int vw=slice->assignement[i];
      float *vals=slice->pts_projected+nproj*i;
        
      for(j=0;j<nproj;j++) {
        long p=vw*nproj+j;
        
        if(vals[j]>thresholds[p]) {
          if(vals[j]<minimaxes[2*p]) minimaxes[2*p]=vals[j];
        } else {
          if(-vals[j]<minimaxes[2*p+1]) minimaxes[2*p+1]=-vals[j];
        }
      }
    }
  }

  if(verbose>2) printf("rank %d sending minimaxes\n",rank);  
  MPI_Reduce(minimaxes,minimaxes_aux,2*k*nproj,MPI_FLOAT,MPI_MIN,ROOT,MPI_COMM_WORLD);
  if(verbose>2) printf("rank %d done\n",rank);
  
}


void slave_process(int rank,
                   int d,int k,int n,
                   int nb_iter_max,
                   FILE *infile,int input_is_fvecs,
                   int nproj) {
  
  float *centroids=fvec_new(k*d);
  float *centroids_new=fvec_new(k*d);
  
  int iter,i;
  
  slice_t *slices=NULL;

  for(iter=0;iter<nb_iter_max;iter++) { 
    
    /* receive centroids */    
    MPI_Bcast(centroids,k*d,MPI_FLOAT,ROOT,MPI_COMM_WORLD);
    
    
    if(verbose>2) printf("rank %d after bcast\n",rank);    
    
    memset(centroids_new,0,sizeof(float)*k*d);

    if(iter==0) {
      
      slice_t *slice=NULL;

      for(;;) {

        slice_t *new=slave_get_new_slice(rank,k,d,infile,input_is_fvecs);
        if(!new) break;

        new->next=slice;
        slice=new;
      
        slave_slice_assign_and_send(rank,k,d,slice,centroids,centroids_new);        

      }
      slices=slice;

    } else {

      slice_t *slice;

      for(slice=slices;slice;slice=slice->next)        
        slave_slice_assign_and_send(rank,k,d,slice,centroids,centroids_new);        

      int stop[]={-1,-1};
      MPI_Send(stop,2,MPI_INT,ROOT,TAG_SLICE,MPI_COMM_WORLD);     

    }

    /* check here because we want clust_assign to be in sync with centroids */ 
    if(iter==nb_iter_max-1) break;

    /* centroids used as a temp buffer */
    if(verbose>2) printf("rank %d reduce...\n",rank);
    MPI_Reduce(centroids_new,centroids,k*d,MPI_FLOAT,MPI_SUM,ROOT,MPI_COMM_WORLD);
    if(verbose>2) printf("... done\n");
        
  }
 
  if(nproj>0) 
    slave_bounds_for_he(rank,slices,d,k,nproj);


}
                   
                   

/**********************************************************************
 * main
 **********************************************************************/



void usage() {
  
  fprintf(stderr,"mpirun [mpirun options] cluster_mpi -i datafile.fvecs -o outfile.fvecs [-k nb_centroids] [-niter #] [-nrlz] [-v [-v]]\n");
  exit(1);

}



int main(int argc,char **argv) {
  
  /* no buffering of stdout */
  setvbuf(stdout,NULL,_IONBF,0);  

  MPI_Init(&argc, &argv);

  int rank,size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  int i;
  int nb_iter_max=40;
  int normalize=0;
  int input_is_fvecs=0;
  const char *infn=NULL;
  const char *outfn=NULL;
  const char *clust_assign_fname=NULL;
  int k=1000;
  int seed;
  const char *he_filename=NULL;
  int nproj=64;
  int he_pca=0;  
  
  for(i=1;i<argc;i++) {
    char *a=argv[i];

    if(!strcmp(a,"-h") || !strcmp(a,"--help")) usage ();
    else if(!strcmp(a,"-i") && i+1<argc) infn=argv[++i];
    else if(!strcmp(a,"-o") && i+1<argc) outfn=argv[++i];
    else if(!strcmp(a,"-niter") && i+1<argc && sscanf(argv[i+1],"%d",&nb_iter_max)==1) i++;
    else if(!strcmp(a,"-nrlz")) normalize=1;
    else if(!strcmp(a,"-infvecs")) input_is_fvecs=1;
    else if(!strcmp(a,"-k") && i+1<argc && sscanf(argv[i+1],"%d",&k)==1) i++;    
    else if(!strcmp(a,"-seed") && i+1<argc && sscanf(argv[i+1],"%d",&seed)==1) {i++; srandom(seed); }
    else if(!strcmp(a,"-v")) verbose++;
    else if(!strcmp(a,"-he") && i+1<argc) he_filename=argv[++i];
    else if(!strcmp(a,"-nproj") && i+1<argc && sscanf(argv[i+1],"%d",&nproj)==1) i++;
    else if(!strcmp(a,"-he_pca")) he_pca=1;
    else if(!strcmp(a,"-ca") && i+1<argc) clust_assign_fname=argv[++i];
    else {
      fprintf(stderr,"unknown arg %s\n",a);      
      usage();
    }
  }
  
  int d,n;

  FILE *infile=open_ptsfile(infn,input_is_fvecs,&d,&n);
  
  if(verbose) 
    printf("rank=%d/%d k=%d n=%d d=%d\n",rank,size,k,n,d);

  if(rank==0) {

    printf("will store result in %s\n",outfn);

    root_process(size,
                 d,k,n,infile,input_is_fvecs,nb_iter_max,normalize,
                 outfn,
                 clust_assign_fname,
                 he_filename,nproj,he_pca);

        
  } else {
    
    slave_process(rank,d,k,n,nb_iter_max,infile,input_is_fvecs,
                  he_filename ? nproj : -1);

  
  }
  fclose(infile); 

  MPI_Finalize();

  return 0;
}
