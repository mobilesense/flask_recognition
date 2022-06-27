#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include <yael/vector.h>
#include <yael/matrix.h>
#include <yael/nn.h>
#include <yael/sorting.h>


#include <siftgeo/siftgeo.h>
#include <siftgeo/siftgeo_binarize.h>




/* blas subroutine for real matrix multiply */

#define real float
#define integer int

int sgemm_ (char *transa, char *transb, integer * m, integer *
            n, integer * k, real * alpha, const real * a, integer * lda,
            const real * b, integer * ldb, real * beta, real * c__,
            integer * ldc);

#undef real
#undef integer



siftgeo_binarize_t * siftgeo_binarize_new (int n, int d, int nbvw, int nproj, 
					   float * points, int * clust_assign)
{
  int i, j, l, npt, vw;
  siftgeo_binarize_t * sb = malloc (sizeof (siftgeo_binarize_t));
  sb->d = d;
  sb->nproj = nproj;
  sb->nbvw = nbvw;
  sb->p = fvec_new(d*nproj);
  sb->medians = fvec_new(nproj*nbvw);  

  /* Generate the random projection matrix */

  if(nproj>d) {
    printf("siftgeo_binarize_new: warn nproj=%d > d=%d: will use several orthogonal matrices\n",
           nproj,d);
  }

  for(l=0;l<(nproj+d-1)/d;l++) {

    float *m = random_orthogonal_basis (d);
    
    for (i = l*d ; i<(l+1)*d && i < nproj ; i++)
      for (j = 0 ; j < d ; j++)
        sb->p[i * d + j] = m[(i-l*d) * d + j];
    
    free (m);

  }

  /* means medians will be handled externally */
  if(!points) return sb;


  /* Get the number of elements per cluster */
  int * cluster_size = ivec_new_histogram (nbvw,clust_assign,n);

  /* Compute the median values */
  for (vw = 0 ; vw < nbvw ; vw++) {
    npt = cluster_size[vw];

    if (npt == 0) {
      fprintf (stderr, "# Warning: no point in cluster %d\n",vw);
      continue;
    }
    
    /* select only the points that correspond to the correct visual word */
    float * pointsvw = fvec_new (d * npt);
    for (i = 0, j = 0 ; i < n ; i++) 
      if (clust_assign[i] == vw) {
	memcpy (pointsvw + j * d, points + i * d, d * sizeof (float));
	j++;
      }

    siftgeo_binarize_fill_medians (sb, vw, pointsvw, npt);
    free (pointsvw);
  }

  free (cluster_size);
  return sb;
}



/* generate a projection matrix using the principal subspace of difference 
   between the vectors and their centroids */
siftgeo_binarize_t * siftgeo_binarize_pca_new (int n, int d, int nbvw, 
					       int nproj, int whiten,
					       float * points, 
					       int * clust_assign, 
					       float * centroids)
{
  int i, j, npt, vw;
  siftgeo_binarize_t * sb = malloc (sizeof (siftgeo_binarize_t));
  sb->d = d;
  sb->nproj = nproj;
  sb->nbvw = nbvw;
  sb->p = fvec_new (d * nproj);
  sb->medians = fvec_new (nproj * nbvw);  

  assert(nproj<=d);

  /* first generate the set of vectors of differences between local descriptors
     and their centroids */
  float * centroid_offset = fvec_new (n * d);
  for (i = 0 ; i < n ; i++)
    for (j = 0 ; j < d ; j++)
      centroid_offset[i * d + j] = points[i * d + j] 
	- centroids[clust_assign[i] * d + j];

  /* compute the PCA on the set of input vectors */
  float *pca_vec=compute_pca (n, d, centroid_offset);
  free (centroid_offset);

  if(whiten) {

    /* Generate the random projection matrix */
    float *randrot = random_orthogonal_basis (nproj);
  
    /* The final projection matrix is obtained by whitening the PCA projection
       matrix with a random (gaussian drawn) rotation matrix */ 
    
    { /* perform m := randrot * pca_vec */
      float one=1,zero=0;
      sgemm_("N","N",&d,&nproj,&nproj,&one,pca_vec,&d,randrot,&nproj,&zero,sb->p,&d);
    }
    free (randrot);
  } else {
    memcpy(sb->p,pca_vec,sizeof(float)*d*nproj);
  }


  {
    printf("############ verif\n");
    int i,j,k;

    double sum2=0;
    for(i=0;i<n;i++) {
      float *v=points+d*i;
      float *vbase=centroids+d*clust_assign[i];

      sum2+=fvec_distance_L2sqr(v,vbase,nproj);

    }

    printf("input sum of dists: %g\n",sum2);

    sum2=0;
    float *vt=fvec_new(nproj);
    for(i=0;i<n;i++) {
      float *v=points+d*i;
      float *vbase=centroids+d*clust_assign[i];
      
      for(j=0;j<nproj;j++) {
        vt[j]=0;
        for(k=0;k<d;k++) 
          vt[j]+=sb->p[k+j*d]*(v[k]-vbase[k]); 
      }
      
      sum2+=fvec_norm2sqr(vt,nproj);
 
    }    
    free(vt);

    printf("output sum of dists: %g\n",sum2);
    

    printf("############ end verif\n");
  }



  free (pca_vec);

  /* means medians will be handled externally */
  if(!points) return sb;

  /* Get the number of elements per cluster */
  int * cluster_size = ivec_new_histogram (nbvw,clust_assign,n);

  /* Compute the median values */
  for (vw = 0 ; vw < nbvw ; vw++) {
    npt = cluster_size[vw];

    if (npt == 0) {
      fprintf (stderr, "# Warning: no point in cluster\n");
      continue;
    }
    
    /* select only the points that correspond to the correct visual word */
    float * pointsvw = (float *) fvec_new (d * npt);
    for (i = 0, j = 0 ; i < n ; i++) 
      if (clust_assign[i] == vw) {
	memcpy (pointsvw + j * d, points + i * d, d * sizeof (float));
	j++;
      }

    siftgeo_binarize_fill_medians (sb, vw, pointsvw, npt);
    free (pointsvw);
  }

  return sb;
}


void siftgeo_binarize_crop_nproj(siftgeo_binarize_t *sb,int new_nproj) {
  int i,j;
  
  assert(new_nproj<=sb->nproj);

  /* compress medians array */

  for(i=0;i<sb->nbvw;i++) 
    for(j=0;j<new_nproj;j++) 
      sb->medians[i*new_nproj+j]=sb->medians[i*sb->nproj+j];

  sb->nproj=new_nproj;

}


/* Compute the binary signatures of a set of descriptors

   P           the matrix used for the projection of descriptors
   nbproj      number of distinct random projection = number of bits in the signature
   siftgeo     the descriptors
   nbdes       the number of descriptors
*/
binsign_t * siftgeo_binarize_full (int nbproj,int nbvw,float *p,
                                   point_t * siftgeo, point_t * vwgeo, 
				   int nbdes, float *medians)
{
  int i, j, des;

  float *proj = calloc(nbdes*nbproj,sizeof(*proj));
  
  if(nbdes==0) return NULL;
  int d=siftgeo[0].dim;

  /* Perform the projections */
  for (des = 0 ; des < nbdes ; des++) {
    
    for (i = 0 ; i < nbproj ; i++)
      for (j = 0 ; j < d ; j++)
	proj[des*nbproj+i] += p[i+j*d] * siftgeo[des].desc[j];
  }
    
  /* Construction of the binary signatures */
  binsign_t * bin = (binsign_t *) malloc (nbdes * sizeof (binsign_t));
  memset (bin, nbdes * sizeof (binsign_t), 0);

  for (des = 0 ; des < nbdes ; des++) {

    for (i = 0 ; i < nbproj ; i++) {
      int vw = vwgeo[des].vw;

      if (vw >= nbvw || vw < 0)   /* When vw is invalid -> TO BE MODIFIED */
	vw = 0;   

      binsign_t bit = (proj[des*nbproj+i] < medians[vw*nbproj+i] ? 0 : 1);
      bin[des] += bit << i;
    }
  }

  free (proj);

  return bin;     
}


/* read the projection matrix and the parameters from a file */
siftgeo_binarize_t * read_sb_common (FILE *f) 
{
  siftgeo_binarize_t *sb = malloc(sizeof(siftgeo_binarize_t));

  int buf[3];
  if(fread (buf, sizeof(int), 3, f) != 3) {
    fprintf (stderr, "read_sb_common error 1\n");
    return NULL;
  }
  /* fwrite(f,[d,randproj,nbvw],"int32"); 
   * nbvw is 0 for only projection files */

  int d = sb->d = buf[0];
  int nproj = sb->nproj = buf[1];
  int nbvw = sb->nbvw = buf[2];

  sb->p = fvec_new (d * nproj);
  if (fread (sb->p, sizeof (float), d * nproj, f) != d * nproj) {
    fprintf (stderr, "read_sb_common error 2\n");
    return NULL;
  }
  sb->medians = fvec_new (nproj *nbvw);  
  if (fread (sb->medians, sizeof(float), nproj * nbvw, f) != nproj * nbvw) {
    fprintf (stderr, "siftgeo_binarize_read error 3\n");
    return NULL;
  }  

  return sb;
}


void write_sb_common (FILE *f, const siftgeo_binarize_t *sb) 
{
  int d = sb->d;
  int nproj = sb->nproj;
  int nbvw = sb->nbvw;
  int buf[3] = {d, nproj, nbvw};

  if (fwrite (buf, sizeof(int), 3, f) != 3) {
    fprintf (stderr, "siftgeo_binarize_write error 1\n");
    return;
  }

  if (fwrite (sb->p, sizeof (float), d * nproj, f) != d * nproj) {
    fprintf (stderr, "siftgeo_binarize_read error 2\n");
    return;
  }

  if (fwrite (sb->medians, sizeof (float), nproj * nbvw, f) != nproj * nbvw) {
    fprintf (stderr, "siftgeo_binarize_read error 3\n");
    return;
  }  
}


/* read params */
siftgeo_binarize_t * siftgeo_binarize_read (const char*fname) 
{
  FILE *f = fopen(fname,"r");
  if (!f) {
    perror ("siftgeo_binarize_read");
    return NULL;
  }

  siftgeo_binarize_t * sb = read_sb_common(f);

  fclose (f);
  return sb;
}


void siftgeo_binarize_write (const char*fname, const siftgeo_binarize_t *sb) 
{
  FILE *f = fopen(fname, "w");
  if(!f) {
    perror ("siftgeo_binarize_write");
    return;
  }

  write_sb_common (f, sb);

  fclose (f);
}


/* scalar product with a unsigned char vector */
static float dotprod (float *a, unsigned char *b, int n) 
{
  double accu = 0;
  while (n--) accu += a[n] * b[n];
  return accu;  
}


/* compute binary signature */
void siftgeo_binarize (siftgeo_binarize_t *sb,pointset_t *siftgeos,pointset_t *vwgeos) 
{
  int i,j;

  assert(siftgeos->n==vwgeos->n);

  int n=siftgeos->n;
  int nproj=sb->nproj;
  int d=sb->d;

  for(i=0;i<n;i++) {
    point_t *siftgeo = &siftgeos->pts[i];
    point_t *vwgeo = &vwgeos->pts[i];
    unsigned char *desc = siftgeo->desc;
    binsign_t b = 0;
    float *med = sb->medians+vwgeo->vw*nproj;
    
    for (j = 0 ; j < nproj ; j++) {
      float v = dotprod (sb->p + j * d, desc, d);
      if (v >= med[j]) 
        b |= 1LL<<j;      
    }
    
    vwgeo->binsign = b;
  }

}


float * he_transform_points(siftgeo_binarize_t *sb,const float *points,int npt) 
{
  int nproj = sb->nproj;
  int d = sb->d;
  float *ret = fvec_new (npt * nproj);
  float one = 1.0, zero = 0.0;

  sgemm_("Transposed", "Not trans", &nproj, &npt, &d, &one, sb->p, &d, 
	 points, &d, &zero, ret, &nproj);

  return ret;
}

/* ditto, but return a transposed matrix (useful to compute medians) */
static float * he_transform_points_transp (siftgeo_binarize_t *sb, float *points, int npt) 
{
  float *f = he_transform_points (sb, points, npt);
  int nproj = sb->nproj;
  float *ret = fvec_new (npt * nproj);
  int i, j;
  for (i = 0 ; i < nproj ; i++)
    for (j = 0 ; j < npt ; j++) 
      ret[j + npt * i] = f[i + nproj * j];
  free (f);
  return ret;
}


void siftgeo_binarize_ffq (siftgeo_binarize_t *sb, const float *points,
			   pointset_t *vwgeos, int k) 
{
  int i, j;

  int n = vwgeos->n / k;
  assert (vwgeos->n % k ==0 );
  int nproj = sb->nproj;

  float *coords = he_transform_points (sb, points, n);

  for (i = 0 ; i < n ; i++) {
    float *cl = coords + i * nproj;
    int l;
    for (l = 0 ; l < k ; l++) {

      point_t *vwgeo = &vwgeos->pts[i * k + l];
      
      if(vwgeo->vw<0) continue; /* will be reomoved later */

      float *med = sb->medians + vwgeo->vw * nproj;
      
      binsign_t b = 0;
      for (j = 0 ; j < nproj ; j++) {
        if (cl[j] >= med[j]) 
          b |= 1LL << j;      
      }    
      vwgeo->binsign = b;
    }
  }

  free (coords);
}


void siftgeo_binarize_ffq_table (siftgeo_binarize_t *sb,
                                int n,int k,
                                 const float *points,const int *assign,
                                unsigned long long *binsigns)
{
  int i, j;

  int nproj = sb->nproj;

  float *coords = he_transform_points (sb, points, n);
  int nl=(nproj+63)/64;  

  for (i = 0 ; i < n ; i++) {
    float *cl = coords + i * nproj;
    int l;
    for (l = 0 ; l < k ; l++) {
      unsigned long long * bs=binsigns+(i*k+l)*nl;
      int vw=assign[i * k + l];
      if(vw<=0) continue; /* not assigned, should be handled later */
      float *med = sb->medians + vw * nproj;

      int j0;
      for (j0 = 0 ; j0 < nproj ; j0+=64) {
        binsign_t b = 0;
        for (j = 0 ; j < 64 ; j++) {
          if (cl[j+j0] >= med[j+j0])
            b |= 1LL << j;
        }
        bs[j0/64] = b;
      }

    }
  }

  free (coords);

}




/* compute medians from training set of points */
void siftgeo_binarize_fill_medians (siftgeo_binarize_t *sb, int vw, float *points, int npt) 
{
  int i;
  int nproj = sb->nproj;

  /* vector of medians to fill in */  
  float *med = sb->medians + vw * nproj;

  if (npt ==0 ) {
    fprintf (stderr, "warning, no points to estimate vw %d\n", vw);
    memset (med, 0, sizeof (float) * nproj);
    return;
  }

  float *coords = he_transform_points_transp (sb, points, npt);

  if (npt == 1) 
    memcpy (med, coords, sizeof (float) * nproj);

  else {    
    for (i = 0 ; i < nproj ; i++) {
      /* c[0..npt-1] = coordinates of all points along projection i */
      float *c = coords + i * npt;

      med[i]=fvec_median(c,npt);

    }
  }
  free (coords);
}

void siftgeo_binarize_display (siftgeo_binarize_t *bs) 
{
  int i, j;
  printf ("siftgeo_binarize_t {\n");
  printf ("  d=%d nproj=%d nbvw=%d\n", bs->d, bs->nproj, bs->nbvw);
  printf ("  p=[");
  for (i = 0 ; i < bs->nproj ; i++) {
    for (j = 0 ; j < bs->d ; j++) 
      printf ("%10g ", bs->p[j + i * bs->d]);    
    printf (";\n     ");
  }
  printf ("  ]\n");
  printf ("  medians=[");
  for (i = 0 ; i < bs->nbvw ; i++) {
    for (j = 0 ; j < bs->nproj ; j++) 
      printf ("%10g ", bs->medians[j + i *bs->nproj]);
    printf (";\n     ");
  }
  printf ("  ]\n");
  printf ("}\n");
}



void siftgeo_binarize_delete(siftgeo_binarize_t *bs) 
{
  free (bs->p);
  free (bs->medians);
  free (bs);
}


/* display a signature */
void siftgeo_binarize_binsign_display (binsign_t bs)
{
  int k;
  for (k = 0 ; k < sizeof (bs) * 8 ; k++)
    printf ("%d", (int) ((bs >> k) & 1));
  printf (" ");
}






