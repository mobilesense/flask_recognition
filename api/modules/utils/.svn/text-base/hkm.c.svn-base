#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include <yael/vector.h>
#include <yael/nn.h>
#include <yael/clustering.h>
#include <yael/sorting.h>

#include "hkm.h"



/*--------------------------------------------------------------
 * balancing clusters       
 --------------------------------------------------------------*/

void quantize_with_balfact (const float * centroids, const float * v, 
			     const float * balfact, int k, int n, int d,
			     int * idx_out)
{
  int i, j;
  float * dis = (float *) malloc (n * k * sizeof (*dis));
  float * dismod = (float *) malloc (k * sizeof (*dismod));
  compute_cross_distances (d, k, n, centroids, v, dis);
  
  /* compute and sort the distances */
  for (i = 0 ; i < n ; i++) {
    for (j = 0 ; j < k ; j++) 
      dismod[j] = dis[i * k + j] * balfact[j];

    /* find the nearest neighbors and update the number of points assigned to a cell */
    idx_out[i] = fvec_arg_min (dismod, k);
  }

  free (dis);
  free (dismod);
}


void quantize_with_balfact_new (const float * centroids, const float * v, 
			     const float * balfact, int k, int n, int d,
			     int * idx_out)
{
  int i, j;
  float * dis = (float *) malloc (n * k * sizeof (*dis));
  float * dismod = (float *) malloc (k * sizeof (*dismod));
  compute_cross_distances (d, k, n, centroids, v, dis);
  
  /* compute and sort the distances */
  for (i = 0 ; i < n ; i++) {
    for (j = 0 ; j < k ; j++) 
      dismod[j] = dis[i * k + j] + balfact[j];

    /* find the nearest neighbors and update the number of points assigned to a cell */
    idx_out[i] = fvec_arg_min (dismod, k);
  }

  free (dis);
  free (dismod);
}


void compute_balance_factors (const float * centroids, const float * v, 
			      int k, int n, int d, float * balfact)
{
  int i, j, loop;
  float target_nnpercell = n / (float) k;

  float * dis = (float *) malloc (n * k * sizeof (*dis));
  float * dismod = (float *) malloc (k * sizeof (*dismod));
  int * nnpercell = (int *) malloc (k * sizeof (*nnpercell));
  float * celldissqr = (float *) malloc (k * sizeof (*celldissqr));

  for (j = 0 ; j < k ; j++)
    balfact[j] = 1;

  compute_cross_distances (d, k, n, centroids, v, dis);


  for (loop = 0 ; loop < 60 ; loop++) {

    memset (nnpercell, 0, k * sizeof (*nnpercell));
    memset (celldissqr, 0, k * sizeof (*nnpercell));

    /* compute and sort the distances */
    for (i = 0 ; i < n ; i++) {
      for (j = 0 ; j < k ; j++) 
	dismod[j] = dis[i * k + j] * balfact[j];
	
      /* find the nearest neighbors and update the number of points assigned to a cell */
      int best = fvec_arg_min (dismod, k);
      nnpercell[best]++;
    }

    ivec_print (nnpercell, k);

    for (j = 0 ; j < k ; j++) {
      if (nnpercell[j] == 0)
	nnpercell[j] = 1;
      balfact[j] *= pow (nnpercell[j] / target_nnpercell, 0.5 / sqrt (d));
    }
  }

  free (dismod);
  free (dis);
  free (nnpercell);
  free (celldissqr);
}


void compute_balance_factors_new (const float * centroids, const float * v, 
			      int k, int n, int d, float * balfact)
{
  int i, j, loop;
  float target_nnpercell = n / (float) k;

  float * dis = (float *) malloc (n * k * sizeof (*dis));
  float * dismod = (float *) malloc (k * sizeof (*dismod));
  int * nnpercell = (int *) malloc (k * sizeof (*nnpercell));
  float * celldissqr = (float *) malloc (k * sizeof (*celldissqr));
  float meandissqr = 0;


  memset (balfact, 0, k * sizeof (*balfact));
  compute_cross_distances (d, k, n, centroids, v, dis);


  for (loop = 0 ; loop < 100 ; loop++) {

    memset (nnpercell, 0, k * sizeof (*nnpercell));
    memset (celldissqr, 0, k * sizeof (*nnpercell));

    /* compute and sort the distances */
    for (i = 0 ; i < n ; i++) {
      for (j = 0 ; j < k ; j++) 
	dismod[j] = dis[i * k + j] + balfact[j];
	
      /* find the nearest neighbors and update the number of points assigned to a cell */
      int best = fvec_arg_min (dismod, k);
      nnpercell[best]++;
      celldissqr[best] += dismod[best];
      //      printf ("%g ", dismod[best]);
      meandissqr += dismod[best] / n;
    }

    for (j = 0 ; j < k ; j++) 
      celldissqr[j] /= nnpercell[j];

    printf ("---------------------------\n");
    ivec_print (nnpercell, k);
    fvec_print (celldissqr, k);
    fvec_print (balfact, k);

    for (j = 0 ; j < k ; j++) {
      if (nnpercell[j] == 0)
	nnpercell[j] = 1;

      balfact[j] += celldissqr[j] * (pow (nnpercell[j] / target_nnpercell, 0.2 / sqrt (d)) - 1);
    }
  }

  free (dismod);
  free (dis);
  free (nnpercell);
  free (celldissqr);
}


/*--------------------------------------------------------------
 * hierarchical clustering
 --------------------------------------------------------------*/

hkm_t * hkm_new (int d, int nlevel, int bf)
{
  int l, j, k = 1;
  hkm_t *hkm = (hkm_t *) malloc (sizeof (*hkm));
  hkm->nlevel = nlevel;
  hkm->bf = bf;
  hkm->d = d;
  hkm->centroids = (float **) malloc (nlevel * sizeof (*hkm->centroids));
  hkm->balfact = (float **) malloc (nlevel * sizeof (*hkm->balfact));

  for (l = 0 ; l < nlevel ; l++) {
    hkm->centroids[l] = fvec_new (k * bf * d);
    hkm->balfact[l] = fvec_new (k * bf);
    for (j = 0 ; j < k * bf ; j++)
      hkm->balfact[l][j] = 1;
    k *= bf;
  }
  hkm->k = k;
  return hkm;
}


hkm_t *hkm_learn (int n, int d, int nlevel, int bf,
		  float *points, int nb_iter_max, 
		  int balmod, int **clust_assign_out)
{
  int i, l, parent, k = 1;

  hkm_t *hkm = hkm_new (d, nlevel, bf);

  /* the absolute assignement of all points and the sizes of clusters */
  int *node_assign = calloc (sizeof (int), n);

  /* the buffer that receives the vectors gathered by parent node */
  float *v = fvec_new (n * d);

  /* Initialization */
  for (l = 0; l < nlevel; l++) {

    /* sort the vectors depending on which cluster they have been assigned to,
       and compute the number of vectors assigned to each cluster 
       *** NOTE: to replace with the find_k_max function of ivfgeo
       -> put this function in a separate library             */
    int *node_assign_idx = malloc (sizeof (*node_assign_idx) * n);
    ivec_sort_index (node_assign, n, node_assign_idx);
    //    ivec_print (node_assign_idx, n);

    /* Re-order the vectors depending on the previous order */
    for (i = 0; i < n ; i++)
      memmove (v + d * i, points + d * node_assign_idx[i], sizeof (*points) * d);

    /* k is the number of nodes/leaves at this level */
    int pos = 0;
    for (parent = 0; parent < k ; parent++) {

      fprintf (stderr, "\n[Level %d | Parent %d]\n\n", l, parent);

      /* Count the number of vectors assigned to this internal node */
      int nassign = 0;
      while (pos + nassign < n)
        if (node_assign[node_assign_idx[pos + nassign]] == parent)
          nassign++;
        else
          break;

       fprintf (stderr, "parent=%d | nassign=%d | pos=%d\n",
               parent, nassign, pos);

      if (nassign == 0) {
        fprintf (stderr, "# Problem: no enough vectors in a node\n");
        exit (1);
      }

      /* Perform the clustering on this subset of points */
      int *clust_assign;
      float *centroids = clustering_kmeans_assign (nassign, d,
                                                   v + d * pos, bf,
                                                   nb_iter_max, 0,
                                                   &clust_assign);

      if (balmod == 1) {
	/* Compute the balancing factors and re-assign the vectors to centroids */
	compute_balance_factors (centroids, v + d * pos, bf,
				 nassign, d, hkm->balfact[l] + parent * bf);

	quantize_with_balfact (centroids, v + d * pos,
			       hkm->balfact[l] + parent * bf, bf,
			       nassign, d, clust_assign);
      }


      memcpy (hkm->centroids[l] + d * parent * bf, centroids,
              d * bf * sizeof (*centroids));

      /* Update the assignment for those points */
      for (i = 0; i < nassign; i++) {
        int truepos = node_assign_idx[pos + i];
        node_assign[truepos] = node_assign[truepos] * bf + clust_assign[i];
      }

      free (centroids);
      free (clust_assign);
      pos += nassign;
    }

    k *= bf;
    free (node_assign_idx);
  }

  if(clust_assign_out) {
    *clust_assign_out = (int *) malloc (n * sizeof (int));
    memcpy (*clust_assign_out, node_assign, n * sizeof (int));
  } 
  free (node_assign);
  free (v);

  return hkm;
}


void hkm_delete (hkm_t * hkm)
{
  int l;
  for (l = 0; l < hkm->nlevel; l++) {
    free (hkm->centroids[l]);
    free (hkm->balfact[l]);
  }

  free (hkm->centroids);
  free (hkm->balfact);
  free (hkm);
}



/* Quantization usign the hierarchical clustering */
int * hkm_quantize (const hkm_t * hkm, const float * v, int npt)
{
  int i, l, vw, vwtmp;
  int * vidx = ivec_new (npt);

  int nlevel = hkm->nlevel;
  int bf = hkm->bf;
  int d = hkm->d;

  for (i = 0 ; i < npt ; i++) {
    vw = 0;
    for (l = 0 ; l < nlevel ; l++) {
      /* at this point, vw contains the parent node */
      //      quantize_codebook (1, bf, d, hkm->centroids[l] + vw * d * bf,
      //			 v + d * i, &vwtmp, NULL, NULL);

      quantize_with_balfact (hkm->centroids[l] + d * vw * bf, v + d * i, 
			     hkm->balfact[l] + vw * bf, bf, 1, d, &vwtmp);
      vw = vw * bf + vwtmp;
    }
    vidx[i] = vw;
  }
  return vidx;
}



/* retrieve the centroids from a particular level */
float * hkm_get_centroids (const hkm_t * hkm, int l, int no)
{
  return hkm->centroids[l] + hkm->d * hkm->bf * no;
}



/***********************************************************************/
/* I/O function for hkm                                                */

/* Macros to handle the i/O of the ivfgeo structure */
#define HKM_READ_ERROR(ret, expected_ret)				\
  {									\
    if (ret != (expected_ret)) {					\
      fprintf (stderr, "# Unable to read the hkm file %s\n",  filename);\
      free (hkm);							\
      return NULL;							\
    }									\
  }

#define HKM_WRITE_ERROR(ret, expected_ret)				\
  {									\
    if (ret != (expected_ret)) {					\
      fprintf (stderr, "# Unable to write the hkm file %s\n", filename);\
      return;								\
    }									\
  }


void hkm_write (const char *filename, const hkm_t * hkm)
{
  int ret = 0, l, k = hkm->bf;
  FILE *f = fopen (filename, "w");

  fprintf (stderr, "### no operational at the time\n");
  assert (0);

  ret = fwrite (&hkm->nlevel, sizeof (hkm->nlevel), 1, f);
  HKM_WRITE_ERROR (ret, 1);
  ret = fwrite (&hkm->bf, sizeof (hkm->bf), 1, f);
  HKM_WRITE_ERROR (ret, 1);
  ret = fwrite (&hkm->k, sizeof (hkm->k), 1, f);
  HKM_WRITE_ERROR (ret, 1);
  ret = fwrite (&hkm->d, sizeof (hkm->d), 1, f);
  HKM_WRITE_ERROR (ret, 1);

  for (l = 0; l < hkm->nlevel; l++) {
    ret = fvec_fwrite (f, hkm->centroids[l], k * hkm->d);
    k *= hkm->bf;
  }
}


hkm_t *hkm_read (const char *filename)
{
  int ret = 0, l;
  FILE *f = fopen (filename, "w");

  hkm_t *hkm = (hkm_t *) malloc (sizeof (*hkm));

  fprintf (stderr, "### no operational at the time\n");
  assert (0);

  assert (hkm);

  ret = fread (&hkm->nlevel, sizeof (hkm->nlevel), 1, f);
  HKM_READ_ERROR (ret, 1);
  ret = fread (&hkm->bf, sizeof (hkm->bf), 1, f);
  HKM_READ_ERROR (ret, 1);
  ret = fread (&hkm->k, sizeof (hkm->k), 1, f);
  HKM_READ_ERROR (ret, 1);
  ret = fread (&hkm->d, sizeof (hkm->d), 1, f);
  HKM_READ_ERROR (ret, 1);

  hkm->centroids = (float **) malloc (hkm->nlevel * sizeof (*hkm->centroids));
  assert (hkm->centroids);

  int k = hkm->bf;
  for (l = 0; l < hkm->nlevel; l++) {
    /* need to allocate the memory */
    hkm->centroids[l] = (float *) malloc (k * sizeof (**hkm->centroids));
    assert (hkm->centroids[l]);

    ret = fvec_fread (f, hkm->centroids[l]);
    HKM_READ_ERROR (ret, k * hkm->d);

    k *= hkm->bf;
  }

  return hkm;
}
