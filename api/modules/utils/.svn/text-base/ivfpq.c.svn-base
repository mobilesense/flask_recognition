#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include <sys/time.h>


#include <yael/vector.h>
#include <yael/nn.h>
#include <yael/sorting.h>
#include <yael/machinedeps.h>

#include <utils/generic_hash_table.h>

#include "ivfpq.h"

#define NEWA(type,n) (type*)malloc(sizeof(type)*(n))
#define NEWAC(type,n) (type*)calloc(n,sizeof(type))
#define NEW(type) NEWA(type,1)


ivfpq_t *ivfpq_new (int nbvw)
{
  ivfpq_t *ivf = NEWAC (ivfpq_t,1);

  ivf->nbvw = nbvw;
  ivf->pq = NULL;

  ivf->labels = NEWAC (int *, nbvw);
  ivf->codes = NEWAC (pqcode_t *, nbvw);

  ivf->nbelems = NEWAC (int, nbvw);
  ivf->segsize = NEWAC (int, nbvw);

  ivf->flags = IVFPQ_FLAG_ASYM;
  ivf->centroids_codes = NULL;
  ivf->pq_assign = NULL;

  return ivf;
}



void ivfpq_delete (ivfpq_t * ivf)
{
  int i;
  if(ivf->pq) 
    pq_free (ivf->pq);

  for (i = 0; i < ivf->nbvw; i++) {
    free (ivf->labels[i]);
    pqcode_free (ivf->codes[i]);
  }

  free (ivf->labels);
  free (ivf->codes);
  free (ivf->nbelems);
  free (ivf->segsize);
  free (ivf->centroids);

  free (ivf);
}

void ivfpq_display(const ivfpq_t *ivf) {
  int i,j,k;
  printf("ivfpq={\n  pq=");
  pq_display(ivf->pq,0);
  
  printf("  invlists[%d]={\n",ivf->nbvw);
  int codelen=pqcode_sizeof(ivf->pq);
  
  for(i=0;i<ivf->nbvw;i++) {
    printf("    (%d entries)={",ivf->nbelems[i]);
    
    const unsigned char *codes=(void*)ivf->codes[i];    
    const int *labels=ivf->labels[i];    
    for(j=0;j<ivf->nbelems[i];j++) {
      
      printf(" %d,",labels[j]);
      for(k=0;k<codelen;k++) 
        printf("%02x",*codes++);
    }

    printf("}\n");
  }

  printf("  }\n"
         "}\n");

}

/************************************************** building */

void ivfpq_add_vw (ivfpq_t * ivf, const int *vw, const int *labels,
                   const float *val, int n)
{
  int i;
  pq_t *pq = ivf->pq;

  for (i = 0; i < n; i++) {
    int w = vw[i];

    if(w==-1) continue; /* means descriptor contains NaNs */

    assert (0 <= w && w < ivf->nbvw);

    int ne = ivf->nbelems[w]++;

    if (ne >= ivf->segsize[w]) {
      int ss = ivf->segsize[w] = ivf->segsize[w]<8 ? 8 : ivf->segsize[w] * 3 / 2;
      ivf->labels[w] = realloc (ivf->labels[w], sizeof (*ivf->labels) * ss);
      ivf->codes[w] = pqcode_realloc (pq, ivf->codes[w], ss);
    }

    ivf->labels[w][ne] = labels[i];
    pq_encode (pq, val + pq->d * i, 1, pqcode_offset (pq, ivf->codes[w], ne));

  }

}



ivfpq_t *ivfpq_dup(ivfpq_t *ivf) {
  ivfpq_t *ivf2=ivfpq_new(ivf->nbvw);
  ivf2->pq=ivf->pq;
  ivf2->centroids=ivf->centroids;
  ivf2->flags=ivf->flags;
  return ivf2;
}



void ivfpq_merge(ivfpq_t *ivf,ivfpq_t *ivf2) {
  int w;
  assert(ivf->pq == ivf2->pq);
  pq_t *pq=ivf->pq;

  for(w=0;w<ivf->nbvw;w++) {
    
    int ne=ivf->nbelems[w];
    
    ivf->nbelems[w]+=ivf2->nbelems[w];

    if(ivf->nbelems[w]>ivf->segsize[w]) {
      int ss=ivf->segsize[w]=ivf->nbelems[w];
      ivf->labels[w]=realloc(ivf->labels[w],sizeof(*ivf->labels)*ss);
      ivf->codes[w]=pqcode_realloc(pq,ivf->codes[w],ss);
    } 
    
    memcpy(ivf->labels[w]+ne,ivf2->labels[w],sizeof(*ivf->labels)*ivf2->nbelems[w]);
    free(ivf2->labels[w]);
    ivf2->labels[w]=NULL;

    memcpy(pqcode_offset (pq,ivf->codes[w],ne),ivf2->codes[w],ivf2->nbelems[w]*pqcode_sizeof(ivf->pq));
    free(ivf2->codes[w]);
    ivf2->codes[w]=NULL;
  }
  
  ivf2->pq=NULL;
  ivf2->centroids=NULL;
  ivfpq_delete(ivf2);  

}

/************************************************** querying */

ivfpq_query_stats_t ivfpq_query_stats={0,};

double getmillisecs() {
  struct timeval tv;
  gettimeofday(&tv,NULL);
  return tv.tv_sec*1e3 +tv.tv_usec*1e-3;
}



int ivfpq_query_vw (const ivfpq_t * ivf, const int *vw, int n,
                    const float *val, int k, int multiple_val, int *labels,
                    float *dists)
{

  maxheap_t *mh = maxheap_new (k);


  int i;
  for (i = 0; i < n; i++) {
    int w = vw[i];
    assert (0 <= w && w < ivf->nbvw);

    int nd = ivf->nbelems[w];

    ivfpq_query_stats.nvisited+=nd;

    float *dists = NEWA (float, nd);

    if (ivf->flags & IVFPQ_FLAG_ASYM)
      pq_asym_L2sq (ivf->pq, val, ivf->codes[w], nd, dists);
    else
      assert (!"not implemented");

    maxheap_add_multiple_labels (mh, ivf->labels[w], dists, nd);

    free (dists);

    if (multiple_val)
      val += ivf->pq->d;
  }


  maxheap_sort (mh);

  int nout = mh->i;

  for (i = 0; i < nout; i++) {
    labels[i] = mh->elts[i].label;
    dists[i] = mh->elts[i].val;
  }

  maxheap_delete (mh);

  return nout;
}

int ivfpq_query_vw_get_label_map (const ivfpq_t * ivf, const int *vw, int n,
                                  const float *val, int k, int multiple_val, int *labels,
                                  float *dists,int *label_map)
{

  maxheap_t *mh = maxheap_new (k);

  int i;
  int nscan=0;
  for (i = 0; i < n; i++) {
    int w = vw[i];
    assert (0 <= w && w < ivf->nbvw);

    int nd = ivf->nbelems[w];

    ivfpq_query_stats.nvisited+=nd;

    float *dists = NEWA (float, nd);

    if (ivf->flags & IVFPQ_FLAG_ASYM)
      pq_asym_L2sq (ivf->pq, val, ivf->codes[w], nd, dists);
    else
      assert (!"not implemented");

    maxheap_add_multiple (mh, nscan, nd, dists);
    nscan+=nd;

    free (dists);

    if (multiple_val)
      val += ivf->pq->d;
  }

  maxheap_sort_labels (mh);

  int nout = mh->i;
  nscan=0;
  int j=0;
  for (i = 0; i < n; i++) {
    int w = vw[i];    
    int nd = ivf->nbelems[w];
    
    while(j<nout && mh->elts[j].label<nscan+nd) {      
      labels[j] = ivf->labels[w][mh->elts[j].label-nscan];
      dists[j] = mh->elts[j].val;
      j++;
    }
    nscan+=nd;
    label_map[i]=j;
  }

  maxheap_delete (mh);

  return nout;
}



int ivfpq_count_nbelems (const ivfpq_t * ivf)
{
  int tot = 0, vw;

  for (vw = 0; vw < ivf->nbvw; vw++)
    tot += ivf->nbelems[vw];

  return tot;
}

double ivfpq_unbalanced_factor(const ivfpq_t *ivf) {
  int vw;
  double tot = 0, uf = 0;

  for (vw = 0 ; vw < ivf->nbvw ; vw++) {
    tot += ivf->nbelems[vw];
    uf += ivf->nbelems[vw] * (double) ivf->nbelems[vw];
  }

  uf = uf * ivf->nbvw / (tot * tot);

  return uf;
}

/**********************************************************************************
 * centroids processing */



/* fill line i of n*d matrix out with line i from a minus line vw[i]
   from n*d matrix b. If single_a, take a single line from a */
static float *fmat_sub_lines (const float *a, int n, int d, const int *vw,
                              const float *b, int single_a, float *out)
{
  int i, j;
  float *ol = out;

  for (i = 0; i < n; i++) {
    const float *bl = b + d * vw[i];

    if(vw[i]>=0) /* <0 means NaNs */
      for (j = 0; j < d; j++)
        ol[j] = a[j] - bl[j];

    if (!single_a)
      a += d;
    ol += d;
  }
  return out;
}


static float * ivfpq_quantize (const ivfpq_t * ivf, int n, int d, const float * v, 
			    int * vw, int k, int nt)
{
  if (k == 1) {
    quantize_codebook_thread (n, ivf->nbvw, d, ivf->centroids, v, vw, nt, NULL, NULL);
    return NULL;
  }
  else {
    float * dis = quantize_codebook_multiple_thread
        (n, ivf->nbvw, d, k, ivf->centroids, v, vw, nt, NULL, NULL);
    return dis;
  }
}

ivfpq_t *ivfpq_new_learn (int nbvw, const float *centroids,
                          int d, int nsq, int nsq_bits,
                           const real * v, int n, int nb_iter_max,
                          int flags, int nt) {

  return ivfpq_new_learn_with_extract_map (nbvw, centroids,
                                           d, nsq, nsq_bits,
                                           v, NULL, n, nb_iter_max,
                                           flags, nt);

}

static double unbalanced_factor_vw(const int *vw,int n,int nbvw) {
  int *hist=ivec_new_0(nbvw); 
  int i;

  for(i=0;i<n;i++) if(vw[i]>=0) 
    hist[vw[i]]++;
  
  
  double tot=ivec_sum(hist,nbvw);

  return ivec_sum_2(hist,nbvw) * nbvw / (tot * tot);

}


ivfpq_t *ivfpq_new_learn_with_extract_map (int nbvw, const float *centroids,
                          int d, int nsq, int nsq_bits,
                           const real * v, const int* extract_map, int n, int nb_iter_max,
                          int flags, int nt)
{
  ivfpq_t *ivf = ivfpq_new (nbvw);
  ivf->flags = flags;

  ivf->centroids = (float *) malloc (nbvw * d * sizeof (*ivf->centroids));
  
  
  if (flags & IVFPQ_FLAG_CENTROID_PQ) {
    /* learn a specific pq_t for the assignment to the centroids, and encode/decode those */
    ivf->pq_assign = pq_learn (d, nsq, nsq_bits, centroids, nbvw, nb_iter_max);
    
    ivf->centroids_codes = pqcode_new (ivf->pq_assign, nbvw);
    pq_encode (ivf->pq_assign, centroids, nbvw, ivf->centroids_codes);
    
    pq_decode (ivf->pq_assign, ivf->centroids_codes, nbvw, ivf->centroids);
  }
  else memcpy (ivf->centroids, centroids, (nbvw * d * sizeof (*centroids)));


  if (flags & IVFPQ_FLAG_RELATIVE) {
    int *vw = ivec_new (n);
    //    quantize_codebook_thread (n, nbvw, d, centroids, v, vw, nt, NULL, NULL);
    //  -> replaced by:
    ivfpq_quantize (ivf, n, d, v, vw, 1, nt);
    
    printf("ivfpq learning set unbalanced factor %g\n",
          unbalanced_factor_vw(vw,n,nbvw));

    v = fmat_sub_lines (v, n, d, vw, centroids, 0, fvec_new (n * d));

    free (vw);
  }


  pq_t *pq = pq_learn_with_extract_map (d, nsq, nsq_bits, v, extract_map, n, nb_iter_max);

  if (flags & IVFPQ_FLAG_RELATIVE)
    free ((float *) v);

  ivf->pq = pq;

  return ivf;
}


/************************************************** building */



  


int * ivfpq_add_and_get_vw (ivfpq_t * ivf, const int *labels, const float *val, int n,
                int nt)
{
  const pq_t *pq = ivf->pq;
  int d = pq->d;

  int *vw = ivec_new (n);

  //  quantize_codebook_thread (n, ivf->nbvw, d, ivf->centroids, val, vw, nt, NULL, NULL);
  //  -> replaced by:
  ivfpq_quantize (ivf, n, d, val, vw, 1, nt);

  
  ivfpq_add_with_vw (ivf, labels, val, n, vw, nt);
 
  return vw;
} 

void ivfpq_add_with_vw (ivfpq_t * ivf, const int *labels, const float * val, int n, const int *vw, int nt) {

  const pq_t *pq = ivf->pq;
  int d = pq->d;

  if (ivf->flags & IVFPQ_FLAG_RELATIVE)
    val = fmat_sub_lines (val, n, d, vw, ivf->centroids, 0, fvec_new (n * d));

  /* TODO multithread this part also */
  ivfpq_add_vw (ivf, vw, labels, val, n);

  if (ivf->flags & IVFPQ_FLAG_RELATIVE)
    free ((float *) val);

}

void ivfpq_add (ivfpq_t * ivf, const int *labels, const float *val, int n,
                int nt)
{
  int *vw = ivfpq_add_and_get_vw(ivf,labels,val,n,nt);
  free(vw);
}


/************************************************** querying */



void ivfpq_query (const ivfpq_t * ivf, const float *val, int n,
                  int ma_q, double ma_disratio,
                  int k, int nt, int *labels, float *dists)
{

  const pq_t *pq = ivf->pq;
  int d = pq->d;

  int *vw = ivec_new (ma_q * n);


  float *  vw_centroid_dis = ivfpq_quantize (ivf, n, d, val, vw, ma_q, nt);

  float *v_shift =
      ivf->flags & IVFPQ_FLAG_RELATIVE ? fvec_new (ma_q * d) : NULL;

/*  double t1=getmillisecs(); */


  /* TODO multithread this part */
  int i, j;
  for (i = 0; i < n; i++) {
    int *vw_i = vw + i * ma_q;
    int nma = ma_q;

    if (ma_q > 1 && ma_disratio > 1)
      nma = compress_labels_by_disratio (vw_i, vw_centroid_dis + ma_q * i, 
					 nma, ma_disratio);

    const float *vq = val + i * d;
    if (ivf->flags & IVFPQ_FLAG_RELATIVE)
      vq = fmat_sub_lines (vq, nma, d, vw_i, ivf->centroids, 1, v_shift);


    int nres =
        ivfpq_query_vw (ivf, vw_i, nma, vq, k,
                        ivf->flags & IVFPQ_FLAG_RELATIVE,
                        labels + k * i, dists + k * i);

    /* clear out uninitialized results */
    for (j = nres; j < k; j++)
      labels[k * i + j] = -1;

  }

  free (v_shift);
  free (vw_centroid_dis);
  free (vw);
}


void ivfpq_query_and_get_vw (const ivfpq_t * ivf, const float *val, int n,
                             int ma_q, double ma_disratio,
                             int k, int nt, 
                             int *labels, float *dists, 
                             int *vw,int *label_map)
{

  const pq_t *pq = ivf->pq;
  int d = pq->d;

  float *  vw_centroid_dis = ivfpq_quantize (ivf, n, d, val, vw, ma_q, nt);

  float *v_shift =
      ivf->flags & IVFPQ_FLAG_RELATIVE ? fvec_new (ma_q * d) : NULL;

  /* TODO multithread this part */
  int i, j;
  for (i = 0; i < n; i++) {
    int *vw_i = vw + i * ma_q;
    int nma = ma_q;

    if (ma_q > 1 && ma_disratio > 1)
      nma = compress_labels_by_disratio (vw_i, vw_centroid_dis + ma_q * i, 
					 nma, ma_disratio);

    const float *vq = val + i * d;
    if (ivf->flags & IVFPQ_FLAG_RELATIVE)
      vq = fmat_sub_lines (vq, nma, d, vw_i, ivf->centroids, 1, v_shift);


    int nres =
        ivfpq_query_vw_get_label_map (ivf, vw_i, nma, vq, k,
                                      ivf->flags & IVFPQ_FLAG_RELATIVE,
                                      labels + k * i, dists + k * i,
                                      label_map + ma_q * i);

    /* clear out uninitialized results */
    for (j = nres; j < k; j++)
      labels[k * i + j] = -1;

    for (j = nma; j < ma_q; j++)
      label_map[ma_q*i+j]=-1;

  }

  free (v_shift);
  free (vw_centroid_dis);
}

#if 0

void ivfpqs_query(const ivfpq_t ** ivfs, int ndb,
                  const float *val, int n,
                  int ma,double ma_disratio, 
                  int k,  
                  int *labels, float *dists) {
  
  const ivfpq_t *ivf0=ivfs[0];
  int d=ivf0->pq->d;
  int alld=d*ndb;
  
  float *alldists=fvec_new(ndb*alld*ivf0->nbvw);

  /* compute distances to centroids */
  for(db=0;db<ndb;db++) {
    
    compute_cross_distances_nonpacked(d,
    


  }

  for(i=0;i<n;i++) {

    /* make new empty map */
      
    for(db=0;db<ndb;db++) {
      /* find ma closest centroids */

      /* loop over closest centroids */    
      for(c=0;c<ma;c++) {
        int w=;
        /* compute relative location */
        
        /* compute distances with all entries of inverted list */
        
        /* keep all (?) */
        for(i=0;i<ivf->codes[w];i++) {
          label=ivf->labels[w];
          dist=;

          ddist=dist-; /* delta with distance estimated from coarse quantizer */

          /* check hashtable for this label */

          if() { /* it exists */
            /* stored distance += ddist */
          } else {
            /* new entry for (label,ddist) */
          }
        }

      }

    }    

    /* make new maxheap (size k) */

    for(;;) { /* for all entries in map */

      /* sum up all distances from coarse quantizer */

      /* add accumulated ddist */

      /* update maxheap */

    }

    /* output maxheap result to labels and dists table */

  }

  free(alldists);

}

#endif

/************************************************** I/O */


#define WRITEANDCHECK(a,n) if(fwrite(a,sizeof(*a),n,f)!=n) {perror("ifvpq_fwrite"); abort(); }

void ivfpq_fwrite (FILE * f, const ivfpq_t * ivf)
{

  pq_fwrite (f, ivf->pq);

  WRITEANDCHECK (&ivf->nbvw, 1);
  WRITEANDCHECK (&ivf->flags, 1);
  WRITEANDCHECK (ivf->nbelems, ivf->nbvw);

  int i;
  for (i = 0; i < ivf->nbvw; i++) {
    WRITEANDCHECK (ivf->labels[i], ivf->nbelems[i]);
    //    WRITEANDCHECK ((char *) ivf->codes[i], ivf->nbelems[i] * pqcode_sizeof (ivf->pq));
    pqcode_fwrite (f, ivf->pq, ivf->codes[i], ivf->nbelems[i]);
  }

}

#undef WRITEANDCHECK


#define READANDCHECK(a,n) if(fread(a,sizeof(*a),n,f)!=n) {perror("ifvpq_read"); abort(); }

ivfpq_t *ivfpq_fread (FILE * f, const float *centroids)
{
  pq_t *pq = pq_fread (f);
  int nbvw;
  READANDCHECK (&nbvw, 1);

  ivfpq_t *ivf = ivfpq_new (nbvw);
  ivf->pq = pq;

  READANDCHECK (&ivf->flags, 1);

  READANDCHECK (ivf->nbelems, nbvw);

  int i;
  for (i = 0; i < ivf->nbvw; i++) {
    ivf->labels[i] = ivec_new (ivf->nbelems[i]);
    READANDCHECK (ivf->labels[i], ivf->nbelems[i]);

    ivf->codes[i] = pqcode_new_fread (f, ivf->pq, ivf->nbelems[i]);
    //    ivf->codes[i] = pqcode_new (pq, ivf->nbelems[i]);
    //    READANDCHECK ((char *) ivf->codes[i], ivf->nbelems[i] * pqcode_sizeof (pq));
  }
  
  int d=pq->d;
  
  ivf->centroids = (float *) malloc (nbvw * d * sizeof (*ivf->centroids));
  memcpy (ivf->centroids, centroids, (nbvw * d * sizeof (*centroids)));

  return ivf;
}

#undef READANDCHECK
