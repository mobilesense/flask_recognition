#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>


#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <unistd.h>

#include <siftgeo/siftgeo.h>
#include <yael/sorting.h>
#include <yael/machinedeps.h>
#include "ivfgeo.h"



#define IVFGEO_REALLOC_FACTOR 1.2
#define IVFGEO_REALLOC_NEWSIZE(oldsize) (((oldsize * 5) + 4) / 4)

#define NEWA(type,n) (type*)malloc(sizeof(type)*(n))
#define NEWAC(type,n) (type*)calloc(sizeof(type),(n))
#define NEW(type) NEWA(type,1)

#define MIN(a,b) ((a)<(b) ? (a) : (b))
#define MAX(a,b) ((a)>(b) ? (a) : (b))


/* number of distinct scales and angles resulting from the previous quantization scheme */
const int ivfgeo_nb_scales = (1 << NB_BITS_VAL_SCALE);
const int ivfgeo_nb_angles = (1 << NB_BITS_VAL_ANGLE);


static int uint8_nbones[256] = {
  0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
  1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
  1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
  2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
  1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
  2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
  2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
  3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 4, 5, 5, 6, 5, 6, 6, 7, 5, 6, 6, 7, 6, 7, 7, 8
};






/************************************************************************************
 * Matching Structure functions
 */



ivfgeo_match_t * ivfgeo_match_new (int maxn)
{
  ivfgeo_match_t * matches = (ivfgeo_match_t *) malloc (sizeof (ivfgeo_match_t));
  matches->n = 0;
  matches->nmax = maxn;
  matches->matches = (ivfgeo_match_elem_t *) malloc (matches->nmax * sizeof (ivfgeo_match_elem_t));
  
  return matches;
}


void ivfgeo_match_delete (ivfgeo_match_t * matches)
{
  free (matches->matches);
  free (matches);
}

inline void ivfgeo_match_add (ivfgeo_match_t * matches, int id, int dqid, 
			      int da, int ds, float score)
{
  int n = matches->n;

  assert(n==0 || matches->matches[n-1].dqid<=dqid);

  /* First check if the amount of memory is sufficient */
  if (n+1 >= matches->nmax) {
    int dqid=matches->matches[n-1].dqid;

    ivfgeo_match_elem_t * newmatch = (ivfgeo_match_elem_t *) realloc (matches->matches, (matches->nmax * 2)
								    * sizeof (ivfgeo_match_elem_t));
    if (!newmatch) {
      fprintf (stderr, "# Warning: unable to allocate memory for the matching structure\n");
      return;
    }
    assert(dqid==newmatch[n-1].dqid);

    /* A nasty bug somewhere triggers this assertion.
     *
     * The bug appears:
     *   - even without WGC/burstiness treatment 
     *   - even if this realloc is protected with a mutex
     *   - even if malloc and realloc sizes are rounded up to a factor of 32
     *   - with another version of the system (that from /home/clear/local, including Python 2.5.2)
     * It appears more probably:
     *   - if the initial nmax is bigger
     *   - with more threads
     *   - if the realloc increment is smaller (5/4 instead of 2)
     * It *seems to* disappear when:
     *   - realloc is replaced with malloc+memcpy+free
     *   - adding sizeof (ivfgeo_match_elem_t) to each malloc/realloc
     *   - doing the realloc before the array is full
     *
     * This last workaround is implemented here, hence the test n+1 >=
     * matches->nmax instead of n >= matches->nmax.
     */ 

    matches->nmax = matches->nmax * 2;
    matches->matches = newmatch;
  }

  /* Populate the structure the parameter match */
  
  matches->matches[n].id = id;
  matches->matches[n].dqid = dqid;
  matches->matches[n].ds = ds;
  matches->matches[n].da = da;
  matches->matches[n].score = score;

  assert(n==0 || matches->matches[n-1].dqid<=dqid);

  matches->n++;
}


void ivfgeo_match_display (const ivfgeo_match_t * matches)
{
/*
  int i, vw;

  printf ("[ ");
  for (i = 0 ; i < matches->nbvw ; i++) 
    if (matches->vw_start[i] != -1)
      printf ("%d %d %d / ", i, matches->vw_start[i], matches->vw_end[i]);
  printf (" ]");


  for (vw = 0 ; vw < matches->nbvw ; vw++)
    for (i = matches->vw_start[vw] ; i < matches->vw_end[vw] ; i++)

      fprintf (stdout, "i=%d vw=%-6d id=%-6d scale=%d angle=%d score=%.3f\n", i, vw, 
	       matches->matches[i].id, matches->matches[i].ds, 
	       matches->matches[i].da, matches->matches[i].score);
*/
}


/* A function to sort the matches according to their image id */
int ivfgeo_match_elem_compare_by_id_qid_score (const void * match1, const void * match2)
{
  int id_diff = ((ivfgeo_match_elem_t *) match1)->id - ((ivfgeo_match_elem_t *) match2)->id;

  if (id_diff != 0)
    return id_diff;

  int dqid_diff = ((ivfgeo_match_elem_t *) match1)->dqid - ((ivfgeo_match_elem_t *) match2)->dqid;

  if (dqid_diff !=0)
    return dqid_diff;

  return ((ivfgeo_match_elem_t *) match2)->score - ((ivfgeo_match_elem_t *) match1)->score;
}


/* A function to sort the matches according to their image id */
int ivfgeo_match_elem_compare_by_qid_id_score (const void * match1, const void * match2)
{
  int dqid_diff = ((ivfgeo_match_elem_t *) match1)->dqid - ((ivfgeo_match_elem_t *) match2)->dqid;

  if (dqid_diff !=0)
    return dqid_diff;

  int id_diff = ((ivfgeo_match_elem_t *) match1)->id - ((ivfgeo_match_elem_t *) match2)->id;

  if (id_diff != 0)
    return id_diff;

  return ((ivfgeo_match_elem_t *) match2)->score - ((ivfgeo_match_elem_t *) match1)->score;
}


/* reference implementation */
void ivfgeo_match_sort_by_id_qid_ref (ivfgeo_match_t * matches)
{
  qsort (matches->matches, matches->n, sizeof (ivfgeo_match_elem_t), 
	 ivfgeo_match_elem_compare_by_id_qid_score);
}


void ivfgeo_match_sort_by_qid_id (ivfgeo_match_t * matches)
{
  qsort (matches->matches, matches->n, sizeof (ivfgeo_match_elem_t), 
	 ivfgeo_match_elem_compare_by_qid_id_score);
}


void ivfgeo_match_sort_by_id_qid (ivfgeo_match_t * matches,int nbvec,int nqid) {

  long n=matches->n;
  ivfgeo_match_elem_t *tab=matches->matches;  
  ivfgeo_match_elem_t *tab2=NEWA(ivfgeo_match_elem_t,n);
  
  long i,q;
  long *im_i=NEWAC(long,nbvec);  
  long *qid_i=NEWA(long,nqid*2),*qid_end=qid_i+nqid;

  long nq=0;
  int prev_id=10000000;
  for(i=0;i<n;i++) {
    int id=tab[i].id;
    if(id<prev_id) 
      qid_i[nq++]=i;
    im_i[id]++;
    prev_id=id;
  }
  assert(nq<=nqid);
  
  for(q=0;q+1<nq;q++) 
    qid_end[q]=qid_i[q+1];
  qid_end[nq-1]=n;

  /* for q=0..nq-1, qid_i[q]..qid_end[q]-1 is a sequence of increasing
     id in tab (corresponding to a single query id) */

  long accu=0;
  for(i=0;i<nbvec;i++) {
    long prev=im_i[i];
    im_i[i]=accu;
    accu+=prev;
  }
  /* im_i[id] is the index in tab2 where data for image id begins */
  
  

  /* make blocks of ids such that the corresponding entries occupy
     less than blocksize. 
     We try to find a tradeoff between cache misses in tab and tab2.
  */
  const long blocksize=65536*1024/sizeof(ivfgeo_match_elem_t);

  int id_begin,id_end;
  for(id_begin=0;id_begin<nbvec;id_begin=id_end) {
    
    for(id_end=id_begin+1;id_end<nbvec;id_end++) 
      if(im_i[id_begin]+blocksize<im_i[id_end])
        break;

    /* fill-in the blocks with contributions from all segments. There
       should no more than nq cache misses */ 
    for(q=0;q<nq;q++) {
      for(i=qid_i[q];i<qid_end[q];i++) {
        int id=tab[i].id;
        if(id>=id_end) break;
        tab2[im_i[id]++]=tab[i];
      }
      qid_i[q]=i;
    }

  }

  free(im_i);
  free(qid_i);

  free(tab);

  matches->matches=tab2;
  matches->nmax=n;

}


/* suppress the matches that appears two times */
void ivfgeo_match_clean (ivfgeo_match_t * matches)
{
  long i, ii;

  /* first count the number of matches that have a score which is non-zero */
  for (i = 0, ii = 0 ; i < matches->n ; i++)
    if (matches->matches[i].score > 0) {
      matches->matches[ii] = matches->matches[i];
      ii++;
    }

  matches->n = ii;
}


/* Modify the weights of matches depending on their burstiness  */
void ivfgeo_match_weight_burstiness_intra (ivfgeo_match_t * matches,
                                           int burstiness,long nbvec,int nqid,
                                           int *map_dqid_to_vw)
{
  
 /*  ivfgeo_match_sort_by_id_qid_ref (matches); */

  ivfgeo_match_sort_by_id_qid (matches,nbvec,nqid); 
  
  long image_begin,image_end;

  int norm_type=(burstiness & BURSTNORM_INTRA_FUNC_MASK)>>8;
  
  /* loop over images */
  for(image_begin=0; image_begin<matches->n ; image_begin=image_end) {
    int id = matches->matches[image_begin].id;

    /* scan to find end of image */
    long i=image_begin+1;
    while(i<matches->n && matches->matches[i].id == id ) i++;  
    image_end=i;   

    long vw_begin,vw_end;

    /* loop over vw segments */
    for (vw_begin = image_begin ; vw_begin < image_end ; vw_begin = vw_end) {
      
      int vw = map_dqid_to_vw[matches->matches[vw_begin].dqid];
      
      /* find the end of this visual word segment */
      long i = vw_begin+1;      
      while (i<image_end && map_dqid_to_vw[matches->matches[i].dqid] == vw) i++;
      vw_end=i;
      
      if(vw_end==vw_begin+1) 
        continue; /* no burst */

      /* compute the total score and the max received in this segment */
      double max_score=0, tot_score = 0;
      for (i = vw_begin ; i < vw_end ; i++) {
        if(matches->matches[i].score>max_score) 
          max_score = matches->matches[i].score;
        tot_score += matches->matches[i].score;
      }


      /* update the scores */

      if(burstiness & BURSTNORM_CANCEL_MM) {  /* handle multiple matches */
        long dqid_begin,dqid_end;

        /* loop over scores with constant dqid */
        for(dqid_begin=vw_begin; dqid_begin<vw_end; dqid_begin=dqid_end) {
          int dqid = matches->matches[dqid_begin].dqid;
          i=dqid_begin+1;

          if(i<vw_end && matches->matches[i].dqid == dqid) { /* multiple matches */
            
            float best_score=matches->matches[dqid_begin].score;
            long best_score_index=dqid_begin;
                        
            /* find end & best score of multiple match */
            do {
              if(matches->matches[i].score>best_score) {
                best_score=matches->matches[i].score;
                best_score_index=i;
              }
              i++;
            } while(i<vw_end && matches->matches[i].dqid == dqid );
            dqid_end=i;   

            /* set all scores except best one to 0 */
            for(i=dqid_begin;i<dqid_end;i++) 
              if(i!=best_score_index) 
                matches->matches[i].score=0;

            query_stats.n_mm_cancel+=dqid_end-dqid_begin-1;

          } else {
            dqid_end=i;
          }         

        }

      }

      if (burstiness & BURSTNORM_INTRA) { /* downscore bursts */
        for (i = vw_begin ; i < vw_end ; i++) {
          double norm;
          switch(norm_type) {
          case 0: norm=sqrt (matches->matches[i].score / tot_score); break;
          case 1: norm=sqrt (matches->matches[i].score / max_score); break;
          case 2: norm=log (1+matches->matches[i].score / tot_score); break;
          case 3: norm=log (1+matches->matches[i].score / max_score); break;
          case 4: norm=sqrt (1.0 / (vw_end-vw_begin)); break;
          case 5: norm=1.0 / (vw_end-vw_begin); break;
          default: assert(0);
          }
          matches->matches[i].score *= norm;
        }
      }

    }
  }

  //  int old_n = matches->n;
  ivfgeo_match_clean (matches);
}



/* Modify the weights of matches depending on the amount of votes  */
void ivfgeo_match_weight_burstiness_inter (const ivfgeo_t * ivf, ivfgeo_match_t * matches, int npt,
                                           int burstiness)
{
  long i;
  int dqid;

  double * dq_tot_score = calloc (npt, sizeof (*dq_tot_score));
  
  int norm_type=(burstiness & BURSTNORM_INTRA_FUNC_MASK)>>16;

  int norm_base=norm_type & 0xf;
 
  /* Gather the score per query id */
  for (i = 0 ; i < matches->n ; i++) {
    dqid = matches->matches[i].dqid;
    if(norm_base==0) { /* sum */
      dq_tot_score [dqid] += matches->matches[i].score;
    } else if(norm_base==1) { /* max */
      if(matches->matches[i].score>dq_tot_score [dqid])
        dq_tot_score [dqid]=matches->matches[i].score;      
    } else if(norm_base==2) { /* count */
      dq_tot_score [dqid]+=1.0;
    }
  }

  int norm_func=(norm_type>>4)&0xf;
    
  /* Update the scores */
  for (i = 0 ; i < matches->n ; i++) {
    dqid = matches->matches[i].dqid;
    double norm;
    switch(norm_func) {
    case 0: norm=sqrt (matches->matches[i].score / dq_tot_score [dqid]); break;
    case 1: norm=1.0 / dq_tot_score [dqid]; break;
    case 2: norm=log(1+matches->matches[i].score / dq_tot_score [dqid]); break;
    default: assert(0);
    }
    matches->matches[i].score *= norm;
  }
  
  free (dq_tot_score);
}



/************************************************************************************
 * ivfgeo_t bookkeeping + geometry and binary signatures manipulation
 */

/* create a void inverted file in memory */
ivfgeo_t *ivfgeo_new (int nbvw, int maxn)
{
  int i;

  /* Default amount of data (nb of elements) allocated for a given vw */
  int default_seg_size = maxn * 100 / nbvw;

  ivfgeo_t *ivf = (ivfgeo_t *) malloc (sizeof (ivfgeo_t));
  assert (ivf);

  ivf->nbvw = nbvw;
  ivf->nbvec = 0;
  ivf->nbvec_max = maxn;
  ivf->norm_type = 2;

  ivf->segsize = (int *) malloc (sizeof (int) * ivf->nbvw);
  ivf->nbelems = (int *) malloc (sizeof (int) * ivf->nbvw);

  ivf->elems =
    (ivfgeo_elem_t **) malloc (sizeof (ivfgeo_elem_t *) * ivf->nbvw);
  assert (ivf->elems);

  /* a minimum segment size by default */
  if (default_seg_size < 10)
    default_seg_size = 10;

  for (i = 0; i < ivf->nbvw; i++) {
    ivf->segsize[i] = default_seg_size;
    ivf->nbelems[i] = 0;
    ivf->elems[i] =
      (ivfgeo_elem_t *) malloc (sizeof (ivfgeo_elem_t) * ivf->segsize[i]);
  }

  /* The structure that will contains the norms of the ivfgeo frequency vectors */
  ivf->norms = (distype_t *) malloc (ivf->nbvec_max * sizeof (*ivf->norms));
  assert (ivf->norms);

  /* The same for the labels */
  ivf->labels = (int *) malloc (ivf->nbvec_max * sizeof (*ivf->labels));
  assert (ivf->labels);

  /* No tf-idf scheme by default */
  ivf->tfidf = TFIDF_NO;  
  ivf->tfidf_weights =
    (distype_t *) malloc (ivf->nbvw * sizeof (*ivf->tfidf_weights));
  assert (ivf->tfidf_weights);
  ivfgeo_compute_tfidf(ivf);

  ivf->he_thres = NB_BITS_VAL_BINSIGN;
  ivf->he_dis_weights = NULL;

  ivf->wgc_type = AS_MAX_NEIGH_SIZE (3, 3);
  ivf->wgc_weights = NULL;

  ivf->scale_w = 0;
  ivf->scale_weights = NULL;

  ivf->he_dist_w = 0;
  ivf->burstiness = 0;

  ivf->mm_fd = -1;

  ivfgeo_compute_scale_weights (ivf);
  ivfgeo_compute_he_weights (ivf);

  return ivf;
}


void ivfgeo_delete (ivfgeo_t * ivf)
{
  int i;

  if (ivf->mm_fd > 0) {
    munmap (ivf->mm_base, ivf->mm_length);
    close (ivf->mm_fd);
  } else {
    for (i = 0; i < ivf->nbvw; i++)
      free (ivf->elems[i]);
  }
  free (ivf->nbelems);
  free (ivf->segsize);
  free (ivf->norms);
  free (ivf->labels);
  free (ivf->tfidf_weights);
  free (ivf->elems);
  free (ivf->wgc_weights);
  free (ivf->scale_weights);
  free (ivf->he_dis_weights);
  free (ivf);
}


/* The function used to quantize the scale and angle values */
inline int quantize_scale (geom_t * geom)
{
  distype_t scale = geom->scale;

  int ret = (int) (log (scale) / log (LOG_SCALE_QUANTIZATION_STEP) + 0.5);
  if (ret < 0)
    return 0;
  else if (ret > ivfgeo_nb_scales - 1)
    return ivfgeo_nb_scales - 1;

  return ret;
}



inline int quantize_angle (geom_t * geom)
{
  distype_t angle = geom->angle;

  int ret = (int) floor ((angle + M_PI) * ivfgeo_nb_angles / (2 * M_PI));
  /* assert(ret>=0 && ret<ivfgeo_nb_angles); */
  if (ret < 0) {
    return 0;
  } else if (ret >= ivfgeo_nb_angles) {
    return ivfgeo_nb_angles - 1;
  }

  return ret;
}


inline int binsign_hamming (binsign_t bs1, binsign_t bs2)
{
  int i, ham = 0;
  binsign_t diff = (bs1 ^ bs2) & ((1LL << NB_BITS_VAL_BINSIGN) - 1);

  for (i = 0; i < (NB_BITS_VAL_BINSIGN + 7) / 8; i++) {
    ham += uint8_nbones[diff & 255];
    diff >>= 8;
  }

  return ham;
}


void binsign_cross_distances(int sizeof_code,int na,int nb,
                             const char *a,
                             const char *b,
                             float *dists) {
  int i,j,k;

  for(i=0;i<nb;i++) {
    
    for(j=0;j<na;j++) {
      int accu=0;
      const char *bi=b+i*sizeof_code;
      const char *ai=a+j*sizeof_code;
      for(k=0;k<sizeof_code;k++) {
        accu+=uint8_nbones[(*ai ^ *bi) & 0xff];
        ai++; bi++;
      }
      *dists++=accu;
    } 

  }

}


/************************************************************************************
 * adding and merging of inverted files
 */


/*  add a vector an existing inverted file */
int ivfgeo_add (ivfgeo_t * ivf, const vwgeoset_t * vw, int label)
{
  int i, w, j;
  int id = ivf->nbvec;

  if (vw->n <= 0)
    return -1;

  if(ivf->nbvec>=(1<<NB_BITS_VAL_ID)) {
    fprintf(stderr,"ivfgeo_add: ivf too big, increase NB_BITS_VAL_ID\n");
    return -1;
  } 

  /* cannot increase size of mmaped inverted file */
  assert (ivf->mm_fd < 0);

  /* First check if the maximum number of vectors has been reached. 
     If so: geometric reallocation of the variable norms              */
  if (ivf->nbvec == ivf->nbvec_max) {
    ivf->nbvec_max = IVFGEO_REALLOC_NEWSIZE (ivf->nbvec);

    ivf->norms =
      (distype_t *) realloc (ivf->norms,
			     ivf->nbvec_max * sizeof (*ivf->norms));
    assert (ivf->norms);

    ivf->labels =
      (int *) realloc (ivf->labels, ivf->nbvec_max * sizeof (*ivf->labels));
    assert (ivf->labels);
  }

  /* compute the norm */

  if (ivf->norm_type == norm_none)
    ivf->norms[id] = 1.0;
  else if (ivf->norm_type == norm_type_1)
    ivf->norms[id] = vwgeoset_norm1 (vw);
  else if (ivf->norm_type == norm_type_2)
    ivf->norms[id] = vwgeoset_norm2 (vw);
  else
    assert (0);

  if (label == -1)
    ivf->labels[id] = id;
  else
    ivf->labels[id] = label;

  ivf->nbvec++;

  /* add descriptors from the vector */  
  for (i = 0; i < vw->n; i++) {
    w = vw->pts[i].vw;
    assert (w < ivf->nbvw);

    /* First check if a realloc is required or not */
    if (ivf->nbelems[w] + 1 == ivf->segsize[w]) {       /* -> YES, it is */

      ivf->segsize[w] = IVFGEO_REALLOC_NEWSIZE (ivf->segsize[w]);

      ivf->elems[w] = realloc (ivf->elems[w], (int) ivf->segsize[w] * sizeof (*ivf->elems[w]));
      assert (ivf->elems[w]);
    }

    j = ivf->nbelems[w];

    /* Store the id of the image and the *quantized* values of the scale and angle */
    ivf->elems[w][j].id = id;

    /* Quantization of the scales and angles */

    ivf->elems[w][j].scale = quantize_scale (&vw->pts[i].geom);
    ivf->elems[w][j].angle = quantize_angle (&vw->pts[i].geom);

    ivf->elems[w][j].binsign = vw->pts[i].binsign;

    /* update the number of elements */
    ivf->nbelems[w]++;
  }

  return id;
}


void ivfgeo_merge (ivfgeo_t * ivf, ivfgeo_t * ivf2, int label2ofs)
{
  int i, w;

  /* cannot increase size of mmaped inverted file */
  assert (ivf->mm_fd < 0);

  assert (ivf->nbvw == ivf2->nbvw);
  assert (ivf->norm_type == ivf2->norm_type);

  int n0 = ivf->nbvec;

  ivf->nbvec += ivf2->nbvec;

  if(ivf->nbvec>(1<<NB_BITS_VAL_ID)) {
    fprintf(stderr,"ivfgeo_merge: ivf too big, increase NB_BITS_VAL_ID\n");
    exit(1);
  } 

  /* First check if the maximum number of vectors has been reached. 
     If so: geometric reallocation of the variable norms              */
  if (ivf->nbvec >= ivf->nbvec_max) {
    ivf->nbvec_max = (int) (ivf->nbvec * IVFGEO_REALLOC_FACTOR);

    ivf->norms =
      (distype_t *) realloc (ivf->norms,
			     ivf->nbvec_max * sizeof (*ivf->norms));
    assert (ivf->norms);

    ivf->labels =
      (int *) realloc (ivf->labels, ivf->nbvec_max * sizeof (*ivf->labels));
    assert (ivf->labels);
  }

  for (i = 0; i < ivf2->nbvec; i++) {
    ivf->norms[i + n0] = ivf2->norms[i];
    ivf->labels[i + n0] = ivf2->labels[i] + label2ofs;
  }

  for (w = 0; w < ivf->nbvw; w++) {

    int to = ivf->nbelems[w];

    ivf->nbelems[w] += ivf2->nbelems[w];
    ivf->segsize[w] = ivf->nbelems[w];

    ivf->elems[w] =
      realloc (ivf->elems[w], ivf->segsize[w] * sizeof (*ivf->elems[w]));

    memcpy (ivf->elems[w] + to, ivf2->elems[w],
            sizeof (ivfgeo_elem_t) * ivf2->nbelems[w]);
    
    /* add an offset to vec indices in ivf2->elems */
    for (i = 0; i < ivf2->nbelems[w]; i++)
      ivf->elems[w][i+to].id += n0;

    if (ivf2->mm_fd < 0) {       /* free some precious memory */
      free (ivf2->elems[w]);
      ivf2->elems[w] = NULL;
    }
  }

  ivfgeo_delete (ivf2);
}

void ivfgeo_crop (ivfgeo_t * ivf, int maxlabel)
{
  int i, w;

  /* First check if the maximum number of vectors has been reached. 
     If so: geometric reallocation of the variable norms              */

  int maxid=-1;

  for (i = 0; i < ivf->nbvec; i++) if(ivf->labels[i]>=maxlabel) {
    maxid=i;
    break;
  }

  if(maxid<0) return; /* nothing to do */

  for(i=maxid;i<ivf->nbvec; i++)
    assert(ivf->labels[i]>=maxlabel);  

  for (w = 0; w < ivf->nbvw; w++) {

    ivfgeo_elem_t *elts=ivf->elems[w];

    int nkeep=0;
    for(i=0;i<ivf->nbelems[w];i++)
      if(elts[i].id<maxid) {
        elts[nkeep++]=elts[i];
      }
    ivf->nbelems[w]=nkeep;
  }

  ivf->nbvec=maxid;  
}



int ivfgeo_max_label (const ivfgeo_t * ivf)
{
  int i;
  int max_label = -1;

  for (i = 0; i < ivf->nbvec; i++)
    if (ivf->labels[i] > max_label)
      max_label = ivf->labels[i];

  return max_label;
}

/************************************************************************************
 * ivfgeo_t querying
 */




/* weighting of scales: bigger scales get higher weights */
void ivfgeo_compute_scale_weights_tab(int nb_scales,double scale_w,double **scale_weights_out) {

  int i;

  double  *scale_weights=*scale_weights_out;

  if (!scale_weights)
    scale_weights = NEWA (double, nb_scales);

  for (i = 0 ; i < nb_scales; i++) 
    scale_weights[i] = exp (i * scale_w);

  *scale_weights_out=scale_weights;
}

void ivfgeo_compute_scale_weights(ivfgeo_t * ivf) 
{
  ivfgeo_compute_scale_weights_tab(ivfgeo_nb_scales,ivf->scale_w,&ivf->scale_weights);  
}


/* weighting of hamming distances: smaller distances get higher weights */
void ivfgeo_compute_he_weights(ivfgeo_t * ivf) 
{
  int i;
  int n = NB_BITS_VAL_BINSIGN;

  if (!ivf->he_dis_weights)
    ivf->he_dis_weights = NEWA (double, n + 1);

  if (ivf->he_dist_w == 1) { /* weights depend on the distance */

    /* Theoretical values if the hypercube is uniformly filled */

    /* translated Matlab code
       a=zeros(1,65);
       for i=0:64 ;
       a(i+1)=nchoosek(64,i);
       end ;
       a=-log2(  );
    */
    unsigned long long cumsum=0L, nchoosek=1L;
      
    // "Standard" version using Shannon Information content
    for (i = 0 ; i < n + 1 ; i++) {
      cumsum += nchoosek;
      ivf->he_dis_weights[i] = n - log2 (cumsum);
      assert (finite (ivf->he_dis_weights[i]));
      nchoosek = nchoosek * (n - i) / (i + 1);
    }
  } else if (ivf->he_dist_w == 2) { /* Gaussian weighting */
    for (i = 0 ; i < n + 1 ; i++) {
      ivf->he_dis_weights[i] = 64.0 * exp (-i * (double) i / (double) (16 * 16));
      assert (finite (ivf->he_dis_weights[i]));
    }
  } else /* -> no weighting */
    for (i = 0; i < n + 1; i++)
      ivf->he_dis_weights[i] = 1.0;
}


query_stats_t query_stats={
  0,0,0,0,0,
  NULL,
  sizeof(ivfgeo_elem_t),sizeof(ivfgeo_match_elem_t)  
};



/* histograms are always 0D (= a single value) or 1D. 2D is not more
   efficient than 2 1D histograms */
#define HIST_DIM(as_max_type) ((as_max_type)==0 ? 0 : 1)

distype_t synth_wgc_tab (distype_t * tab, int synth_type,
				distype_t * weights);


/* querying using geometric consistancy */
int *ivfgeo_query (const ivfgeo_t * ivf,
                   const vwgeoset_t * vw,
                   distype_t ** retdis, int n)
{
  return ivfgeo_query_peek (ivf, vw, retdis, n, NULL, NULL);
}

/* estimation of query times for the various steps */

typedef struct {
  double b; /* begin */
  double a; /* total for this step */
} query_time_fraction_step_t;

static double query_step_end(const query_time_fraction_step_t*qs) {
  return qs->a+qs->b;
}

/* timings for the steps that will move the progress bar */
typedef struct {
  query_time_fraction_step_t invfile;
  query_time_fraction_step_t matches_vote;
  query_time_fraction_step_t synth_wgc;  
} query_time_fraction_t;

/* estimate timings */
static void estimate_query_time_fraction(const ivfgeo_t * ivf,int nqpt,
                                         query_time_fraction_t *qf) {
  int histd = HIST_DIM (ivf->wgc_type);
  
  /* reference timings for these sizes */
  int ref_nqpt=1402;
  int ref_nbvec=9000000;
  
  double nnratio=(nqpt*(double)ivf->nbvec)/(ref_nqpt*(double)ref_nbvec);
  double nbvecratio=ivf->nbvec/(double)ref_nbvec;

  /* ratios */
  qf->invfile.a=1025*nnratio;
  qf->matches_vote.a=ivf->burstiness ?  2200*nnratio : 0;
  qf->synth_wgc.a=histd ? 4100*nbvecratio : 0;

  /* normalize to 0.9 */
  double sum_a=qf->invfile.a+qf->matches_vote.a+qf->synth_wgc.a;
  double norm_a=0.9/sum_a;
  
  qf->invfile.a*=norm_a;
  qf->matches_vote.a*=norm_a;
  qf->synth_wgc.a*=norm_a;

  /* accumulate bs */
  qf->invfile.b=0.05;
  qf->matches_vote.b=query_step_end(&qf->invfile);
  qf->synth_wgc.b=query_step_end(&qf->matches_vote);

}


int *ivfgeo_query_peek (const ivfgeo_t * ivf,
                        const vwgeoset_t * vw,
                        distype_t ** retdis, int n,
                        void (*peek_fun) (void *arg, double frac),
                        void *peek_arg)
{
  int i, j, w, id;
  binsign_t binsign;
  float cornerness;
  int qscale, qangle;

  /* check if it is the null vector */
  if (vw->n == 0) {
    *retdis = NULL;
    return NULL;
  }

  query_time_fraction_t qtf;
  estimate_query_time_fraction(ivf,vw->n,&qtf);    

  int used_npt=vw->n;

  if(ivf->burstiness && used_npt >= (1<<NB_BITS_VAL_DQID)) {
    used_npt=(1<<NB_BITS_VAL_DQID);
    fprintf(stderr,"ivfgeo_query_peek: Too many (%d) query points, will use only %d. Increase NB_BITS_VAL_DQID\n",
            vw->n,used_npt);
  }

  /* the vector that will receive the distance values */
  const int na = ivfgeo_nb_angles, ns = 2 * ivfgeo_nb_scales - 1;
  int he_thres = ivf->he_thres;

  distype_t *dis = NEWAC (distype_t, ivf->nbvec);

  int histd = HIST_DIM (ivf->wgc_type);
  long ld_hist = histd == 0 ? 0 : na + ns;

  distype_t *dis_wgc = histd ? NEWAC (distype_t, ivf->nbvec * ld_hist) : NULL;

  assert (!histd || dis_wgc);
  assert (dis);

  /* ensure that the size of short-list to be returned is smaller 
     that the number of vectors in the ivfgeo structure */
  if (n > ivf->nbvec)
    n = ivf->nbvec;

  query_stats.n_images++;
  query_stats.n_points+=used_npt;

  /* create the structure that will receive the matched elements */
  ivfgeo_match_t * matches = ivf->burstiness ? ivfgeo_match_new (ivf->nbvec * 10) : NULL;

  long n_he_success=0;

  for (i = 0; i < used_npt ; i++) {

    /* Which visual word it is ? */
    w = vw->pts[i].vw;

    distype_t w2 = ivf->tfidf_weights[w] * ivf->tfidf_weights[w];

    binsign = vw->pts[i].binsign;
    cornerness = vw->pts[i].geom.cornerness;

    /* Quantization of the query scale and angle for this visual word */
    qscale = quantize_scale (&vw->pts[i].geom);
    qangle = quantize_angle (&vw->pts[i].geom);


    for (j = 0; j < ivf->nbelems[w]; j++) {

      int hd = binsign_hamming (binsign, ivf->elems[w][j].binsign);

      if (hd <= he_thres) {

	     id = ivf->elems[w][j].id;
       
        int ref_scale = MIN(qscale,ivf->elems[w][j].scale);

        double w3 = w2 * ivf->scale_weights[ref_scale] * ivf->he_dis_weights[hd];

        int da = (ivf->elems[w][j].angle - qangle) & (ivfgeo_nb_angles - 1);
        int ds = ivf->elems[w][j].scale - qscale + ivfgeo_nb_scales - 1;

        if(ivf->burstiness) {
          /* add this vote to the voting structure */
          ivfgeo_match_add (matches, id, i, da, ds, w3);
        } else if(histd!=0) {
          dis_wgc[ld_hist * id + da] += w3;
          dis_wgc[ld_hist * id + na + ds] += w3;          
        } else {
          dis[id] += w3;
        }

        n_he_success++;
      }
    }

    query_stats.n_visited_elements+=ivf->nbelems[w];

    if (peek_fun && i % 64 == 0)
      (*peek_fun) (peek_arg, i / (double) vw->n * qtf.invfile.a + qtf.invfile.b);
  }
  
  query_stats.n_he_success+=n_he_success;

/*
  for(i=0;i<matches->n;i++) 
    assert(matches->matches[i].dqid<used_npt);
*/

  if(ivf->burstiness) {

    query_stats.n_he_success+=matches->n;  

    /* Apply the match score update procedure */
    if (ivf->burstiness & BURSTNORM_INTRA) {
      int *map_dqid_to_vw=NEWA(int,vw->n);
      
      for(i=0;i<vw->n;i++) 
        map_dqid_to_vw[i]=vw->pts[i].vw;
      
      ivfgeo_match_weight_burstiness_intra (matches, ivf->burstiness, ivf->nbvec, vw->n, map_dqid_to_vw);
      
      free(map_dqid_to_vw);
    }
    
    if (ivf->burstiness & BURSTNORM_INTER)
      ivfgeo_match_weight_burstiness_inter (ivf, matches, vw->n, ivf->burstiness);

    /* Use the matching structure to vote */
    
    for (i = 0 ; i < matches->n ; i++) {
      
      unsigned int id = matches->matches[i].id;
      float sc = matches->matches[i].score;
      
      if (histd == 0) {
        dis[id] += sc;        
      } else { /* histd == 1 */
        
        unsigned char da = matches->matches[i].da;
        unsigned char ds = matches->matches[i].ds;
        
        dis_wgc[ld_hist * id + da] += sc;
        dis_wgc[ld_hist * id + na + ds] += sc;
      }

      if (peek_fun && i % 65536 == 0)
        (*peek_fun) (peek_arg, i / (double) matches->n * qtf.matches_vote.a + qtf.matches_vote.b);
    }
    
    ivfgeo_match_delete (matches);    

  }

  if(query_stats.as_hist && dis_wgc) {
    int tot=ld_hist*ivf->nbvec;
    for(i=0;i<tot;i++) 
      query_stats.as_hist[i]+=dis_wgc[i];
  }

  /* smooth the scale and angle assignation & compute scores */
  if (histd != 0) {
    for (i = 0; i < ivf->nbvec; i++) {
      dis[i] = synth_wgc_tab (dis_wgc + i * ld_hist, ivf->wgc_type, ivf->wgc_weights);
      
      if (peek_fun && i % 8192 == 0)
        (*peek_fun) (peek_arg, i / (double) ivf->nbvec * qtf.synth_wgc.a + qtf.synth_wgc.b);
      
    }
    free (dis_wgc);
  }

  /* compute the norm */
  if (ivf->norm_type != norm_none) {
    distype_t qnorm;

    if (ivf->norm_type == norm_type_1)
      qnorm = vwgeoset_norm1 (vw);
    else if (ivf->norm_type == norm_type_2)
      qnorm = vwgeoset_norm2 (vw);
    else if (ivf->norm_type == norm_type_1_sqrt) 
      qnorm = sqrt(vwgeoset_norm1 (vw));      
    else assert (0);

    assert(qnorm>0);

    for (i = 0; i < ivf->nbvec; i++)
      dis[i] = dis[i] / (ivf->norms[i] * qnorm);
  }

  int *ret=NULL;

  if(n<0) { /* return raw distance table */

    float *dislabel=NEWAC(float,ivfgeo_max_label(ivf)+1);
    
    /* remap labels */    
    for(i=0;i<ivf->nbvec;i++)
      dislabel[ivf->labels[i]]=dis[i];
    
    *retdis=dislabel;
    
  } else { /* return top n results */

    ret = NEWA(int,n);

    fvec_find_k_max (dis, ivf->nbvec, ret, n);
    *retdis = (distype_t *) malloc (n * sizeof (distype_t));
    
    for (i = 0; i < n; i++)
      (*retdis)[i] = dis[ret[i]];
    
    /* Now, replace the internal ids by the labels */
    for (i = 0; i < n; i++)
      ret[i] = ivf->labels[ret[i]];
    
  }

  free (dis);
    
  if (peek_fun)
    (*peek_fun) (peek_arg, 1);
  
  return ret;
}





/************************************************************************************
 * synthesis of geometric Hough tables 
 */

static double sqr (double x)
{
  return x * x;
}


void ivfgeo_set_wgc_type (ivfgeo_t * ivf, int wgc_type)
{
  ivf->wgc_type = wgc_type;
  ivfgeo_compute_wgc_tab(ivfgeo_nb_angles,ivfgeo_nb_scales,ivf->wgc_type,&ivf->wgc_weights);  
}

void ivfgeo_compute_wgc_tab (int nb_angles,int nb_scales,int wgc_type,
                             float **wgc_weights_out) 
{
  free(*wgc_weights_out);

  int premult_type=wgc_type & AS_MAX_PREMULTIPLY;

  if (premult_type) {

    /* Heuristic weighting */

    int na = nb_angles;
    int ns = 2 * nb_scales - 1;
    distype_t *bayes_coeffs = NEWA (distype_t, na * ns);

    int a, s;
    float ka1 = 4, ka2 = 0.4, ka3 = 0.2, k3 = 5, k4 = 1, k5 = 40;

    if(premult_type == AS_MAX_PREMULTIPLY_3)  
      ka2=0;

    if(premult_type == AS_MAX_PREMULTIPLY_1)
      ka2=ka3=0;

    for (a = 0; a < na; a++) {
      for (s = 0; s < ns; s++)
        bayes_coeffs[s + ns * a] =
	  1 +                                                          /* base coefficient */
	  ka1 * (exp (-sqr (a - 0) / k3) + exp (-sqr (a - na) / k3)) + /* peek at 0 & 2*pi */
          ka2 * exp (-sqr (a - na / 2) / k3) +                        /* peek at pi */
          ka3 * (exp (-sqr (a - na / 4) / k3) +                       /* peeks at +-pi/2 */
                 exp (-sqr (a - 3 * na / 4) / k3)) +
	  k4 * exp (-sqr (s - 31) / k5);                               /* peek at 0 for scale */

    }

    /* accumulate coefficients along the 2 dimensions */

    float *wgc_weights = NEWAC (distype_t, na + ns);
    for (a = 0; a < na; a++)
      for (s = 0; s < ns; s++) {
        distype_t v = bayes_coeffs[s + ns * a];
        wgc_weights[a] += v;
        wgc_weights[na + s] += v;
      }
    free (bayes_coeffs);
    *wgc_weights_out=wgc_weights;
  } 
  else
    *wgc_weights_out = NULL;
}


static int tab_is_0(distype_t * tab, int n) {
  int *ti=(int*)tab;
  int orv=0,i;
  n=n*sizeof(distype_t)/sizeof(int);

  for(i=0;i<n;i++) orv|=ti[i];
    
  return !orv;
}


static distype_t max_tab (distype_t * tab, int n)
{
  distype_t m = tab[0];
  int i;
  for (i = 1; i < n; i++)
    if (tab[i] > m)
      m = tab[i];
  return m;
}


/* m: neighbourhood size */
static distype_t *sum_blocks_1D (distype_t * tab, int n, int m)
{
  /* pass 1: make integral image */
  distype_t accu = 0;
  int i;
  for (i = 0; i < n; i++) {
    accu += tab[i];
    tab[i] = accu;
  }
  /* pass 2: replace last n-m values with sum over m */
  for (i = n - 1; i - m >= 0; i--) {
    tab[i] = tab[i] - tab[i - m];
  }
  /* values 0..m-2 are invalid */
  return tab + m - 1;
}




distype_t synth_wgc_tab (distype_t * dis_as,int wgc_type,distype_t * weights)
{

  distype_t dist;

  int na = ivfgeo_nb_angles;
  int ns = 2 * ivfgeo_nb_scales - 1;

  int histd = HIST_DIM (wgc_type);
  assert (histd == 1);
  int ld_hist = na + ns;

  /* quick check for imqges with no votes */
  if(tab_is_0(dis_as,ld_hist)) 
    return 0;

  if (wgc_type & AS_MAX_PREMULTIPLY) {
    int i;
    for (i = 0; i < ld_hist; i++)
      dis_as[i] *= weights[i];
  }

  int ad = (wgc_type >> 8) & 0xff;    /* angle smooth neighbourhood */
  int sd = (wgc_type >> 16) & 0xff;   /* scale smooth neighbourhood  */

  distype_t *sum_angles = dis_as;
  distype_t *sum_scales = dis_as + na;
  
  if (ad >= 2) {
    sum_angles = sum_blocks_1D (sum_angles, na, ad);
    na = na - ad + 1;
  }
  if (sd >= 2) {
    sum_scales = sum_blocks_1D (sum_scales, ns, sd);
    ns = ns - sd + 1;
  }
  
  distype_t max_a = max_tab (sum_angles, na);
  distype_t max_s = max_tab (sum_scales, ns);
  
  dist = MIN (max_a, max_s);

  return dist;
}



/************************************************************************************
 * Normalization coefficients
 */



/* Compute the norms associated with the elements of the inverted file */
void ivfgeo_compute_norms (ivfgeo_t * ivf)
{
  int i,w;

  if (ivf->norm_type == norm_none) {
    for (i = 0; i < ivf->nbvec; i++)
      ivf->norms[i] = 1.0;
    return;
  } 

  

  for (i = 0; i < ivf->nbvec; i++)
    ivf->norms[i] = 0;

#if 0

  /* reference implementation */

  distype_t *occ = (distype_t *) malloc (ivf->nbvec * sizeof (*ivf->norms));

  for (w = 0; w < ivf->nbvw; w++) {
    for (i = 0; i < ivf->nbvec; i++)
      occ[i] = 0;

    /* count the occurences of each vector for this visual word */
    for (i = 0; i < ivf->nbelems[w]; i++)
      occ[ivf->elems[w][i].id]++;

    if (ivf->norm_type == norm_type_1 || ivf->norm_type == norm_type_1_sqrt)
      for (i = 0; i < ivf->nbvec; i++)
        ivf->norms[i] += fabs (occ[i]); // * ivf->tfidf_weights[w];
    else if (ivf->norm_type == norm_type_2)
      for (i = 0; i < ivf->nbvec; i++)
        ivf->norms[i] += occ[i] * occ[i]; // * ivf->tfidf_weights[w] * ivf->tfidf_weights[w];
  }
  free (occ);

#else


  /* this implementation assumes that the elems table is ordered by id */
  for (w = 0; w < ivf->nbvw; w++) {

    /* count the occurences of each vector for this visual word */
    
    for (i = 0; i < ivf->nbelems[w]; ) {
      int id=ivf->elems[w][i].id;
      int i0=i++;
      
      while(i<ivf->nbelems[w] && ivf->elems[w][i].id==id) 
        i++;

      int occ=i-i0;

      if (ivf->norm_type == norm_type_1 || ivf->norm_type == norm_type_1_sqrt)
        ivf->norms[id] += fabs (occ);
      else if (ivf->norm_type == norm_type_2)        
        ivf->norms[id] += occ * (float)occ; // * ivf->tfidf_weights[w] * ivf->tfidf_weights[w];
      
    }

  }

#endif

  if (ivf->norm_type == norm_type_2 || ivf->norm_type == norm_type_1_sqrt)
    for (i = 0; i < ivf->nbvec; i++)
      ivf->norms[i] = sqrt (ivf->norms[i]);
}




void ivfgeo_compute_tfidf (ivfgeo_t * ivf)
{
  int w;
  if(ivf->tfidf==0) {
    for (w = 0; w < ivf->nbvw; w++)
      ivf->tfidf_weights[w] = 1.0;
  } else if(ivf->tfidf==1) {

/* Compute the tf-idf weights according to the criterion log(N/N_w),
   where N is the total number of descriptors and N_w is the number of 
   descriptors assigned to visual word w                               */

    long totelems = 0;
    
    for (w = 0; w < ivf->nbvw; w++)
      totelems += ivf->nbelems[w];
    
    for (w = 0; w < ivf->nbvw; w++)
      if (ivf->nbelems[w] == 0)
        ivf->tfidf_weights[w] = log (totelems);
      else
        ivf->tfidf_weights[w] = log (totelems / (distype_t) ivf->nbelems[w]);
  } else
    assert(0);
}


/************************************************************************************
 * ivfgeo_t I/O
 */

/* Macros to handle the i/O of the ivfgeo structure */
#define IVFGEO_READ_ERROR(ret, expected_ret)				\
  {									\
    if (ret != (expected_ret)) {					\
      fprintf (stderr, "# Unable to read the elements from the inverted file %s\n", \
	       filename);						\
      free (ivf);							\
      return NULL;							\
    }									\
  }

#define IVFGEO_WRITE_ERROR(ret, expected_ret)				\
  {									\
    if (ret != (expected_ret)) {					\
      fprintf (stderr, "# Unable to write the elements of the inverted file %s\n", \
	       filename);						\
      return;								\
    }									\
  }



/* open or close an existing inverted file for reading or insertion */
ivfgeo_t *ivfgeo_read (const char *filename, int fmt)
{
  int i, j, ret;
  ivfgeo_t *ivf = (ivfgeo_t *) malloc (sizeof (ivfgeo_t));
  assert (ivf);

  FILE *f = fopen (filename, "r");
  if (!f) {
    fprintf (stderr, "# Unable to open inverted file %s for reading\n",
             filename);
    return NULL;
  }

  ret = fread (&(ivf->nbvw), sizeof (ivf->nbvw), 1, f);
  IVFGEO_READ_ERROR (ret, 1);

  ret = fread (&(ivf->nbvec), sizeof (ivf->nbvec), 1, f);
  IVFGEO_READ_ERROR (ret, 1);

  if(ivf->nbvec>(1<<NB_BITS_VAL_ID)) {
    fprintf(stderr,"siftgeo_read: nbvec too big for NB_BITS_VAL_ID\n");
    assert(0);
  }    

  ret = fread (&(ivf->norm_type), sizeof (ivf->norm_type), 1, f);
  IVFGEO_READ_ERROR (ret, 1);

  ivf->nbvec_max = IVFGEO_REALLOC_NEWSIZE (ivf->nbvec);

  ivf->segsize = (int *) malloc (sizeof (*ivf->segsize) * (ivf->nbvw));
  ivf->nbelems = (int *) malloc (sizeof (*ivf->segsize) * (ivf->nbvw));

  ivf->elems =
    (ivfgeo_elem_t **) malloc (sizeof (ivfgeo_elem_t *) * ivf->nbvw);
  assert (ivf->elems);

  /* Read the number of elements stored for each visual word */
  ret = fread (ivf->nbelems, sizeof (*ivf->nbelems), ivf->nbvw, f);
  IVFGEO_READ_ERROR (ret, ivf->nbvw);

  /* Read the norms of the frequency vectors (=number of descriptors).
     Note that no extra memory is used here */
  ivf->norms = (distype_t *) malloc (ivf->nbvec_max * sizeof (*ivf->norms));
  assert (ivf->norms);

  ret = fread (ivf->norms, sizeof (*ivf->norms), ivf->nbvec, f);
  IVFGEO_READ_ERROR (ret, ivf->nbvec);


  /* Read the labels associated with the frequency vectors */
  ivf->labels = (int *) malloc (ivf->nbvec_max * sizeof (*ivf->labels));
  assert (ivf->labels);

  ret = fread (ivf->labels, sizeof (*ivf->labels), ivf->nbvec, f);
  IVFGEO_READ_ERROR (ret, ivf->nbvec);

  char *mm = NULL;

  if (fmt == 1 || fmt==3) {
    size_t totlen = 0;
    int i;
    for (i = 0; i < ivf->nbvw; i++)
      totlen += ivf->nbelems[i];

    int fd = open (filename, O_RDONLY);
    if (fd < 0) {
      perror ("open before mmap failed");
      exit (1);
    }
    /* map whole file */
    struct stat sb;
    fstat (fd, &sb);
    mm = mmap (NULL, sb.st_size, PROT_READ, MAP_SHARED, fd, 0);
    if (mm == MAP_FAILED) {
      perror ("map failed");
      exit (1);
    }
    ivf->mm_fd = fd;
    ivf->mm_base = mm;
    ivf->mm_length = sb.st_size;
    mm += ftell (f);
  } else {
    ivf->mm_fd = -1;
  }

  int fmt_3_sum=0;

  /* Read the data corresponding to visual words */
  for (i = 0; i < ivf->nbvw; i++) {
    int tmpint;
    unsigned char tmpuchar;
    binsign_t binsign;

    /* When reading from a file, by default behavior is to adjust the size 
       of the segment to the exact number of elements to be loaded          */

    /* Now, read the elements in a raw */
    if (fmt == 0) {
      ivf->segsize[i] = (int) (ivf->nbelems[i] * IVFGEO_REALLOC_FACTOR) + 1;
      ivf->elems[i] =
	(ivfgeo_elem_t *) malloc (sizeof (ivfgeo_elem_t) * ivf->segsize[i]);
      ret = 0;
      for (j = 0; j < ivf->nbelems[i]; j++) {
        ret += fread (&tmpint, sizeof (tmpint), 1, f);
        ivf->elems[i][j].id = tmpint;
        ret += fread (&tmpuchar, sizeof (tmpuchar), 1, f);
        ivf->elems[i][j].scale = tmpuchar;
        ret += fread (&tmpuchar, sizeof (tmpuchar), 1, f);
        ivf->elems[i][j].angle = tmpuchar;
        ret += fread (&binsign, sizeof (binsign), 1, f);
        ivf->elems[i][j].binsign = binsign;
      }
      IVFGEO_READ_ERROR (ret, ivf->nbelems[i] * 4);
    } else if (fmt == 1 || fmt == 3) {
      ivf->segsize[i] = ivf->nbelems[i];
      ivf->elems[i] = (void *) mm;
      mm += ivf->nbelems[i] * sizeof (ivfgeo_elem_t);
      assert ((void*)mm <= ivf->mm_base + ivf->mm_length);
      if (fmt==3) {
        for(j=0;j<ivf->nbelems[i];j++) 
          fmt_3_sum+=ivf->elems[i][j].id;
      }
    } else if (fmt == 2) {
      ivf->segsize[i] = (int) (ivf->nbelems[i] * IVFGEO_REALLOC_FACTOR) + 1;
      ivf->elems[i] =
	(ivfgeo_elem_t *) malloc (sizeof (ivfgeo_elem_t) * ivf->segsize[i]);
      ret = fread (ivf->elems[i], sizeof (ivfgeo_elem_t), ivf->nbelems[i], f);
      IVFGEO_READ_ERROR (ret, ivf->nbelems[i]);
    }
  }

  if(fmt==3) {
    /* make sure the compiler does not optimize out the computed value... */
    printf("format 3 sum=%d\n",fmt_3_sum);
  }

  /* Initialize the variable associated with the tf-idf weighting */
  ivf->tfidf = TFIDF_NO;
  ivf->tfidf_weights =
    (distype_t *) malloc (ivf->nbvw * sizeof (*ivf->tfidf_weights));
  ivfgeo_compute_tfidf(ivf);


  fclose (f);

  ivf->he_dist_w = 0;
  ivf->burstiness = 0;

  ivf->he_thres = NB_BITS_VAL_BINSIGN;
  ivf->he_dis_weights = NULL;

  ivf->wgc_type = AS_MAX_NEIGH_SIZE (3, 3);
  ivf->wgc_weights = NULL;

  ivf->scale_w = 0;
  ivf->scale_weights = NULL;

  ivfgeo_compute_scale_weights (ivf);
  ivfgeo_compute_he_weights (ivf);

  return ivf;
}


/* write the inverted file in a file */
void ivfgeo_write (const char *filename, const ivfgeo_t * ivf, int fmt)
{
  int i, j, ret;

  FILE *f = fopen (filename, "w");
  if (!f) {
    fprintf (stderr, "# Unable to create the new inverted file %s\n",
             filename);
    return;
  }

  ret = fwrite (&(ivf->nbvw), sizeof (ivf->nbvw), 1, f);
  IVFGEO_WRITE_ERROR (ret, 1);

  printf ("nbvw = %d\n", ivf->nbvw);

  ret = fwrite (&(ivf->nbvec), sizeof (ivf->nbvec), 1, f);
  IVFGEO_WRITE_ERROR (ret, 1);

  ret = fwrite (&(ivf->norm_type), sizeof (ivf->norm_type), 1, f);
  IVFGEO_WRITE_ERROR (ret, 1);

  /* Write the number of elements stored for each visual word */
  ret = fwrite (ivf->nbelems, sizeof (*ivf->nbelems), ivf->nbvw, f);
  IVFGEO_WRITE_ERROR (ret, ivf->nbvw);

  /* Write the norms and the labels of the frequency vectors */
  ret = fwrite (ivf->norms, sizeof (*ivf->norms), ivf->nbvec, f);
  IVFGEO_WRITE_ERROR (ret, ivf->nbvec);

  ret = fwrite (ivf->labels, sizeof (*ivf->labels), ivf->nbvec, f);
  IVFGEO_WRITE_ERROR (ret, ivf->nbvec);

  for (i = 0; i < ivf->nbvw; i++) {
    binsign_t binsign;
    int tmpint;
    unsigned char tmpuchar;


    if (fmt == 0) {
      ret = 0;
      for (j = 0; j < ivf->nbelems[i]; j++) {
        tmpint = ivf->elems[i][j].id;
        ret += fwrite (&tmpint, sizeof (tmpint), 1, f);

        tmpuchar = ivf->elems[i][j].scale;
        ret += fwrite (&tmpuchar, sizeof (tmpuchar), 1, f);

        tmpuchar = ivf->elems[i][j].angle;
        ret += fwrite (&tmpuchar, sizeof (tmpuchar), 1, f);

        binsign = ivf->elems[i][j].binsign;
        ret += fwrite (&binsign, sizeof (binsign), 1, f);
      }
      IVFGEO_WRITE_ERROR (ret, ivf->nbelems[i] * 4);
    } else {
      ret =
	fwrite (ivf->elems[i], sizeof (ivfgeo_elem_t), ivf->nbelems[i], f);
      IVFGEO_WRITE_ERROR (ret, ivf->nbelems[i]);
    }
  }

  fclose (f);
}



/* display the contents of an inverted file  */
void ivfgeo_display (const ivfgeo_t * ivf)
{
  int i, j, k;
  binsign_t id, scale, angle, binsign;
  printf ("Nb visual words   %d\n", ivf->nbvw);
  printf ("Nb vectors        %d\n", ivf->nbvec);


  printf ("Norms ");
  for (i = 0; i < ivf->nbvec; i++)
    printf ("%.3f ", ivf->norms[i]);
  printf ("\n\n");

  printf ("Labels ");
  for (i = 0; i < ivf->nbvec; i++)
    printf ("%d ", ivf->labels[i]);
  printf ("\n\n");


  /* for each segment, display the contents */
  for (i = 0; i < ivf->nbvw; i++) {

    if (ivf->nbelems[i] > 0)
      fprintf (stdout, "[ Visual word %d ] %d elements (segsize: %d)\n", i,
               ivf->nbelems[i], ivf->segsize[i]);
    else {
      continue;
    }

    for (j = 0; j < ivf->nbelems[i]; j++) {
      id = ivf->elems[i][j].id;
      scale = ivf->elems[i][j].scale;
      angle = ivf->elems[i][j].angle;
      binsign = ivf->elems[i][j].binsign;

      printf ("   %d %d %d ", (int) id, (int) scale, (int) angle);


      for (k = 0; k < NB_BITS_VAL_BINSIGN; k++)
        printf ("%d", (int) ((binsign >> k) & 1));
      printf ("\n");

      /* printf ("%llx\n", binsign); */

    }
    printf ("\n");
  }
}

void ivfgeo_summary (const ivfgeo_t * ivf)
{
  int i;

  int nb   = 0;
  int size = 0;

  printf ("Nb visual words   %d\n", ivf->nbvw);
  printf ("Nb vectors        %d\n", ivf->nbvec);
  printf ("Max label         %d\n", ivfgeo_max_label(ivf));

  for (i = 0; i < ivf->nbvw; i++) {
    if (ivf->nbelems[i] > 0) {
      fprintf (stdout, "[ Visual word %d ] %d elements (segsize: %d)\n", i, ivf->nbelems[i], ivf->segsize[i]);
      nb   += ivf->nbelems[i];
      size += ivf->segsize[i];
    }
  }

  printf("Total elements     %d (size %d)\n", nb, size);
}

void ivfgeo_check(ivfgeo_t * ivf) {
  int i,j;
  assert(ivf->nbvw>=0);
  assert(ivf->nbvec>=0);
  assert(ivf->nbvec_max>=ivf->nbvec);
  assert(ivf->norm_type>=0 && ivf->norm_type<=3);
  for(i=0;i<ivf->nbvw;i++) {
    assert(ivf->nbelems[i]>=0 && ivf->nbelems[i]<=ivf->segsize[i]);
    for(j=0;j<ivf->nbelems[i];j++) 
      assert(ivf->elems[i][j].id>=0 &&
             ivf->elems[i][j].id<ivf->nbvec);
  }  
}


void ivfgeo_filter_duplicate_vw (ivfgeo_t * ivf)
{
  int i, j, id;

  for (i = 0; i < ivf->nbvw; i++) {

    int n = ivf->nbelems[i];
    int jj = 1;

    if (n < 2)
      continue;

    for (j = 1; j < n ; j++) {
      id = ivf->elems[i][j].id;

      /* same image as before? */
      if (id == ivf->elems[i][j - 1].id) 
	continue;

      ivf->elems[i][jj].id = ivf->elems[i][j].id;
      jj++;
    }
    ivf->nbelems[i] = jj;
  }
}

void ivfgeo_mask_binsigns (ivfgeo_t * ivf,unsigned long long mask) {
  int i, j;
  
  for (i = 0; i < ivf->nbvw; i++) {

    int n = ivf->nbelems[i];

    for (j = 1; j < n ; j++) {
      ivf->elems[i][j].binsign &= mask;
    }
  } 

}



/*----------------- Statistical Information ----------------*/

/* Count the total number of elements (descriptors) in the inverted file */
int ivfgeo_count_nbelems (const ivfgeo_t * ivf)
{
  int tot = 0, vw;
  
  for (vw = 0 ; vw < ivf->nbvw ; vw++)
    tot += ivf->nbelems[vw];

  return tot;
}


/* Compute the "unbalanced factor" of the inverted file */
double ivfgeo_unbalanced_factor (const ivfgeo_t * ivf)
{
  int vw;
  double tot = 0, uf = 0;

  for (vw = 0 ; vw < ivf->nbvw ; vw++) {
    tot += ivf->nbelems[vw];
    uf += ivf->nbelems[vw] * (double) ivf->nbelems[vw];
  }

  uf = uf * ivf->nbvw / (tot * tot);

  return uf;
}



/* compute the natural burstiness of visual words (use, e.g., maxcount=1000)  */
int * ivfgeo_compute_burstiness_hist (const ivfgeo_t * ivf, int maxcount)
{
  int i, j, id, totseg = 0;

  /* histogram of the number of visual words that appears */
  int * count = (int *) calloc (maxcount + 1, sizeof (int));

  for (i = 0; i < ivf->nbvw; i++) {

    int lastid = -1, curcount = 1;
    int n = ivf->nbelems[i];

    for (j = 0; j < n ; j++) {
      id = ivf->elems[i][j].id;

      /* same image as before? */
      if (id == lastid) 
	curcount++;

      /* otherwise, populate the histogram */
      else {
	if (curcount < maxcount)
	  count[curcount]++;
	else
	  count[maxcount]++;

	curcount = 1;
	lastid = id;
	totseg++;
      }
    }

    if (curcount < maxcount)
      count[curcount]++;
    else
      count[maxcount]++;
    totseg++;
  }

  return count;
}


/* Return an histogram of the difference between images id */
int * ivfgeo_histo_diff_im_id (const ivfgeo_t * ivf)
{
  int i, j, d;
  int * histo = (int *) calloc (ivf->nbvec, sizeof (*histo));

  for (i = 0; i < ivf->nbvw; i++) {
    int pos = -1;
    for (j = 0 ; j < ivf->nbelems[i] ; j++) {
      d = ivf->elems[i][j].id - pos;
      histo[d]++;
      pos = ivf->elems[i][j].id;
    }
  }
  return histo;
}
