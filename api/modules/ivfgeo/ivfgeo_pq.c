#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <stdio.h>

#include <yael/vector.h>
#include <yael/sorting.h>

#include <siftgeo/siftgeo_and_nn.h>


 
#include "ivfgeo_pq.h"



#define NEWA(type,n) (type*)malloc(sizeof(type)*(n))
#define NEWAC(type,n) (type*)calloc(n,sizeof(type))
#define NEW(type) NEWA(type,1)

#define MIN(a,b) ((a)<(b) ? (a) : (b)) 

ivfgeo_pq_t *ivfgeo_pq_new(ivfpq_t *ivfpq) {
  ivfgeo_pq_t *ivf=NEWAC(ivfgeo_pq_t,1);

  ivf->ivfpq=ivfpq; 

  /* clear out contents of ivfpq */ 
  int i;
  for(i=0;i<ivfpq->nbvw;i++) 
    ivfpq->nbelems[i]=0;

  ivf->norm_type=norm_type_2;
  ivf->tfidf=TFIDF_NO;  

  /* decent defaults for query params */
  
  ivf->ma=1;
  ivf->ma_disratio=1.2;
  ivf->k=10;
  ivf->nthread=1; 
  
  ivfgeo_pq_compute_tfidf(ivf);
  ivfgeo_pq_compute_scale_weights(ivf);  
  ivfgeo_pq_set_wgc_type(ivf,0);

  return ivf;
}

void ivfgeo_pq_display(const ivfgeo_pq_t * ivf) {
  int i;
  ivfpq_display(ivf->ivfpq); 
  
  printf("nim=%d tot_pts=%ld tfidf=%d\n",ivf->nim,ivf->tot_pts,ivf->tfidf); 

  printf("norms=[ "); 
  for(i=0;i<ivf->nim;i++) printf("%g ",ivf->norms[i]);
  printf("]\nim_ends=[ "); 
  for(i=0;i<ivf->nim;i++) printf("%d ",ivf->im_ends[i]);
  printf("]\npoints=[ "); 
  for(i=0;i<ivf->tot_pts;i++) printf("%d,%d ",ivf->points[i].angle,ivf->points[i].scale);
  printf("]\n");

}

/******************************************************* adding */

/* computes norm & changes set! */
static double vw_tab_norm1(int *vw, long n) {
  if(n==0) return 1.0;
  return n;
}

static double vw_tab_norm2(int *vw, long n) {
  if(n==0) return 1.0;
  int i;
  ivec_sort(vw,n);

  /* step over -1's from the multiple assignement */
  for(i=0;i<n;i++) if(vw[i]>=0) break;
  
  int prev_vw=-1,nocc=0;
  double accu=0;
  for(;i<n;i++) {
    if(vw[i]==prev_vw) 
      nocc++;
    else {
      accu+=nocc*nocc;
      prev_vw=vw[i];
      nocc=1;
    }
  }
  accu+=nocc*nocc;

  return sqrt(accu); 
}


void ivfgeo_pq_compute_tfidf (ivfgeo_pq_t * ivf)
{
  int nbvw=ivf->ivfpq->nbvw;

  if(!ivf->tfidf_weights) 
    ivf->tfidf_weights=NEWA(float,nbvw);

  int w;
  if(ivf->tfidf==0) {
    for (w = 0; w < nbvw; w++)
      ivf->tfidf_weights[w] = 1.0;
  } else if(ivf->tfidf==1) {
    
    const int *nbelems=ivf->ivfpq->nbelems;

/* Compute the tf-idf weights according to the criterion log(N/N_w),
   where N is the total number of descriptors and N_w is the number of 
   descriptors assigned to visual word w                               */

    long totelems = 0;
    
    for (w = 0; w < nbvw; w++)
      totelems += nbelems[w];
    
    for (w = 0; w < nbvw; w++)
      if (nbelems[w] == 0)
        ivf->tfidf_weights[w] = log (totelems);
      else
        ivf->tfidf_weights[w] = log (totelems / (distype_t) nbelems[w]);
  } else
    assert(0);
}




void ivfgeo_pq_compute_norms (ivfgeo_pq_t * ivf) {

  assert (ivf->norm_type == norm_none || !"not implemented");

  int i;

  for(i=0;i<ivf->nim;i++) 
    ivf->norms[i]=1.0;
}


void ivfgeo_pq_compute_scale_weights(ivfgeo_pq_t * ivf)
{
  ivfgeo_compute_scale_weights_tab(ivfgeo_nb_scales,ivf->scale_w,&ivf->scale_weights); 
}

void ivfgeo_pq_set_wgc_type (ivfgeo_pq_t * ivf, int wgc_type)
{
  ivf->wgc_type = wgc_type;
  ivfgeo_compute_wgc_tab(ivfgeo_nb_angles,ivfgeo_nb_scales,ivf->wgc_type,&ivf->wgc_weights);  
}



int ivfgeo_pq_add (ivfgeo_pq_t * ivf, const pointset_t * ps) { 
  return ivfgeo_pq_add_with_vw (ivf, ps, NULL);
}


static void ivfgeo_pq_add_pts(ivfgeo_pq_t * ivf, int add_im,int add_pts) {
  
  ivf->nim+=add_im;
  ivf->tot_pts+=add_pts;

  if(ivf->nim>=ivf->na_img) {
    ivf->na_img=ivf->nim<8 ? 8 : ivf->nim*3/2;
    ivf->norms=realloc(ivf->norms,sizeof(*ivf->norms)*ivf->na_img);
    ivf->im_ends=realloc(ivf->im_ends,sizeof(*ivf->im_ends)*ivf->na_img);
  }


  if(ivf->tot_pts>=ivf->na_pts) {
    ivf->na_pts=ivf->tot_pts<8 ? 8 : ivf->tot_pts*3/2;
    ivf->points=realloc(ivf->points,sizeof(*ivf->points)*ivf->na_pts);
  }

}

int ivfgeo_pq_add_with_vw (ivfgeo_pq_t * ivf, const pointset_t * ps, const int *vw_in) {
  int i;
  int d=ivf->ivfpq->pq->d;

  int id=ivf->nim;
  int pt0=ivf->tot_pts;  

  ivfgeo_pq_add_pts(ivf,1,ps->n);

  int *vw;

  if(ps->n==0) {
    vw=NULL;
  } else { 

    float *descs=fvec_new(d*ps->n);
    int d2=pointset_into_fvecs(ps,descs);
    
    assert(d==d2);
    
    /* int *labels=ivec_new_set(ps->n, id); */
    int *labels=ivec_new_range(pt0,ivf->tot_pts+ps->n);
    
    int nthread=ivf->nthread;

    if(!vw_in) 
      vw=ivfpq_add_and_get_vw (ivf->ivfpq, labels, descs, ps->n, nthread);
    else {
      vw=(int*)vw_in;
      ivfpq_add_with_vw (ivf->ivfpq, labels, descs, ps->n, vw, nthread);
    } 
    free(labels);
    free(descs);
  } 


  /* handle norms */

  if (ivf->norm_type == norm_none)
    ivf->norms[id] = 1.0;
  else if (ivf->norm_type == norm_type_1)
    ivf->norms[id] = vw_tab_norm1 (vw, ps->n);
  else if (ivf->norm_type == norm_type_2)
    ivf->norms[id] = vw_tab_norm2 (vw, ps->n);
  else
    assert (0);
 
  /* fill in quantized points */

  ivf->im_ends[id]=ivf->tot_pts;
  
  for(i=0;i<ps->n;i++) {
    ivfgeo_pq_point_t*pt=ivf->points+pt0+i;
    pt->scale=quantize_scale(&ps->pts[i].geom);
    pt->angle=quantize_angle(&ps->pts[i].geom);
  }  

  if(!vw_in) 
    free(vw);

  return id;
}

ivfgeo_pq_t *ivfgeo_pq_dup(ivfgeo_pq_t * ivf) {
  assert(ivf->nim==0);

  ivfgeo_pq_t *ivf2=NEW(ivfgeo_pq_t);

  /* works because im and pts arrays are NULL */
  *ivf2=*ivf; 

  ivf2->ivfpq=ivfpq_dup(ivf->ivfpq);

  return ivf2;
} 


void ivfgeo_pq_merge(ivfgeo_pq_t * ivf,ivfgeo_pq_t *ivf2) {
  int i,j;

  int nim=ivf->nim;
  int tot_pts=ivf->tot_pts;

  /* shift labels */  
  ivfpq_t *ivfpq=ivf->ivfpq;
  ivfpq_t *ivfpq2=ivf2->ivfpq;
  
  for(i=0;i<ivfpq->nbvw;i++) {
    for(j=0;j<ivfpq2->nbelems[i];j++) 
      ivfpq2->labels[i][j]+=tot_pts;
  }

  /* shift ends */
  for(i=0;i<ivf2->nim;i++) 
    ivf2->im_ends[i]+=tot_pts;  

  ivfpq_merge(ivf->ivfpq,ivf2->ivfpq);

  ivfgeo_pq_add_pts(ivf,ivf2->nim,ivf2->tot_pts);

  /* copy norms */
  memcpy(ivf->norms+nim,ivf2->norms,ivf2->nim*sizeof(*ivf->norms));
  memcpy(ivf->im_ends+nim,ivf2->im_ends,ivf2->nim*sizeof(*ivf->im_ends));
  
  /* copy pts */  
  memcpy(ivf->points+tot_pts,ivf2->points,ivf2->tot_pts*sizeof(*ivf->points));

  free(ivf2);
}

/******************************************************* querying */

static double normal(double x,double sigma) {
  return exp(-x*x/(sigma*sigma));
}

static double fvec_min_non0(float *tab,int n) {
  float min=1e30;
  int i;
  for(i=0;i<n;i++) 
    if(tab[i]>1e-5 && tab[i]<min) min=tab[i];

  return min;
}

static float fvec_quantile_non0(const float *tab,int n,int q) {
  float *t=fvec_new(n);
  int i,j=0;
  for(i=0;i<n;i++) 
    if(tab[i]>1e-4) t[j++]=tab[i];
  float ret=fvec_quantile(t,j,j*q/n);
  free(t);
  return ret;
}



static int ivfgeo_pq_pt_to_imno(const ivfgeo_pq_t * ivf,
                                 int label) {
  int i0=-1,i1=ivf->nim-1;
  const int *im_ends=ivf->im_ends;
  /* */
  while(i0+1<i1) {
    int imed=(i0+i1+1)/2;

    if(label<im_ends[imed]) i1=imed;
    else                    i0=imed;    
  }
/*  printf("%d -> %d\n",label,i1);*/
  return i1;
}

static ivfgeo_pq_point_t* pointset_to_ivfgeo_pq_points(const pointset_t *ps) {
  ivfgeo_pq_point_t*pts=NEWA(ivfgeo_pq_point_t,ps->n); 
  int i;
  for(i=0;i<ps->n;i++) {
    pts[i].scale=quantize_scale(&ps->pts[i].geom); 
    pts[i].angle=quantize_angle(&ps->pts[i].geom); 
  }  
  return pts;
}


static float fvec_exp_minus_and_sum(float *tab,int n,float sigma) {
  float fm=fvec_max(tab,n);
  double sum=0;
  int i;

  for(i=0;i<n;i++) if(tab[i]>0) 
    sum+=tab[i]=exp(-tab[i]/(fm*sigma));
  
  return sum;
}

static float synth_im_scores(const ivfgeo_pq_t * ivf,int imno,const ivfgeo_pq_point_t *qpts,
                             int nmatch,const int *qnos,const int *bnos,const float *match_scores);


void ivfgeo_pq_query (const ivfgeo_pq_t * ivf,
                      const pointset_t * ps,
                      float *im_dists) {
  
  /* query params */
  int k=ivf->k;
  int ma=ivf->ma;
  double ma_disratio=ivf->ma_disratio;
  int nthread=ivf->nthread;
  int dist_w_type=ivf->dist_w_type;
  double sigma=ivf->sigma;


  fvec_set(im_dists,ivf->nim,0.0);

  if(ps->n==0) {
    return;
  } 

  int d=ivf->ivfpq->pq->d;
  
  float *descs=fvec_new(d*ps->n);

  int d2=pointset_into_fvecs(ps,descs);

  ivfgeo_pq_point_t* qpts=pointset_to_ivfgeo_pq_points(ps);

  assert(d2==d);

  int *labels=ivec_new(k*ps->n);
  float *dists=fvec_new_0(k*ps->n); /* some distances are not filled in so set to 0 */

  int *vw=ivec_new(ma*ps->n);
  int *label_map=ivec_new(ma*ps->n);


  /* do all point queries */

  ivfpq_query_and_get_vw (ivf->ivfpq, descs, ps->n, ma, ma_disratio, k, nthread, 
                          labels, dists, vw, label_map);
  
  /* convert to scores, apply TF-IDF weighting */

  double thr;

  if(dist_w_type==0) 
    thr=sqrt((1+sigma)*fvec_max(dists,k*ps->n));

  int i,j;
  for(i=0;i<ps->n;i++) { /* loop over query points */ 
    float *dists_i=dists+k*i;
    int *vw_i=vw+i*ma;
    int *label_map_i=label_map+i*ma;

    

    switch(dist_w_type) {
    case 1: case 2: 
      thr=sqrt(sigma*fvec_max(dists_i,k));
      break;
    case 3: 
      thr=sqrt(fvec_min_non0(dists_i,k)*sigma);
      break;
    case 4: 
      thr=fvec_max(dists_i,k);
      break;
    case 6: 
      thr=fvec_quantile_non0(dists_i,k,k/4);
      break;
    case 8: 
      thr=fvec_exp_minus_and_sum(dists_i,k,sigma); 
      break;
    }

    int vw_begin=0;
    for(j=0;j<ma;j++) { /* loop over vw of the multiple assignement */ 
      int w=vw_i[j];
      if(w<0) break;

      double tfidf_w=ivf->tfidf_weights[w] * ivf->tfidf_weights[w]; 
      int vw_end=label_map_i[j];
      int l;
      for(l=vw_begin;l<vw_end;l++) {
        double d=dists_i[l];
        double score; 

        switch(dist_w_type) {
        case 0: case 1: 
          score=1-sqrt(d)/thr; 
          break;
        case 2: 
          score=normal(sqrt(d),2*thr);
          break;
        case 4: case 6:
          score=exp(-d/(thr*sigma));
          break;
        case 5: 
          score=1/(d+sigma);
          break;
        case 8: 
          score=d/thr;
          break;
        default: assert(0);
        }

        double score_w=score*tfidf_w;

        dists_i[l]=score_w;
      }
      vw_begin=vw_end;
    }
  }


  /* transpose score tables (qpt, vw, bpt) -> (bpt, qpt) */

  free(label_map); /* not interested in labels any more */
  
  int nmatch;
  int *qpt_map;
  {

    /* TODO implement with merge_ordered_sets */
    int kn=k*ps->n;
    int *perm=ivec_new(kn);

    ivec_sort_index(labels,kn,perm); 

    int skip=0;
    while(skip<kn && labels[perm[skip]]<0) skip++;
    
    nmatch=kn-skip;

    int *labels2=ivec_new(nmatch); 

    for(i=0;i<nmatch;i++) 
      labels2[i]=labels[perm[i+skip]];
    
    free(labels);
    labels=labels2;

    float *dists2=fvec_new(nmatch); 

    for(i=0;i<nmatch;i++) 
      dists2[i]=dists[perm[i+skip]];
    
    free(dists);
    dists=dists2;
    
    qpt_map=ivec_new(nmatch); 
    for(i=0;i<nmatch;i++) 
      qpt_map[i]=perm[i+skip]/k;
    
   
    free(perm);
  } 
  
  int i0=0;
  while(i0<nmatch) {
    
    int imno=ivfgeo_pq_pt_to_imno(ivf,labels[i0]);
    int i1=i0+1;
    while(i1<nmatch && labels[i1]<ivf->im_ends[imno]) i1++;

    /* db points are in i0..i1 */

    im_dists[imno]=synth_im_scores(ivf,imno,qpts,
                                   i1-i0,qpt_map+i0,labels+i0,dists+i0);

    i0=i1;
  }

  free(qpt_map);
  free(qpts);
  free(labels);
  free(dists);

  
  /* per-image normalization */
  
  if (ivf->norm_type != norm_none) {
    distype_t qnorm;

    if (ivf->norm_type == norm_type_1)
      qnorm = vw_tab_norm1 (vw, ps->n * ma);
    else if (ivf->norm_type == norm_type_2)
      qnorm = vw_tab_norm2 (vw, ps->n * ma);
    else if (ivf->norm_type == norm_type_1_sqrt) 
      qnorm = sqrt(vw_tab_norm1 (vw, ps->n * ma));      
    else assert (0);


    assert(qnorm>0);

    for (i = 0; i < ivf->nim; i++) {
      im_dists[i] = im_dists[i] / (ivf->norms[i] * qnorm);
    }
  }

  free(vw);
  
}


static float synth_im_scores(const ivfgeo_pq_t * ivf,int imno,const ivfgeo_pq_point_t *qpts,
                             int nmatch,const int *qnos,const int *bnos,const float *match_scores) {
  float score=0;
  int i;


  const int na = ivfgeo_nb_angles, ns = 2 * ivfgeo_nb_scales - 1;

  float *dis_wgc = ivf->wgc_type == 0 ? NULL : fvec_new_0(na+ns);  
  
  for(i=0;i<nmatch;i++) {
    
    const ivfgeo_pq_point_t *qpt=qpts+qnos[i];
    const ivfgeo_pq_point_t *bpt=ivf->points+bnos[i];
    
    int ref_scale = MIN(qpt->scale,bpt->scale);

    double w3 = match_scores[i] * ivf->scale_weights[ref_scale];    

    if(!ivf->wgc_type) 
      score+=w3;
    else {
      int da = (bpt->angle - qpt->angle) & (ivfgeo_nb_angles - 1);
      int ds = bpt->scale - qpt->scale + ivfgeo_nb_scales - 1;
      dis_wgc[da] += w3;
      dis_wgc[na + ds] += w3;          
    }
  }

  if(ivf->wgc_type) 
    score=synth_wgc_tab (dis_wgc,ivf->wgc_type,ivf->wgc_weights); 


  return score;
}


/******************************************************* I/O */


#define WRITEANDCHECK(a,n) if(fwrite(a,sizeof(*a),n,f)!=n) {perror("ifvpq_fwrite"); abort(); }



void ivfgeo_pq_fwrite(const ivfgeo_pq_t * ivf, FILE *f) {

  ivfpq_fwrite(f,ivf->ivfpq);

  WRITEANDCHECK(&ivf->nim,1);
  WRITEANDCHECK(&ivf->tot_pts,1);

  WRITEANDCHECK(&ivf->norm_type,1);
  WRITEANDCHECK(&ivf->tfidf,1);

  WRITEANDCHECK(ivf->norms,ivf->nim);
  WRITEANDCHECK(ivf->im_ends,ivf->nim);
  WRITEANDCHECK(ivf->points,ivf->tot_pts);

}

#undef WRITEANDCHECK

#define READANDCHECK(a,n) if(fread(a,sizeof(*a),n,f)!=n) {perror("ifvpq_read"); abort(); }


ivfgeo_pq_t *ivfgeo_pq_fread(FILE *f,float *centroids) {

  ivfpq_t *ivfpq=ivfpq_fread(f,centroids);
  
  ivfgeo_pq_t *ivf=NEWAC(ivfgeo_pq_t,1);

  ivf->ivfpq=ivfpq;

  int nim; 
  long tot_pts;



  READANDCHECK(&nim,1);
  READANDCHECK(&tot_pts,1);

  ivfgeo_pq_add_pts(ivf,nim,tot_pts); 

  READANDCHECK(&ivf->norm_type,1);
  READANDCHECK(&ivf->tfidf,1);
   
  READANDCHECK(ivf->norms,nim);
  READANDCHECK(ivf->im_ends,nim);
  READANDCHECK(ivf->points,tot_pts);

  /* decent defaults for query params */
  
  ivf->ma=1;
  ivf->ma_disratio=1.2;
  ivf->k=10;
  ivf->nthread=1; 

  ivfgeo_pq_compute_tfidf(ivf);
  ivfgeo_pq_compute_scale_weights(ivf);  
  ivfgeo_pq_set_wgc_type(ivf,0);


  return ivf;
}
