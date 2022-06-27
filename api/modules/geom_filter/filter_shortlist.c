
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>

#include <pthread.h>

#include <yael/vector.h>
#include <yael/nn.h>
#include <utils/generic_hash_table.h>
#include <utils/fast_alloc.h>
#include <siftgeo/siftgeo_and_nn.h>

#include "filter_shortlist.h"
#include "geometry.h"


/* for compatibility with C++ version */
#define DESC_BYTE_TO_FLOAT(d) (((d)+0.5)/512.0)


/* #define DESC_BYTE_TO_FLOAT(d) ((float)(d)) */ 


#define NEWA(type,n) (type*)malloc(sizeof(type)*(n))
#define NEWAC(type,n) (type*)calloc(sizeof(type),(n))
#define NEW(type) NEWA(type,1)

#define MAX(a,b) ((a)>(b) ? (a) : (b))
#define MIN(a,b) ((a)<(b) ? (a) : (b))




/*-----------------------------------------------------------------
 * shortlist functions
 *-----------------------------------------------------------------*/


typedef struct {
  int npt;
  int ptofs;
  imdesc_t *imdesc;
} shortlist_elt_t;

struct shortlist_t {
  /* set in shortlist_new */

  int nim;
  int npt,d;
  
  /* find image number that corresponds to a point in 0..npt 
   * could be done by bissection, but who cares about memory nowadays?
   */
  int *map_pt_to_im;

  shortlist_elt_t *elts;
  
  pthread_mutex_t fpoints_lock;

  float *fpoints; /* cached data matrix */

};



shortlist_t* shortlist_new(imdesc_t **imdescs,int nim) {
  shortlist_t *sl=NEW(shortlist_t);

  sl->nim=nim;

  sl->elts=NEWA(shortlist_elt_t,nim);
  
  int totnpt=0;
  {  
    int i;
    for(i=0;i<nim;i++) {
      int npt=imdescs[i]->n;
      sl->elts[i].imdesc=imdescs[i];
      if(npt!=0) 
        sl->d=imdescs[i]->pts[0].dim;      
      sl->elts[i].npt=npt;
      sl->elts[i].ptofs=totnpt;      
      totnpt+=npt;
    }
  }
  if(totnpt==0) {
    fprintf(stderr,"shortlist_new: no points to index nim=%d!!\n",nim);
    /* should dealloc... */
    return NULL;
  }
  sl->npt=totnpt;
    
  {
    sl->map_pt_to_im=NEWA(int,totnpt);
    int *wp=sl->map_pt_to_im;
    int i,j;

    for(i=0;i<nim;i++) 
      for(j=0;j<imdescs[i]->n;j++) 
        *wp++=i;
    assert(sl->map_pt_to_im+totnpt==wp);
  }

  sl->fpoints=NULL;
  pthread_mutex_init(&sl->fpoints_lock,NULL);

  return sl;
}

/* delete functions are cascaded: if data for a stage is deleted 
 * data for all subsequent stages is deleted as well.
 */

void shortlist_delete(shortlist_t *sl) {

  free(sl->map_pt_to_im);

  free(sl->elts);
  free(sl->fpoints);
  
  free(sl);
    
}


/*! free n imagematch_t's */
void imagematches_delete(imagematch_t *imms,int n) {
  int i;
  for(i=0;i<n;i++)
    fa_free(imms[i].fa_pointmatch);
  free(imms);
}


imagematch_t* imagematches_new(int n) {
  imagematch_t *imms=NEWAC(imagematch_t,n);
  int i;
  for(i=0;i<n;i++)
    imms[i].fa_pointmatch=fa_alloc(sizeof(pointmatch_t),512);
  return imms;
}


/*
static double dist2_u8(unsigned char *a,unsigned char *b,int n) {
  double accu=0;
  while(n--) {
    double d=DESC_BYTE_TO_FLOAT(a[n])-DESC_BYTE_TO_FLOAT(b[n]);
    accu+=d*d;
  }
  return accu;
}
*/



static void descs_convert(float *tab,int n) {
  while(n--) tab[n]=DESC_BYTE_TO_FLOAT(tab[n]);
}


// match points for single image
static void imm_addmatch(imagematch_t *imm,int qi,int dbi,double d2,double thr) {
  assert(d2>=0);
  double score=1-sqrt(d2)/thr;
          
  imm->stage2_nmatch++;
  imm->stage2_votes+=score;
    
  pointmatch_t *pm=fa_allocunit(imm->fa_pointmatch);
  
  pm->score=score;
  pm->qpt=qi;
  pm->dbpt=dbi;
  
  pm->next=imm->ptmatches;
  imm->ptmatches=pm;      

}

// match points for short list images
imagematch_t *shortlist_match_points_exact(shortlist_t *sl,imdesc_t *qim,
                                           int method,double thr,
                                           void (*peek_fun)(void *arg,double frac),
                                           void *peek_arg) {

  int k=(int)thr; /* k nearest neighbours */

  imagematch_t *imms=imagematches_new(sl->nim);


  int d,i,j;
  float *qpts=siftgeo_to_fvecs(qim->pts,qim->n,&d);

  assert(d==sl->d);

  descs_convert(qpts,d*qim->n);

  /* put evrything in a big matrix */

  

  float *dbpts=NULL;
  
  pthread_mutex_lock(&sl->fpoints_lock);
  
  if(!sl->fpoints) {

    sl->fpoints=dbpts=fvec_new(sl->npt*d);

    for(i=0;i<sl->nim;i++)
      pointset_into_fvecs(sl->elts[i].imdesc,dbpts+sl->elts[i].ptofs*d);
    
    descs_convert(dbpts,d*sl->npt);
  }
  pthread_mutex_unlock(&sl->fpoints_lock);

  if(method==2) {    
    int *indices=NEWA(int,k*qim->n);

    float *dist2=quantize_codebook_multiple(qim->n,sl->npt,d,k,
                                            dbpts,qpts,indices,NULL,NULL);
    
    thr=0;    
    for(i=0;i<k*qim->n;i++) 
      if(dist2[i]>thr) thr=dist2[i];
    thr=sqrt(thr*1.1);

    for(i=0;i<qim->n;i++) {
      float *disti=dist2+i*k;
      int *indi=indices+i*k;

      for(j=0;j<k;j++) {
        int imno=sl->map_pt_to_im[indi[j]];
        assert(0<=imno && imno<sl->nim);

        imagematch_t *imm=&imms[imno];
        int ptno=indi[j]-sl->elts[imno].ptofs;
        assert(0<=ptno && ptno<sl->elts[imno].npt);

        imm_addmatch(imm,i,ptno,disti[j],thr);

      }

    }

    free(dist2);
    free(indices);

  } else 
    assert(!"not implemented");
  
  free(qpts);
  return imms;
}


static int binsign_hamming (binsign_t bs1, binsign_t bs2);

/*
 * find common VWs between qim and dbim, which should be sorted by vw
 */
void imagematch_align_vw(imdesc_t *qim,
                     imdesc_t *dbim,
                     imagematch_t *imm,
                     int hamming_thresh) {
  if(qim->n==0 || dbim->n==0) return;

  int vw=-1;
  int j0=0,j1=0; 
  int i,j;
  int prev_qvw=-1;

  for(i=0;i<qim->n;i++) {    

    int qvw=qim->pts[i].vw;

    assert(qvw>=prev_qvw);
    prev_qvw=qvw;

    while(qvw>vw) {
      /* step db */
      j0=j1;
      if(j0>=dbim->n) goto stoploop;
      vw=dbim->pts[j0].vw;
      for(j1=j0+1;j1<dbim->n;j1++) {
        assert(dbim->pts[j1].vw>=vw);
        if(dbim->pts[j1].vw!=vw) break;
      }
    } 
    if(qvw<vw) {
      /* step q */
      continue;
    } 
    
    for(j=j0;j<j1;j++) {

      point_t *pt=&dbim->pts[j];
    
      int hd=binsign_hamming(pt->binsign,qim->pts[i].binsign);
      if(!(hd<=hamming_thresh)) continue;
      
      double score=1-hd/(double)(hamming_thresh+1);
        
      imm->stage2_nmatch++;      
      imm->stage2_votes+=score;
      
      pointmatch_t *pm=fa_allocunit(imm->fa_pointmatch);
      
      pm->score=score;
      pm->qpt=i;
      pm->dbpt=j;
      
      pm->next=imm->ptmatches;
      imm->ptmatches=pm;      
      
    }
  }
 stoploop:;
}


imagematch_t *shortlist_match_points_vw(shortlist_t *sl,imdesc_t *query,
                                        int hamming_thresh,
                                        void (*peek_fun)(void *arg,double frac),
                                        void *peek_arg) {

  int i;

  /* find matches */
  imagematch_t *imms=imagematches_new(sl->nim);

  for(i=0;i<sl->nim;i++) {

    imagematch_align_vw(query,
             sl->elts[i].imdesc,&imms[i],
             hamming_thresh);

    if(peek_fun) 
      (*peek_fun)(peek_arg,i/(double)sl->nim);
  }
  return imms;
}


/*-----------------------------------------------------------------
 * geometric verification 1: the Hough table
 *-----------------------------------------------------------------*/



#define LHDIM 4




static void compute_hough_coords(geom_t *qpt,
                                 geom_t *dbpt,
                                 double h[LHDIM]) {

  double scale_ratio=dbpt->scale/qpt->scale;
  double angle_diff=dbpt->angle-qpt->angle;

  h[0]=log(scale_ratio);
  h[1]=angle_diff;

  double c=cos(angle_diff),s=sin(angle_diff);

  h[2]=dbpt->x+scale_ratio*(-qpt->x*c+qpt->y*s);
  h[3]=dbpt->y+scale_ratio*(-qpt->x*s-qpt->y*c);

}


static int params_1aff(double i0,double o0,
                       double i1,double o1,
                       double h[2]) {
  double a=(o1-o0)/(i1-i0);

  if(!finite(a)) return 0;

  h[0]=log(a);
  h[1]=o0-a*i1;

  return 1;
}


static int compute_hough_coords_2aff(geom_t *q0,geom_t *db0,
                                     geom_t *q1,geom_t *db1,
                                     double h[LHDIM]) {
  
  if(!params_1aff(q0->x,db0->x,q1->x,db1->x,h)) 
    return 0;

  if(!params_1aff(q0->y,db0->y,q1->y,db1->y,h+2)) 
    return 0;
  
  return 1;
}



static long hash_lh_cell(void *key) {
  unsigned int *k=(unsigned int*)key;
  unsigned int b    = 378551;
  unsigned int a    = 63689;
  unsigned int hash = 0;
  int i;
  for(i = 0; i < LHDIM; i++) {
    hash = hash * a + k[i];
    a    = a * b;
  }
  return hash;
} 


static int cmp_lh_cell(void *key0,void *key1) {
  return memcmp(key0,key1,sizeof(int)*LHDIM);
}


typedef struct {
  int hi[LHDIM];
  float vote;
  int n,na;
  pointmatch_t *matches[0];
} hough_bin_t;



/* insertion is recursive because of soft-assign: 
 * there are 2^LHDIM cells to assign
 * */
static void hough_table_insert(hash_table_t *hough_table,
                               pointmatch_t *pm,
                               pointmatch_t *pm2,
                               double h[LHDIM],
                               int hi[LHDIM],
                               double vote,int dim) {
  if(dim<LHDIM) {    
    /* insert in upper and lower bound of bin */
    double hf=floor(h[dim]);
    double dh=h[dim]-hf;

    hi[dim]=(int)hf;
    if(dh!=1)
      hough_table_insert(hough_table,pm,pm2,h,hi,vote*(1-dh),dim+1);

    hi[dim]++;
    if(dh!=0)
      hough_table_insert(hough_table,pm,pm2,h,hi,vote*dh,dim+1);

  } else { /* we can add the match to the cell */

    int npm=pm2 ? 2 : 1;

    hough_bin_t *bin=hash_table_lookup(hough_table,hi);

    if(!bin) { /* allocate a new bin */
      int init_alloc=8;
      bin=malloc(sizeof(hough_bin_t)+sizeof(pointmatch_t*)*init_alloc);
      memcpy(bin->hi,hi,sizeof(int)*LHDIM);
      bin->vote=0;
      bin->n=0;
      bin->na=init_alloc;
      hash_table_insert(hough_table,bin->hi,bin);
    } else if(bin->n+npm>=bin->na) {      
      /* !! cannot do realloc because mem still referenced by hashtable */
      hough_bin_t *prev_bin=bin;      
      bin=malloc(sizeof(hough_bin_t)+sizeof(pointmatch_t*)*prev_bin->na*2);
      memcpy(bin,prev_bin,
             sizeof(hough_bin_t)+sizeof(pointmatch_t*)*prev_bin->na);
      bin->na*=2;      
      hash_table_insert(hough_table,bin->hi,bin);
      free(prev_bin);
    } 
    
    bin->vote+=vote;
    bin->matches[bin->n++]=pm;
    if(pm2) 
      bin->matches[bin->n++]=pm2; 

  }

}



static void hough_table_to_array(hash_table_t *hough_table,
                                 hough_bin_t ***bins_out,int *nbin_out) {
  int nbin;
  
  hash_table_stats(hough_table,NULL,&nbin);
  
  hough_bin_t **bins=NEWA(hough_bin_t*,nbin);
  
  hash_table_it_t *it;
  int i=0;
  for(it=hash_table_iterate(hough_table,NULL); it;
      it=hash_table_iterate(hough_table,it)) {
    hough_bin_t *bin=it->val;
    bins[i++]=bin;
  }
  assert(i==nbin);

  *nbin_out=nbin;
  *bins_out=bins;
  
}

static void build_hough_transform(point_t *qimdesc,point_t *dbimdesc,
                                  pointmatch_t *match,
                                  hough_bin_t ***bins_out,int *nbin_out,
                                  double bin_sizes[LHDIM]) {
  
  /* maps int[LHDIM] to hough_bin_t* (key is included in value) */
  hash_table_t *hough_table=hash_table_init(hash_lh_cell,cmp_lh_cell);
  pointmatch_t *it;

  for(it=match;it;it=it->next) {

    /* compute hough coordinates */
    double h[LHDIM];

    compute_hough_coords(&qimdesc[it->qpt].geom,
                         &dbimdesc[it->dbpt].geom,
                         h);
    /* scale coords */
    int i;
    for(i=0;i<LHDIM;i++) 
      h[i]/=bin_sizes[i];
                       
    double vote=it->score;

    assert(vote>0);

    int hi[LHDIM];
    hough_table_insert(hough_table,
                       it,NULL,h,hi,
                       vote,0);

  }
    
  hough_table_to_array(hough_table,bins_out,nbin_out);

  hash_table_delete(hough_table);

}



//static int count_ptmatches(pointmatch_t *pm) {
int count_ptmatches(pointmatch_t *pm) {
  int n=0;
  for(;pm;pm=pm->next) n++;
  return n;
}


static double geom_dist2(geom_t *g1,geom_t *g2) {
  double dx=g1->x-g2->x;
  double dy=g1->y-g2->y;
  return dx*dx+dy*dy;
}


static int cmp_pointmatch(const void *av, const void *bv) {
  pointmatch_t *a=*(pointmatch_t**)av;
  pointmatch_t *b=*(pointmatch_t**)bv;
  return a->score<b->score ? 1 : -1;
}

#define BIGPRIME 100003

/* build a Hough transform of a model */

static void build_hough_transform_2aff
    (point_t *qimdesc,point_t *dbimdesc,
     pointmatch_t *match,
     hough_bin_t ***bins_out,int *nbin_out,
     double bin_sizes[LHDIM],
     double distinct_tolerance,
     int max_match_per_i_plus_j) {
  
  /* maps int[LHDIM] to hough_bin_t* (key is included in value) */
  hash_table_t *hough_table=hash_table_init(hash_lh_cell,cmp_lh_cell);

  int nmatch=count_ptmatches(match);

  /* matches ordered by decreasing score */
  pointmatch_t **ordered_ptmatches=NEWA(pointmatch_t *,nmatch);
  { 
    int i;
    pointmatch_t *it;
    
    for(it=match,i=0;it;it=it->next,i++) ordered_ptmatches[i]=it;
    assert(i==nmatch);
    
    qsort(ordered_ptmatches,nmatch,sizeof(pointmatch_t *),cmp_pointmatch);
  }

  int i_plus_j;

  double dtol2=distinct_tolerance*distinct_tolerance;

  /* we want pairs (i,j) \in [0..nmatch-1] of matches, but there are
    too many of them. So we order the pairs by relevance and
    (randomly) take max_match_per_i_plus_j different pairs for every
    i+j. This way, there are fewer pairs for bigger i's and j's.
  */

  for(i_plus_j=1;i_plus_j<=2*nmatch-3;i_plus_j++) {
    int i_min=i_plus_j/2+1;    

    int i=i_min;
    int n;

    for(n=0;n<=i_plus_j-i_min && n<max_match_per_i_plus_j;n++) {

      /* step to next element */
      i=i_min+(i-i_min+BIGPRIME)%(i_plus_j-i_min+1);
      
      if(i>=nmatch) continue;

      int j=i_plus_j-i;

      pointmatch_t *iti=ordered_ptmatches[i];
      pointmatch_t *itj=ordered_ptmatches[j];

      geom_t *qi=&qimdesc[iti->qpt].geom;
      geom_t *dbi=&dbimdesc[iti->dbpt].geom;

      geom_t *qj=&qimdesc[itj->qpt].geom;
      geom_t *dbj=&dbimdesc[itj->dbpt].geom;
      
      /* check distance between src and dest pts */
      
      if(geom_dist2(qi,qj)<dtol2 ||
         geom_dist2(dbi,dbj)<dtol2) continue;
      
      /* compute hough coordinates */
      double h[LHDIM];
      
      int ret=compute_hough_coords_2aff(qi,dbi,qj,dbj,h);
      if(!ret) continue;

      /* scale coords */
      int i;
      for(i=0;i<LHDIM;i++) 
        h[i]/=bin_sizes[i];
      
      double vote=iti->score+itj->score;
      
      assert(vote>0);
      
      int hi[LHDIM];
      hough_table_insert(hough_table,iti,itj,h,hi,vote,0);
      
    }

  }

    
  hough_table_to_array(hough_table,bins_out,nbin_out);

  hash_table_delete(hough_table);
  
  free(ordered_ptmatches);
}
                           

/*-----------------------------------------------------------------
 * geometric verification 2: estimation
 *-----------------------------------------------------------------*/

static int cmp_bins(const void *av, const void *bv) {
  hough_bin_t *a=*(hough_bin_t**)av;
  hough_bin_t *b=*(hough_bin_t**)bv;
  return a->vote<b->vote ? 1 : -1;
}

static void order_bins(hough_bin_t **bins,int nbin) {
  qsort(bins,nbin,sizeof(hough_bin_t*),&cmp_bins);  
}

static double compute_votes_subset(pointmatch_t **pms,int npm) {
  double accu=0;
  int i;
  for(i=0;i<npm;i++) accu+=pms[i]->score;

  return accu;
}

/*
eigenvals of [a c ; c b]
*/

static double sqr(double x) {return x*x; }

static void eigen_values(double a,double c,double b,
                         double *ev1_out,double *ev2_out) {
  double sdelta=sqrt(sqr(a-b)+4*c*c);

  *ev1_out=0.5*(a+b-sdelta);
  *ev2_out=0.5*(a+b+sdelta);

}


static int singular_values(double a11,double a12,double a21,double a22,
                            double *sv1_out,double *sv2_out) {
  double ev1,ev2;

  eigen_values(a11*a11+a12*a12,
               a11*a21+a12*a22,a22*a22+a21*a21,
               &ev1,&ev2);

  /* happens with numerical instabilities */
  if(!(ev1>=0 && ev2>=0))
    return 0;

  *sv1_out=sqrt(ev1);
  *sv2_out=sqrt(ev2);

  return 1;
}



static int aff_acceptable(double aff[6],
                          lowehough_parameters_t *params) {
  
  
  double det=fabs(aff[0]*aff[4]-aff[3]*aff[1]);
 
  if(det>params->max_det || det<1/params->max_det) 
    return -1;
  
  double sv1,sv2;

  if(!singular_values(aff[0],aff[1],aff[3],aff[4],&sv1,&sv2)) 
    return -4;  
  
  /* aspect ratio of an ellipse transformed by aff */
  double ar = (sv1 < sv2) ? sv1/sv2 : sv2/sv1;

  if(ar<params->min_ar) 
    return -2;

  return 1;
}


static int aff_bin_acceptable(double aff[6],
                              int *bin,double *bin_sizes,
                              lowehough_parameters_t *params) {

  int ret=aff_acceptable(aff,params);

  if(ret<0) return ret;

  /* translation close enough to bin? */
  
  if(!(fabs(bin[2]-aff[2]/bin_sizes[2])<params->max_bin_deviation &&
       fabs(bin[3]-aff[5]/bin_sizes[3])<params->max_bin_deviation))
    return -3;

  return 1;
}



/* coding of a 2*2 transformation matrix that attempts to produce a perceptual
 * notion of deformation */
typedef struct {
  /* deformations are in order of increasing (perceptual) badness */
  double scale;  /* scaling */
  double ar;     /* aspect ratio change */
  double a1;     /* rotation angle of horizontal axis */
  double a12;    /* difference of rotation of the vertical axis */
} deformation_t;


#if 0
/* to avoid "unused" warning */

static void deformation_to_aff(deformation_t *df,double aff[6]) {
  /* deformation is the composition of:
   * 1. scaling of factor df->scale
   * 2. anisotropic scaling of sqrt(df->ar) horizontally and 1/sqrt(df->ar)
   *    vertically
   * 3. rotation of angle a1 horizonatally and (a1+a12) vertically
   */

  double f0=df->scale*sqrt(df->ar);
  double f1=df->scale/sqrt(df->ar);

  aff[0]=f0*cos(df->a1); aff[1]=-f0*sin(df->a1+df->a12); 
  aff[3]=f1*sin(df->a1); aff[4]=f1*cos(df->a1+df->a12); 
  
}

#endif 

static const double overflow_double=1e10,underflow_double=1e-10;

#define DECENT_DOUBLE(d) ((fabs(d)<overflow_double && fabs(d)>underflow_double))


static int aff_to_deformation(double aff[6],deformation_t *df) {
  
  double ar2=(aff[1]*aff[1]-aff[0]*aff[0])/(aff[3]*aff[3]-aff[4]*aff[4]);
  
  if(!DECENT_DOUBLE(ar2)) return -1;

  df->ar=sqrt(ar2);

  double f1=sqrt(aff[0]*aff[0]/ar2+aff[3]*aff[3]);
  double f0=f1*df->ar;

  df->scale=f1*sqrt(df->ar);

  if(!DECENT_DOUBLE(df->scale)) return -4;

  if(aff[3]==0 && aff[0]==0) return -2;
  if(aff[1]==0 && aff[4]==0) return -3;

  df->a1=atan2(aff[3]/f1,aff[0]/f0);
  df->a12=atan2(-aff[1]/f0,aff[4]/f1)-df->a1;

  return 0;
}

static double gaussian_weight(double sigma,double x) {
  if(sigma<0) return 1.0;
  assert(finite(x));
  return exp(-x*x/(sigma*sigma));  
}


static double deformation_weight(lowehough_parameters_t *params,
                                 double aff[6]) {
  if(!params->weight_deformation) 
    return 1.0;
  
  deformation_t df;
  
  int cv_ret=aff_to_deformation(aff,&df);

  if(cv_ret<0) return 0.0;

  double a1=df.a1;
  if(params->weight_portrait) {
    /* find closest 1/4 turn to a1 */
    double da=fabs(a1);
    if(fabs(a1-M_PI/2)<da) 
      {a1=a1-M_PI/2; da=fabs(a1-M_PI/2); }
    if(fabs(a1+M_PI/2)<da) 
      a1=a1+M_PI/2;
  }

  return 
    gaussian_weight(params->sigma_logscale,log(df.scale))*
    gaussian_weight(params->sigma_logar,log(df.ar))*
    gaussian_weight(params->sigma_a1,a1)*
    gaussian_weight(params->sigma_a12,df.a12);
}


void lowehough_parameters_default(lowehough_parameters_t *params) {

  params->bin_sizes[0]=log(3.0); /* scale */
  params->bin_sizes[1]=M_PI/3.0; /* angle */
  params->bin_sizes[2]=13.0;     /* translation x */
  params->bin_sizes[3]=13.0;     /* translation y */

  params->max_nbin=1024;
  params->pos_error_in_affine=sqrt(15);
  params->min_match_after_affine=4;
  params->min_match_before_affine=3;
  params->distinct_tolerance=sqrt(5);

  params->max_det=1e4;
  params->min_ar=sqrt(0.15);
  params->max_bin_deviation=3.0;

  params->weight_deformation=0;
  params->sigma_logscale=1.0; /* weight for scale 0.5 = 0.8 */
  params->sigma_logar=0.4;    /* weight for ar 0.75 = 0.6  */
  params->sigma_a1=0.4;
  params->sigma_a12=0.1;
  params->weight_portrait=1;  
  params->min_weight=0.05;
  
  params->verbose=0;

  params->use_2aff_model=0;

  params->bin_sizes_2aff[0]=log(1.5); /* scale x */
  params->bin_sizes_2aff[1]=13.0;     /* translation x */
  params->bin_sizes_2aff[2]=log(1.5); /* scale y */
  params->bin_sizes_2aff[3]=13.0;     /* translation y */

  memset(&params->delta_a_stats,0,sizeof(params->delta_a_stats));
  
}


// ----------------------------------------------------------------

// key geom verif function
//filter matched points
void match_lowehough(imdesc_t *qimdesc,imdesc_t *dbimdesc,
                     imagematch_t *match,
                     lowehough_parameters_t *params) {

  
  
  lowehough_parameters_t default_params;
  if(!params) {
    lowehough_parameters_default(&default_params);
    params=&default_params;
  }

  int verbose=params->verbose;


  hough_bin_t **bins;
  int nbin;
  int nmatch=count_ptmatches(match->ptmatches);    

  if(!params->use_2aff_model) {

    build_hough_transform(qimdesc->pts,dbimdesc->pts,match->ptmatches,
                          &bins,&nbin,
                          params->bin_sizes);
  } else {

    build_hough_transform_2aff(qimdesc->pts,dbimdesc->pts,match->ptmatches,
                               &bins,&nbin,
                               params->bin_sizes_2aff,
                               params->distinct_tolerance,
                               10);

  }

  order_bins(bins,nbin); 
  
  if(verbose>1)
    printf("  db image of %d pts -> %d matches -> %d Hough bins\n",
           dbimdesc->n,nmatch,nbin);

  typedef struct {
    int n_agree;
    pointmatch_t **agreeing_matches;    
    double votes;    
    double raw_votes;
    double aff[6];
    int* ids; //added, indexes of matched and passed stage 3 descs 
  } aff_hypothesis_t;

  aff_hypothesis_t best_hyp;

  memset(&best_hyp,0,sizeof(best_hyp));
  
  pointmatch_t **agreeing_matches_0=NEWA(pointmatch_t*,nmatch);
  pointmatch_t **agreeing_matches_1=NEWA(pointmatch_t*,nmatch);
  
  int b;
  for(b=0;b<nbin && b<params->max_nbin;b++) {
        hough_bin_t *bin=bins[b];
        
        aff_hypothesis_t hyp;                
        memset(&hyp,0,sizeof(hyp));
        hyp.ids = ivec_new_0(nmatch);
        
        hyp.agreeing_matches=
          best_hyp.agreeing_matches==agreeing_matches_0 ?
          agreeing_matches_1 : agreeing_matches_0;

        if(verbose>2) 
          printf("    Bin %d [%d %d %d %d] (%d/%d pts, vote=%g)\n",
                 b,bin->hi[0],bin->hi[1],bin->hi[2],bin->hi[3],
                 bin->n,params->min_match_before_affine,bin->vote);

        if(bin->n<params->min_match_before_affine) 
          continue;    

        /* estimate affine transform from point matches in bin */
        
        int aff_ret=estimate_affine_transform(qimdesc->pts,dbimdesc->pts,
                                              bin->matches,bin->n,
                                              hyp.aff);
        if(verbose>2) 
          printf("    affine estimation returns %d\n",aff_ret);
        
        if(aff_ret<0) continue;      

        /* find matches that "agree" with transform */            
        hyp.n_agree=find_agreeing_matches(qimdesc->pts,dbimdesc->pts,
                                          hyp.aff,
                                          match->ptmatches,
                                          params->pos_error_in_affine,
                                          hyp.agreeing_matches,hyp.ids);
        assert(hyp.n_agree<=nmatch);
        if(verbose>2) {
              printf("      affine=[%g %g %g ; %g %g %g]\n",
                     hyp.aff[0],hyp.aff[1],hyp.aff[2],hyp.aff[3],hyp.aff[4],hyp.aff[5]);
              printf("      n_agree=%d/%d (best so far=%d)\n",
                     hyp.n_agree,params->min_match_after_affine,best_hyp.n_agree);
              printf("ids:");
              ivec_print(hyp.ids, nmatch);
        }
        if(hyp.n_agree<params->min_match_after_affine /*  || 
           hyp.n_agree<=best_hyp.n_agree */) 
          continue;
        if(!params->weight_deformation) {
          int acceptable=aff_bin_acceptable(hyp.aff,bin->hi,params->bin_sizes,params);
          if(verbose>2) 
            printf("      acceptable=%d\n",acceptable);      
          if(acceptable<0) 
            continue;    

        }
        /* count nb of distinct query & db points */
        int n_distinct_q,n_distinct_db;
        n_distinct_q=count_distinct_points(qimdesc->pts,
                                           hyp.agreeing_matches,hyp.n_agree,1,
                                           params->min_match_after_affine,
                                           params->distinct_tolerance);
        n_distinct_db=count_distinct_points(dbimdesc->pts,
                                            hyp.agreeing_matches,hyp.n_agree,0,
                                            params->min_match_after_affine,
                                            params->distinct_tolerance);
        if(verbose>2) 
              printf("      n_distinct=%d,%d/%d\n",
                     n_distinct_q,n_distinct_db,params->min_match_after_affine);    
        if(!(n_distinct_db>=params->min_match_after_affine &&
             n_distinct_q>=params->min_match_after_affine)) /* are there enough? */
          continue;     
        hyp.raw_votes=compute_votes_subset(hyp.agreeing_matches,hyp.n_agree);
        double weight=deformation_weight(params,hyp.aff);
        if(verbose>2) 
             printf("      deformation weight %g\n",weight);     
        if(weight<params->min_weight) weight=0;
        hyp.votes=hyp.raw_votes*weight;
        if(verbose>2) 
              printf("      votes=%g (best=%g)\n",
                     hyp.votes,best_hyp.votes);     
        if(hyp.votes<best_hyp.votes){ 
          free(hyp.ids);   
          continue; //return to top of the loop
        }
        /* update best estimate so far */    
        best_hyp=hyp;
        if(verbose>2) 
            printf("      keep\n");    
            
  }// end for bins

  /* deallocate  */
  for(b=0;b<nbin;b++) free(bins[b]);
  free(bins);
  
  /* fill in fields of *match */
  if(verbose>1) {
    printf("  Result: best n_agree=%d ",best_hyp.n_agree);
    if(best_hyp.n_agree>0) {
      printf("votes=%g affine=[%g %g %g ; %g %g %g]\n",
             best_hyp.votes,
             best_hyp.aff[0],best_hyp.aff[1],best_hyp.aff[2],
             best_hyp.aff[3],best_hyp.aff[4],best_hyp.aff[5]);
    } else printf(" (no image match)\n");
    
  }
  
  match->stage3_nmatch=best_hyp.n_agree;

  if(best_hyp.n_agree>0) {
    memcpy(match->affine,best_hyp.aff,sizeof(double)*6);
    match->stage3_votes=best_hyp.votes;
    /* set non-agreeing entries in match to 0 */
    int i;
    for(i=0;i<best_hyp.n_agree;i++) { /* mark correct matches with a score<0 */
      pointmatch_t *pm=best_hyp.agreeing_matches[i]; //agreeing_matches contains pointers to win interest points
      pm->score=-pm->score;
    }
    pointmatch_t *pm;
    for(pm=match->ptmatches;pm;pm=pm->next) { // reset scores of pos and 0 negs 
      pm->score=pm->score<0 ? -pm->score : 0.0;
    } 

    /* collect stats */
    if(params->delta_a_stats.bins || params->delta_a_stats.hamming_dist_bins) {
      for(i=0;i<best_hyp.n_agree;i++) {
        pointmatch_t *pm=best_hyp.agreeing_matches[i];
        point_t *qpt=&qimdesc->pts[pm->qpt];
        point_t *dbpt=&dbimdesc->pts[pm->dbpt];
        

        if(params->delta_a_stats.bins) {
          double delta_a=qpt->geom.angle-dbpt->geom.angle;
          
          delta_a/=2*M_PI;
          delta_a-=floor(delta_a);
          
          int b=(int)(delta_a*params->delta_a_stats.nbin);
          
          params->delta_a_stats.bins[b]+=1.0; /* pm->score; */
        } else {
          binsign_t mask=(binsign_t)(-1UL)<< (64-params->delta_a_stats.hamming_dist_nbit);
          params->delta_a_stats.hamming_dist_bins[binsign_hamming(qpt->binsign&mask,dbpt->binsign&mask)]++;
        }
      }

    }
        
  } else {
    pointmatch_t *pm;
    for(pm=match->ptmatches;pm;pm=pm->next) 
      pm->score=0;   
  }

  free(agreeing_matches_0);  
  free(agreeing_matches_1);  
    
}

void shortlist_filter_lowehough(shortlist_t *sl,imdesc_t *query,imagematch_t *imms,
                                lowehough_parameters_t *params,
                                void (*peek_fun)(void *arg,double frac),
                                void *peek_arg) {

  int verbose=params->verbose;
  int i;
  if(verbose>0) 
    printf("shortlist_filter_lowehough: query image (%d pts), %d db images\n",
           query->n,sl->nim);
  for(i=0;i<sl->nim;i++) {    
    if(verbose>1) 
      printf("image %d\n",i);

    match_lowehough(query,sl->elts[i].imdesc,&imms[i],params);

    if(peek_fun) 
      (*peek_fun)(peek_arg,i/(double)sl->nim);   
  }

}

#if 0
/* to avoid "unused" warning */

static void compute_similarity(geom_t *qpt,geom_t *dbpt,
                               double aff[6]) {
  double da=dbpt->angle-qpt->angle;
  double rs=dbpt->scale/qpt->scale;
  double c=cos(da),s=sin(da);  

  aff[0]=c*rs; aff[1]=-s*rs; 
  aff[3]=s*rs; aff[4]=c*rs; 

  aff[2]=dbpt->x-(aff[0]*qpt->x+aff[1]*qpt->y);
  aff[5]=dbpt->y-(aff[3]*qpt->x+aff[4]*qpt->y);

}

#endif

/*************************************************************************************
 * image pairs
 */






image_pairs_t *image_pairs_new() {

  image_pairs_t *ip=NEWAC(image_pairs_t,1);
  ip->im0=pointset_new();
  ip->im1=pointset_new();
  ip->imm.fa_pointmatch=fa_alloc(sizeof(pointmatch_t),1024);
  ip->imm.ptmatches=NULL;


  return ip;
}


void image_pairs_dump(image_pairs_t *ip,char *fname) {
  FILE *f=fopen(fname,"w");

  if(!f) {
    perror("image_pairs_dump");
    return ;
  }

  int sizes[3]={
    ip->im0->n,
    ip->im1->n,
    count_ptmatches(ip->imm.ptmatches)
  };
  fwrite(sizes,sizeof(int),3,f);
  
  write_points_file(f,ip->im0->pts,ip->im0->n,2);
  write_points_file(f,ip->im1->pts,ip->im1->n,2);
  
  pointmatch_t *it;
  for(it=ip->imm.ptmatches;it;it=it->next) {
    fwrite(&it->qpt,sizeof(it->qpt),1,f);
    fwrite(&it->dbpt,sizeof(it->dbpt),1,f);
    fwrite(&it->score,sizeof(it->score),1,f);
  } 

  fclose(f);

}



void image_pairs_delete(image_pairs_t *ip) {

  free(ip->pairs);

  pointset_delete(ip->im0);
  pointset_delete(ip->im1);
  fa_free(ip->imm.fa_pointmatch);
  free(ip);

}



static int n_to_na(int n) {
  int na=0;
  while(n>na) na=(na+1)*3/2;
  return na;
}


static int seq_pt(imdesc_t *imdesc,int *map,int p,imdesc_t *imdesc_src) {

  if(map[p]>=0) 
    return map[p];

  int na=n_to_na(imdesc->n);
  
  if(imdesc->n>=na) {
    na=(na+1)*3/2;
    imdesc->pts=realloc(imdesc->pts,sizeof(point_t)*na);
  }
  int mp=imdesc->n++;
  map[p]=mp;
  imdesc->pts[mp]=imdesc_src->pts[p];
  imdesc->pts[mp].dim=0;
  imdesc->pts[mp].desc=NULL;
  return mp;
}


/* query and db images should be vwgeo's */


int image_pairs_add_match_vw(image_pairs_t *ip,
                              imdesc_t *im0,imdesc_t *im1,int hamming_thresh) {

  
  pointmatch_t *pm_end=ip->imm.ptmatches;

  ip->pairs=realloc(ip->pairs,sizeof(two_image_t)*(ip->n_pair+1));
  
  imagematch_align_vw(im0,im1,&ip->imm,hamming_thresh);
  
  ip->pairs[ip->n_pair].pm0=ip->imm.ptmatches;

  int *im0_map=NEWA(int,im0->n+im1->n),*im1_map=im0_map+im0->n;
  
  /* set all to -1 */
  memset(im0_map,0xff,sizeof(im0_map[0])*im0->n);
  memset(im1_map,0xff,sizeof(im1_map[0])*im1->n);
  
  pointmatch_t *it;
  /* handle added matches (they are at the beginning of the list) */
  for(it=ip->imm.ptmatches;it!=pm_end;it=it->next) {
    /* remap pointmatches and 
     * add points to im0 and im1 lists */
    it->qpt=seq_pt(ip->im0,im0_map,it->qpt,im0);
    it->dbpt=seq_pt(ip->im1,im1_map,it->dbpt,im1);    
  }

  free(im0_map);

  return ip->n_pair++;
}

static void ptset_bbox_init(ptset_bbox_t *bb) {
  bb->xmin=bb->ymin=1e30;
  bb->xmax=bb->ymax=-1e30;  
}

static void ptset_bbox_update(ptset_bbox_t *bb,geom_t *geom) {
  if(geom->x<bb->xmin) bb->xmin=geom->x;
  if(geom->x>bb->xmax) bb->xmax=geom->x;
  if(geom->y<bb->ymin) bb->ymin=geom->y;
  if(geom->y>bb->ymax) bb->ymax=geom->y;
}

static int ptset_bbox_isinside(ptset_bbox_t *bb,geom_t *geom) {
  return (bb->xmin<=geom->x && bb->ymin<=geom->y &&
          geom->x<=bb->xmax && geom->y<=bb->ymax);
}


void image_pairs_filter_lowehough(image_pairs_t *ip,
                                  lowehough_parameters_t *params) {

  match_lowehough(ip->im0,ip->im1,&ip->imm,params);  

  /* compute pairwise scores */

  ptset_bbox_init(&ip->bbox0);
  ptset_bbox_init(&ip->bbox1);

  int i;
  for(i=0;i<ip->n_pair;i++) {
    two_image_t *ti=&ip->pairs[i];
    double score=0;
    int n_match=0;
    pointmatch_t *pm,*pm_end=i-1>=0 ? ip->pairs[i-1].pm0 : NULL;
    
    for(pm=ti->pm0;pm!=pm_end;pm=pm->next) {
      if(pm->score>0) {
        score+=pm->score;
        n_match++;

        ptset_bbox_update(&ip->bbox0,&ip->im0->pts[pm->qpt].geom);
        ptset_bbox_update(&ip->bbox1,&ip->im1->pts[pm->dbpt].geom);
      }
    }
    ti->score=score;
    ti->n_match=n_match;
  }
  
}

int ptset_bbox_count_pts(ptset_bbox_t *bb,imdesc_t *imdesc) {
  int i,count=0;
  for(i=0;i<imdesc->n;i++) {
    if(ptset_bbox_isinside(bb,&imdesc->pts[i].geom))
      count++;    
  }
  return count;
}


/* put ugly constants table at end of file */

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


#define NB_BITS_VAL_BINSIGN 64

static int binsign_hamming (binsign_t bs1, binsign_t bs2)
{
  int i, ham = 0;
  binsign_t diff = (bs1 ^ bs2); /* & ((1LL << NB_BITS_VAL_BINSIGN) - 1); */

  for (i = 0 ; i < (NB_BITS_VAL_BINSIGN + 7) / 8 ; i++) {
    ham += uint8_nbones [diff & 255];
    diff >>= 8;
  }

  return ham;
}

