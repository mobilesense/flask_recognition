#include <stdio.h>
#include <stdlib.h>
#include <assert.h>


#include "geometry.h"
#include "filter_shortlist.h"

#define NEWA(type,n) (type*)malloc(sizeof(type)*(n))
#define NEWAC(type,n) (type*)calloc(sizeof(type),(n))
#define NEW(type) NEWA(type,1)

#define MAX(a,b) ((a)>(b) ? (a) : (b))
#define MIN(a,b) ((a)<(b) ? (a) : (b))

/* for compatibility with C++ version */
#define DESC_BYTE_TO_FLOAT(d) (((d)+0.5)/512.0)


/* #define DESC_BYTE_TO_FLOAT(d) ((float)(d)) */ 




/*-----------------------------------------------------------------
 * Geometry related functions
 *-----------------------------------------------------------------*/


#define integer int

/* LAPACK routine for least squares solution of linear system */
extern int dgels_(char *trans, integer *m, integer *n, integer *
                  nrhs, double *a, integer *lda, double *b, integer *ldb,
                  double *work, integer *lwork, integer *info);


/* return <0 for bad estimation */
int estimate_affine_transform(point_t *qpts,point_t *dbpts,
                              pointmatch_t **pms,int npm,
                              double aff[6]) {
  if(npm<3) return -1;
  
  int maxpm=2000;
  int npmused=MIN(npm,maxpm);
  int ret=0;

  double *mat=NEWA(double,npmused*12);
  double *vec=NEWA(double,npmused*2);

  { /* fill the linear system */

    int i;
    double *m=mat, *v=vec;

    for (i=0;i<npmused;i++) { 
      int j= npm<=maxpm ? i : i*npm/maxpm;
      pointmatch_t *pm=pms[j];
      geom_t *qpt=&qpts[pm->qpt].geom;
      geom_t *dbpt=&dbpts[pm->dbpt].geom;

      m[0] = m[8] = qpt->x;
      m[1] = m[9] = qpt->y;
      m[2] = m[3] = m[5] = m[6] = m[7] = m[10] = 0.0; 
      m[4] = m[11] = 1.0; 
      
      v[0] = dbpt->x;
      v[1] = dbpt->y;

      v+=2;
      m+=12;
    } 

  }
  
  { /* solve system */ 
    integer info;
    integer m=6, nrhs=1, lda=6, lwork=-1;
    integer n=2*npmused,ldb=n;
    double work_sz;

    dgels_("Transposed", &m, &n, &nrhs, mat, &lda, 
           vec, &ldb, &work_sz, &lwork, &info); 
    
    double *work;    
    lwork=(int)work_sz;
    work=NEWA(double,lwork);
    
    dgels_("Transposed", &m, &n, &nrhs, mat, &lda, 
           vec, &ldb, work, &lwork, &info); 
    
    free(work);
  
    assert(info>=0); /* there is always a result for coherent input */

    /* Not documented in LAPACK: info>0 for (some?) rank deficient matrices */

    ret=-info;
  }
  
  aff[0]=vec[0]; aff[1]=vec[1]; aff[2]=vec[4];
  aff[3]=vec[2]; aff[4]=vec[3]; aff[5]=vec[5];

  free(mat);
  free(vec);

  return ret;
}


static void transform_affine(double aff[6],double x,double y,
                             double *x_out,double *y_out) {
  *x_out=aff[0]*x+aff[1]*y+aff[2];
  *y_out=aff[3]*x+aff[4]*y+aff[5];
}



int find_agreeing_matches(point_t *qpts,point_t *dbpts,
                          double aff[6],
                          pointmatch_t *pm,
                          double thresh,
                          pointmatch_t **agreeing_matches,int* ids) {
  double t2=thresh*thresh;
  int na=0; int i=0;

  while(pm) {
    geom_t *qpt=&qpts[pm->qpt].geom;
    geom_t *dbpt=&dbpts[pm->dbpt].geom;

    double mqx,mqy;
    transform_affine(aff,qpt->x,qpt->y,&mqx,&mqy);

    double dx=mqx-dbpt->x,dy=mqy-dbpt->y;

    if(dx*dx+dy*dy<t2){
      agreeing_matches[na++]=pm;
      ids[i]=1;      
    }  
    else{
      ids[i]=0;
    }

    pm=pm->next;
    i++;
    
  }

  return na;
}

int count_distinct_points(point_t *pts,
                          pointmatch_t **pms,int npm,
                          int use_query,int stopat,
                          double radius) {
  typedef struct {int q; double x; double y;} qxy;
  qxy *qp=NEWA(qxy,stopat);  
  int count=0;
  double r2=radius;
  int i;

  for(i=0;i<npm && count<stopat;i++) {
    pointmatch_t *pm=pms[i];
    int ptno=use_query ? pm->qpt : pm->dbpt;    
    double x=pts[ptno].geom.x,y=pts[ptno].geom.y;

    int j;
    for(j=0;j<count;j++) {
      double dx=qp[j].x-x;
      double dy=qp[j].y-y;

      if(ptno==qp[j].q || dx*dx+dy*dy<r2) 
        goto dontkeep;
    }   
    
    qp[count].x=x;
    qp[count].y=y;
    qp[count].q=ptno;
    count++;

  dontkeep: ;
  }

  free(qp);

  return count;
}



#if 0

struct {
  double min,max;
} minmax_t;

static void init_minmax(minmax_t *mm) {
  mm->min=1.0/0.0;
  mm->max=mm->min;
}





/* minimum range of the points when projected on a line, over line orientations
 *  TODO: translate from ~/test/narrow/narrow.py
 */ 

static double narrowest_diameter(point_t *pts,
                                 pointmatch_t **pms,int npm,
                                 int use_query) {
  
  for(i=0;i<npm;i++) {
    pointmatch_t *pm=pms[i];
    int ptno=use_query ? pm->qpt : pm->dbpt;    
    double x=pts[ptno].x,y=pts[ptno].y;

    /* TODO in O(n*log(n)) */

  }  

}

#endif
