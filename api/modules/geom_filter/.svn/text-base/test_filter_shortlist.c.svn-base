#include <string.h>
#include <assert.h>
#include <stdlib.h>


#include "filter_shortlist.h"
#include "geometry.h"

float *random_floats(int n) {
  int i;
  float *pts=malloc(sizeof(float)*n);
  for(i=0;i<n;i++) pts[i]=drand48();
  return pts;
}

static int *mask;

static void set_mask(void *arg,int n,float *pt2) {
  mask[n]|=2;
}



static double dist2(float *a,float *b,int n) {
  double accu=0;
  while(n--) {
    double d=a[n]-b[n];
    accu+=d*d;
  }
  return accu;
}

void test_match(const char *filelist,const char *pcafile,
                const char *queryfile,double thr,double fthr) {

  imdesc_t **imdescs=NULL;
  int nim=0;
  {
    printf("loading db images\n");
    char buf[1024];
    FILE *f=fopen(filelist,"r");
    assert(f);
    
    while(fgets(buf,1024,f)) {
      *strchr(buf,'\n')=0;
      printf("loading %s\n",buf);
      
      imdescs=realloc(imdescs,sizeof(imdesc_t*)*(nim+1));
      imdesc_t *imdesc=imdescs[nim]=pointset_new();
      
      imdesc->n=read_points(buf,&imdesc->pts,0);    
      
      nim++;
    }
    
    fclose(f);
  }  

  printf("shortlist_new\n");

  shortlist_t *sl=shortlist_new(imdescs,nim);

  printf("loading query\n");
  imdesc_t *query=pointset_new();
  query->n=read_points(queryfile,&query->pts,0);

  printf("shortlist_match_points\n");
  imagematch_t *imms=shortlist_match_points_exact(sl,query,0,thr,NULL,NULL);

  printf("dump point matches dbim qpt dbpt dist^2\n");
  int i;
  for(i=0;i<nim;i++) {
    imagematch_t *imm=&imms[i];
    pointmatch_t *it;
    for(it=imm->ptmatches;it;it=it->next) {
//      printf("%d %d %d %g\n",i,it->qpt,it->dbpt,it->dist);
    }
  }

  lowehough_parameters_t lhp;
  lowehough_parameters_default(&lhp);
  lhp.verbose=10;
 
  shortlist_filter_lowehough(sl,query,imms,&lhp,NULL,NULL);

  for(i=0;i<nim;i++) {
    imagematch_t *imm=&imms[i];

    printf("match with db im %d\n",i);
    printf(
           "  stage1_nmatch=%d\n" 
           "  stage2_nmatch=%d\n" 
           "  stage2_votes=%g\n"
           "  stage3_nmatch=%d\n"
           "  stage3_votes=%g\n"
           "  final_votes=%g\n"
           "  affine[6]=[%g %g %g ; %g %g %g]\n",
           imm->stage1_nmatch, 
           imm->stage2_nmatch, 
           imm->stage2_votes,
           imm->stage3_nmatch,
           imm->stage3_votes,
           imm->final_votes,
           imm->affine[0], imm->affine[1], imm->affine[2],
           imm->affine[3], imm->affine[4], imm->affine[5]);
  }
  
  imagematches_delete(imms,nim);
  
  shortlist_delete(sl);
  
  for(i=0;i<nim;i++) 
    pointset_delete(imdescs[i]);
  free(imdescs);

  pointset_delete(query);
}


void usage(char *progname) {
  fprintf(stderr,"%s \n",progname);
  fprintf(stderr," match <descfilelist> <pcafile> <query> <thr> <fthr>\n");
  exit(1);
}

int main(int argc,char**args) {
 
#define A(cond,...) {                          \
  if(!(cond)) {                                 \
    fprintf(stderr,__VA_ARGS__);                       \
    fprintf(stderr,"\n");                       \
    usage(args[0]);                             \
  }                                             \
}
  
  if(argc<2 || !strcmp(args[1],"-h")  || !strcmp(args[1],"--help")) {
    A(0," ")
  } else if(!strcmp(args[1],"match")) {
    A(argc==7,"wrong # of args");
    test_match(args[2],args[3],args[4],
               atof(args[5]),atof(args[6]));
  } else {
    A(0,"unknown command %s",args[1]);
  } 


  return 0;
}
