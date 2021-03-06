#ifndef _feature_h_
#define _feature_h_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

using namespace std;

#include <iostream>
#include <fstream>
#include <vector>
#include "../ImageContent/imageContent.h"


#define SCALE_FACTOR  1

#define FALSE 0
#define TRUE 1
#define  JLA 2
#define  SHAPE 3 
#define  CC 4 
#define  CF 5
#define  MOM 6
#define  KOEN 7
#define  SPIN 8
#define  PCA 9
#define  SIFT 10
#define  MSIFT 11
#define  GLOH 1
#define  HARRIS 5
#define  HESSIAN 1
#define  HARHES 2
#define  DENSE 100
#define  EDGE 3
#define  SEDGE 4
#define  AFFINE 16
#define  SIFTNONORM 17

extern DARY *patch_mask;
extern float PATCH_SUM;
#define  PATCH_SIZE  21

void initPatchMask(int size);
class  FeatureDescriptor;
typedef FeatureDescriptor FD;

class Occurence{
 public:
  float scale,angle,fangle,fscale,dist,weight,x,y,area,el,ea; 
  int obj;
  FD *feat;
  Occurence(int x_in, int y_in, int scale_in,float angle_in, float weight_in, int obj_in){
    x=x_in;
    y=y_in;
    scale=scale_in;
    angle=angle_in; //absolut angle
   weight=weight_in;
    obj=obj_in;
  }

  Occurence(Occurence *occ){
    x=occ->x;
    y=occ->y;
    scale=occ->scale;
    fscale=occ->fscale;
    angle=occ->angle; //absolut angle
    fangle=occ->fangle; //absolut angle
    weight=occ->weight;
    obj=occ->obj;
  }

  Occurence(float x_in, float y_in, float scale_in,float angle_in,float fangle_in, float area_in, float dist_in, float weight_in, int obj_in){
    x=x_in;
    y=y_in;
    scale=scale_in;
    angle=angle_in; //absolut angle
    fangle=fangle_in; // angle related to the object center
    dist=dist_in;
    weight=weight_in;
    obj=obj_in;
    area=area_in;
    fscale=1;
  }
  ~Occurence(){}
  void Cout(){
    cout << x << " "<< y << " s " << fscale << " " << weight<< endl;
  }

  void out2(ofstream &out, int width_2, int height_2, float image_scale){
    //cout << x << " x "<< y << " "  << width_2 << " "<< fscale << " " << image_scale<< endl;
    out  << obj << " " << (x-(width_2*fscale))/image_scale << " "<< (y-height_2*fscale)/image_scale << " " << (x+(width_2*fscale))/image_scale << " "<< (y+height_2*fscale)/image_scale << " " << weight << endl;
  }
};




class  FeatureDescriptor{

 protected:
 public:
    float *vec;  
    char *imagename;
    DARY *cluster_dist;
    float mean_dist;
    uint nbf;
    int size, tree_lev, obj;
    float d_scale, weight, var, sim, radius, area;    
    float x,y,l1,l2,lap, mi11,mi12,mi21,mi22,c_scale, int_sig, der_sig,featureness;
    float angle,eangle; //dominant orientation, ellipse angle
    int int_lev, der_lev, extr, type;
    vector<FeatureDescriptor *> features;
    vector<Occurence *> occurences;
    vector<float> neg_match;
    vector<float> obj_weight;

 public:
     FeatureDescriptor(){init();} 
     void init();
     
     FeatureDescriptor(float xin, float yin, float scale_in, float featureness_in);

     void copy(FeatureDescriptor* ds);

    /*****READ WRITE***/
    void readCommon( ifstream &input, int size_in);
    void writeCommon( ofstream &output);
    void write( ofstream &output);
    void read( ifstream &input, int size_in);	
    void writeBin(FILE *fid, float *buf);
    void readBin(FILE *fid, int size, float *buf);

 
    void allocVec(int);
    inline float * getVec(void){return vec;}
    inline float getV(int i){ if(i<size)return (vec[i]);else return 0;}    
    inline void setV(int i, float val){if(i<size && i>=0)vec[i]=val;}

    inline int const getSize() const {return size;} 
    inline void setSize(int size_in){size=size_in;} 
  
    inline float const getSim() const {return sim;} 
    inline void setSim(float sim_in){sim=sim_in;}   
 
    inline float const getRadius() const {return radius;} 
    inline void setRadius(float radius_in){radius=radius_in;}    

    inline int const getTreeLevel() const {return tree_lev;} 
    inline void setTreeLevel(int tree_lev_in){tree_lev=tree_lev_in;}    

    inline float const getWeight() const {return weight;} 
    inline void setWeight(float weight_in){weight=weight_in;}    

    inline float const getArea() const {return area;} 
    inline void setArea(float area_in){area=area_in;}    

    inline float const getMeanDist() const {return mean_dist;} 
    inline void setMeanDist(float d_in){mean_dist=d_in;}    

    inline uint const getNbf() const {return nbf;} 
    inline void setNbf(uint nbf_in){nbf=nbf_in;}    

    inline int   const  getType() const { return type;}
    inline void     setType(int type_in)  {type=type_in;}


    inline int   const  getObj() const { return obj;}
    inline void     setObj(int obj_in)  {obj=obj_in;}
    
  
    /****CORNERNESS***/
    inline float const getFeatureness(void) const { return featureness;}
    inline void setFeatureness(float featureness_in) {featureness=featureness_in;}
    
    inline void setMi(float m11,float m12,float m21,float m22) 
    {mi11=m11;mi12=m12;mi21=m21;mi22=m22;}
    inline void setMi(float m11,float m12,float m21,float m22,float e1,float e2, float ea) 
      {mi11=m11;mi12=m12;mi21=m21;mi22=m22;l1=e1;l2=e2;eangle=ea;}
    inline float   const  getMi11() const { return mi11;}
    inline float   const  getMi12() const { return mi12;}
    inline float   const  getMi21() const { return mi21;}
    inline float   const  getMi22() const { return mi22;}
    inline float   const  getL1() const { return l1;}
    inline float   const  getL2() const { return l2;}
    inline void   setL2(float e2) {  l2=e2;}
    
    /****LOCALISATION***/
    inline float   const  getX() const { return x;}
    inline float   const  getY() const { return y;}
    inline void     setX_Y(float x_in,float y_in)  {x=x_in;y=y_in;}


    /****LAPLACIAN VALUE***/
    inline float   const  getLap() const { return lap;}
    inline void setLap(float lap_in) {lap=lap_in;}
    inline void setExtr(int ex) {extr=ex;}
    inline void setMax() {extr=1;}
    inline void setMin() {extr=-1;}
    inline int   const  getExtr() const { return extr;}
    

    /****ANGLE*****/
    inline float   const  getAngle() const { return angle;}
    inline void    setAngle(const float angle_in){ angle=angle_in;}

    inline float   const  getEAngle() const { return eangle;}
    inline void    setEAngle(const float eangle_in){ eangle=eangle_in;}

    inline float   const  getVar() const { return var;}
    inline void    setVar(const float var_in){ var=var_in;}

    /*****SCALE******/
    inline float const  getScale() const { return c_scale;}
    inline void setScale(const float scale_in){c_scale=scale_in;}
    inline void setIntSig(const float sig_in){int_sig=sig_in;}
    inline float const getIntSig() const { return int_sig;}
    inline void setDerSig(const float sig_in){der_sig=sig_in;}
    inline float const  getDerSig() const { return der_sig;}
    inline int const  getIntLev() const { return int_lev;}
    inline int const  getDerLev() const { return der_lev;}
    inline void setIntLev(int lev_in) {int_lev=lev_in;}
    inline void setDerLev(int lev_in) {der_lev=lev_in;}
    
    
    inline const char*   getImageName() const { return imagename;}
    void setImageName(const char* name);
    
    void addOccurence(float x, float y, float scale,float angle,float fangle, float area, float dist, float sim, int obj){
      occurences.push_back(new Occurence(x,y,scale,angle,fangle,area,dist,sim,obj));
    }  
    void findNearestNeighbor(vector<FeatureDescriptor*> features, int dim, int &ind, float &sim);

    void changeBase(float* mat);
    void pca(int dim,  float *avg, float *base);
    void Cout(int dim=3);


    ~FeatureDescriptor();
   
};

typedef vector <FeatureDescriptor*> FeatureDescList;
typedef vector <FeatureDescList> FeatureDescSequence;

inline float square(float a){return a*a;}

int normalize(DARY * img,int x, int y, float radius);

void deleteDescriptors(vector<FeatureDescriptor *> &desc);

DARY* normalizeAffine(DARY *image, float x, float y, float c_scale, 
		      float angle, float mi11, float mi12,float , float mi22, float &scal,
		      DARY *patch, float DESC_CALE);
void computeAffineDescriptor( DARY *imgbs, DARY *patch, float scal, int DESC,
			FeatureDescriptor *cor, vector<FeatureDescriptor *> &desc);
void computeAffineDescriptor( DARY *image, DARY *patch, int DESC, 
			FeatureDescriptor *cor, vector<FeatureDescriptor *> &desc, float DESC_CALE);
void computeAffineDescriptors(DARY *image,  vector<FeatureDescriptor *> &desc, int DESC, float DESC_CALE);

void computeDescriptor( DARY *gpatch, DARY *opatch, int DESC,
			FeatureDescriptor *cor, vector<FeatureDescriptor *> &desc);

void loadBinFeatures(const char* points_out, vector<FeatureDescriptor*> &cor);
void saveBinFeatures(vector<FeatureDescriptor*> cor, const char* points_out);


void writeFeatures(vector<FeatureDescriptor*> cor, const char* points_out, int format=0);
void writeFeaturesHerve(vector<FeatureDescriptor*> cor, const char * points_out);
void writeFeaturesSameh(const vector<FeatureDescriptor*> & cor, const char * points_out);
void writeCoordinates(vector<FeatureDescriptor*> cor, const char * coord_file);

void loadFeatures( const char* points1, vector<FeatureDescriptor*> &cor1, int format=0);
void loadFredFeatures( const char* points1, vector<FeatureDescriptor*> &cor1);





int computeSift(DARY *dx, DARY *dy , DARY *grad, DARY *gori ,  FeatureDescriptor *ds, int do_normalization);
void computeESift(DARY *dx, DARY *dy ,DARY *grad, DARY *gori , FeatureDescriptor *ds);
void computeShape(DARY *dx, DARY *dy ,DARY *grad, DARY *gori, FeatureDescriptor *ds);


void computeJLA(DARY *patch, FeatureDescriptor *ds);
void computeSift(DARY *patch, FeatureDescriptor *ds, int do_normalization);
void computeESift(DARY *patch, FeatureDescriptor *ds);
void computeMoments(DARY *patch, FeatureDescriptor *ds);
void computeKoen(DARY *patch, FeatureDescriptor *ds);
void computeCF(DARY *patch, FeatureDescriptor *ds);
void computeShape(DARY *patch, FeatureDescriptor *ds);
void computeSpin(DARY *patch, FeatureDescriptor *ds);
void computePcaSift(DARY *patch, FeatureDescriptor *ds);
void computeCC(DARY *patch, FeatureDescriptor *ds);
 


void displayFeatures(DARY *image, vector<FeatureDescriptor*> cor, char *filename, float color, int circle=0);
void displayFeatures(const char *filein, vector<FeatureDescriptor*> cor, char *filename, float color, int circle=0);


void cannyEdges(DARY *img, DARY *edge,  float scale, float lower_threshold, float higher_threshold);
void cannyEdges(DARY *dx, DARY *dy, DARY *edge,  float lower_threshold, float higher_threshold);
void cannyEdges(DARY *img, DARY *edge,  float lower_threshold, float higher_threshold);
void cannyEdges(DARY *dx, DARY *dy, DARY *grad, DARY *edge,  float lower_threshold, float higher_threshold);

void fastSEdge(DARY* img, vector<FeatureDescriptor*>&features, float threshold,float step, int DESC, int noangle);
void fastEdge(DARY* img, vector<FeatureDescriptor*>&features, float threshold,float step, int DESC, int noangle);
void fastHarrisHessian(DARY *img, vector<FeatureDescriptor*>&corners, float threshold, float step, int DESC, int aff, int noangle);
void fastHessian(DARY *img, vector<FeatureDescriptor*>&features, float threshold, float step, int DESC, int aff, int noangle, int max_desc);
void fastHarris(DARY *img, vector<FeatureDescriptor*>&features, float threshold, float step, int DESC, int aff, int noangle);
void pca(vector<FD*> &features,float *mvec, float *base, int newdim);
void pcaBase(vector<FD*> &features,float *mvec, float *base, int newdim);
void savePCAVectors(const char *filename, float *vmean, uint size1 ,float *vbase, uint size2);
int loadPCAVectors(const char *filename, float *&vmean, uint &mean_dim, float *&vbase, uint &base_dim);

void extractFeatures(DARY *image, vector<FD*> &features, int detector, int descriptor, int affine, int noangle, float thres, int max_desc);
void detectFeatures(const char *filein, vector<FD *> &features, int mask, vector<pair<int,int> > &dims, int detector, int DECSRIPTOR, float image_scale, int image_shift, float det_threshold, int max_desc);
void detectFeatures(const char *name, DARY *img, DARY *mask, vector<FD *> &features,  vector<pair<int,int> > &dims, int detector, int DECSRIPTOR, float image_scale, int image_shift, float det_threshold, int max_desc);

/* low level for 1 feature */
int  compDesc(DARY *dx, DARY *dy, DARY *grad, DARY *ori, FD *ds, int DESC);

void fastDense(DARY *img, vector<FeatureDescriptor*>&features, 
               int DESC, int xstep, int ystep, int max_pyrlevels );

#endif
