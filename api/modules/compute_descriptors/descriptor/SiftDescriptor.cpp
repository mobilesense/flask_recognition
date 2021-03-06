#include "feature.h"
#include "../gauss_iir/gauss_iir.h"

#include "sift_base.h"
#include "sift_pca_base.h"
#include "esift_base.h"
#include "esift_pca_base.h"
#include "sc_base.h"
#include "sc_pca_base.h"
#include "pca_base.h"
#include "pca_eig.h"

const int LocSize=4;
const int OriSize=8;
const int SiftSize=LocSize*LocSize*OriSize;
const int MagFactor=2;
const float MaxIndexVal = 0.2;  /* Good value is 0.2 */

const int ELocSize=5; 
const int EOriSize=12;
const int ESiftSize=ELocSize*ELocSize*EOriSize;

const int SrSize=3;
const int ScSize=4;
const int SLocSize=4;
const int SOriSize=8;
//const int ShapeSize=SLocSize*SLocSize*SOriSize;
const int ShapeSize=SrSize*ScSize*SOriSize;

const int GPLEN=3042;
const int PCALEN=36;

void MakeKeypointSample(DARY *grad, DARY *ori, float *vec, int size, int oriSize, int locSize);
void AddSample(float *index, DARY *grad, DARY *orim, float angle, int r, int c, float rpos, 
	       float cpos,float rx, float cx, int oriSize, int locSize);
void PlaceInIndex(float *index,
		  float mag, float ori, float rx, float cx, int oriSize, int locSize);
void KeySample(float *index, DARY *grad, DARY *ori, float angle, int oriSize, int locSize);
int NormalizeVect(float *vec, int len);
void NormalizeEVect(float *vec, int len, float norm);



/* David Lowe's code*/

/* Increment the appropriate locations in the index to incorporate
   this image sample.  The location of the sample in the index is (rx,cx). 
*/

void PlaceInLogPolIndex(float *index,
			float mag, float ori, float rx, float cx, int oriSize, int rSize, int cSize)
{
   int r, c, ort, ri, ci, oi, rindex, cindex, oindex, rcindex;
   float oval, rfrac, cfrac, ofrac, rweight, cweight, oweight;
   
   oval = oriSize * ori / (M_2PI);
   
   ri = (int)((rx >= 0.0) ? rx : rx - 1.0);  /* Round down to next integer. */
   ci = (int)((cx >= 0.0) ? cx : cx - 1.0);
   oi = (int)((oval >= 0.0) ? oval : oval - 1.0);
   rfrac = rx - ri;         /* Fractional part of location. */
   cfrac = cx - ci;
   ofrac = oval - oi; 
   /*   assert(ri >= -1  &&  ri < IndexSize  &&  oi >= 0  &&  oi <= OriSize  && rfrac >= 0.0  &&  rfrac <= 1.0);*/
   //cout << ri << " " << ci << " " << oi << endl;
   /* Put appropriate fraction in each of 8 buckets around this point
      in the (row,col,ori) dimensions.  This loop is written for
      efficiency, as it is the inner loop of key sampling. */
   for (r = 0; r < 2; r++) {
      rindex = ri + r;
      if (rindex >=0 && rindex < rSize) {
         rweight = mag * ((r == 0) ? 1.0 - rfrac : rfrac);
         
         for (c = 0; c < 2; c++) {
            cindex = ci + c;
	    if(cindex >= cSize)
	      cindex=0;	    
            if (cindex >=0 && cindex < cSize) {
               cweight = rweight * ((c == 0) ? 1.0 - cfrac : cfrac);
               rcindex=(rindex*cSize+cindex)<<3;//remember when you change the orientation number
               for (ort = 0; ort < 2; ort++) {
                  oindex = oi + ort;
                  if (oindex >= oriSize)  /* Orientation wraps around at PI. */
                     oindex = 0;
                  oweight = cweight * ((ort == 0) ? 1.0 - ofrac : ofrac);
                  //cout << rcindex+oindex<< endl;
		  index[rcindex+oindex]+=oweight;
               }
            }  
         }
      }
   } 
}

/* Given a sample from the image gradient, place it in the index array.
*/
void AddLogPolSample(float *index,
		     DARY *grad, DARY *orim, float angle, int r, int c, float rpos, float cpos,
		     float rx, float cx, int oriSize, int rSize, int cSize)
{
    float mag, ori;
    
    /* Clip at image boundaries. */
    if (r < 0  ||  r >= (int)grad->y()  ||  c < 0  ||  c >= (int)grad->x())
       return;
    
    mag = patch_mask->fel[r][c] * grad->fel[r][c];
    //mag = grad->fel[r][c];
    /* Subtract keypoint orientation to give ori relative to keypoint. */
    ori = orim->fel[r][c]-angle;
    
    /* Put orientation in range [0, 2*PI].  If sign of gradient is to
       be ignored, then put in range [0, PI]. */
    
    while (ori > M_2PI)
      ori -= M_2PI;
    while (ori < 0.0)
      ori += M_2PI;     
    PlaceInLogPolIndex(index, mag, ori, rx, cx, oriSize, rSize, cSize);
} 
void KeyLogPolSample(float *index, DARY *grad, DARY *ori, float angle, int oriSize, int rSize, int cSize)
{
   int i, j, iradius;
   float rspacing,cspacing, rpos, cpos, lrpos, lcpos, rx, cx;
   
   /* Radius of index sample region must extend to diagonal feature of
      index patch plus half sample for interpolation. */
   
   iradius = grad->x()>>1;
   rspacing = (rSize) / (iradius);
   cspacing = (cSize) / (M_2PI);
   float sine = sin(-angle);
   float cosine = cos(-angle);

   //   printf("spacing %f, scale %f, radius %d\n",spacing,scale, iradius);
   // Examine all points from the gradient image that could lie within the index square. 
   for (i = -iradius; i <= iradius; i++)
     for (j = -iradius; j <= iradius; j++) {
       

       rpos = (cosine*i + sine*j);// coreect with dominant angle
       cpos = (-sine*i + cosine*j);
       lcpos=(M_PI+atan2(rpos,cpos))*cspacing;
       lrpos=log(1+sqrt((float)i*i+j*j)/iradius)*rSize;
       
       // Compute location of sample in terms of real-valued index array
       //  coordinates.  Subtract 0.5 so that rx of 1.0 means to put full
       // weight on index[1] (e.g., when rpos is 0 and IndexSize is 3. 
       rx = lrpos;// + (rSize - 1) / 2.0;
       cx = lcpos;// + (cSize - 1) / 2.0;
       //cout << rx << " "<< cx << endl;
       
       //cout <<"rx " << rx << " " << cx << endl;
       // Test whether this sample falls within boundary of index patch. 
       if (rx > -1.0 && rx < (float) rSize  &&
	   cx > -1.0 && cx < (float) cSize)
         //cout << "in" << cpos << " " << rpos << endl;
	 AddLogPolSample(index, grad, ori, angle, iradius + i, iradius + j, lrpos, lcpos,
		   rx, cx, oriSize, rSize, cSize);
   } 
}


void PlaceInIndex(float *index,
		  float mag, float ori, float rx, float cx, int oriSize, int locSize)
{
   int r, c, ort, ri, ci, oi, rindex, cindex, oindex, rcindex;
   float oval, rfrac, cfrac, ofrac, rweight, cweight, oweight;
   
   oval = oriSize * ori / (M_2PI);
   
   ri = (int)((rx >= 0.0) ? rx : rx - 1.0);  /* Round down to next integer. */
   ci = (int)((cx >= 0.0) ? cx : cx - 1.0);
   oi = (int)((oval >= 0.0) ? oval : oval - 1.0);
   rfrac = rx - ri;         /* Fractional part of location. */
   cfrac = cx - ci;
   ofrac = oval - oi; 
   /*   assert(ri >= -1  &&  ri < IndexSize  &&  oi >= 0  &&  oi <= OriSize  && rfrac >= 0.0  &&  rfrac <= 1.0);*/
   
   /* Put appropriate fraction in each of 8 buckets around this point
      in the (row,col,ori) dimensions.  This loop is written for
      efficiency, as it is the inner loop of key sampling. */
   for (r = 0; r < 2; r++) {
      rindex = ri + r;
      if (rindex >=0 && rindex < locSize) {
         rweight = mag * ((r == 0) ? 1.0 - rfrac : rfrac);
         
         for (c = 0; c < 2; c++) {
            cindex = ci + c;
            if (cindex >=0 && cindex < locSize) {
               cweight = rweight * ((c == 0) ? 1.0 - cfrac : cfrac);
               rcindex=(rindex*locSize+cindex)<<3;//remember when you change the orientation number
               for (ort = 0; ort < 2; ort++) {
                  oindex = oi + ort;
                  if (oindex >= oriSize)  /* Orientation wraps around at PI. */
                     oindex = 0;
                  oweight = cweight * ((ort == 0) ? 1.0 - ofrac : ofrac);
                  
		  index[rcindex+oindex]+=oweight;
               }
            }  
         }
      }
   } 
}
 
  
/* Given a sample from the image gradient, place it in the index array.
*/
void AddSample(float *index,
	       DARY *grad, DARY *orim, float angle, int r, int c, float rpos, float cpos,
	       float rx, float cx, int oriSize, int locSize)
{
    float mag, ori;
    
    /* Clip at image boundaries. */
    if (r < 0  ||  r >= (int)grad->y()  ||  c < 0  ||  c >= (int)grad->x())
       return;
    
    /* Compute Gaussian weight for sample, as function of radial distance
       from center.  Sigma is relative to half-width of index. */
    //sigma =  0.5*(IndexSize+1)*(IndexSize+1);
    //sigma = (IndexSize+1)*(IndexSize+1);
    //weight = 0.6*exp(- (rpos * rpos + cpos * cpos) / (sigma) );
    //cout << "rpos "<< rpos << " cpos "<< cpos << " weight " << weight<< endl;
    //mag = weight * grad->fel[r][c];
    //mag = grad->fel[r][c];
    mag = patch_mask->fel[r][c] * grad->fel[r][c];
    /* Subtract keypoint orientation to give ori relative to keypoint. */
    ori = orim->fel[r][c]-angle;
    
    /* Put orientation in range [0, 2*PI].  If sign of gradient is to
       be ignored, then put in range [0, PI]. */
    
    while (ori > M_2PI)
      ori -= M_2PI;
    while (ori < 0.0)
      ori += M_2PI;     
    PlaceInIndex(index, mag, ori, rx, cx, oriSize, locSize);
} 

/* Add features to vec obtained from sampling the grad and ori images
   for a particular scale.  Location of key is (scale,row,col) with respect
   to images at this scale.  We examine each pixel within a circular
   region containing the keypoint, and distribute the gradient for that
   pixel into the appropriate bins of the index array.
*/

void KeySample(float *index, DARY *grad, DARY *ori, float angle, int oriSize, int locSize)
{
   int i, j, iradius;
   float spacing, rpos, cpos, rx, cx;
   
   /* Radius of index sample region must extend to diagonal feature of
      index patch plus half sample for interpolation. */
   
   iradius = grad->x()>>1;
   spacing = (locSize + 1) / (2.0*iradius);
   float sine = sin(-angle);
   float cosine = cos(-angle);

   //   printf("spacing %f, scale %f, radius %d\n",spacing,scale, iradius);
   // Examine all points from the gradient image that could lie within the index square. 
   for (i = -iradius; i <= iradius; i++)
     for (j = -iradius; j <= iradius; j++) {
       
       rpos = (cosine*i + sine*j) * spacing;
       cpos = (-sine*i + cosine*j) * spacing;
        
       // Compute location of sample in terms of real-valued index array
       //  coordinates.  Subtract 0.5 so that rx of 1.0 means to put full
       // weight on index[1] (e.g., when rpos is 0 and IndexSize is 3. 
       rx = rpos + (locSize - 1) / 2.0;
       cx = cpos + (locSize - 1) / 2.0;
       
       //cout <<"rx " << rx << " " << cx << endl;
       // Test whether this sample falls within boundary of index patch. 
       if (rx > -1.0 && rx < (float) locSize  &&
	   cx > -1.0 && cx < (float) locSize)
         //cout << "in" << cpos << " " << rpos << endl;
	 AddSample(index, grad, ori, angle, iradius + i, iradius + j, rpos, cpos,
		   rx, cx, oriSize, locSize);
   } 
}

 
/* Normalize length of vec to 1.0.
*/
int NormalizeVect(float *vec, int len)
{
   int i;
   float val, fac, sqlen = 0.0;

   for (i = 0; i < len; i++) {
     val = vec[i];
     sqlen += val * val;
   }
   if(sqlen<=0) return 0;
   fac = 1.0 / sqrt(sqlen);
   for (i = 0; i < len; i++)
     vec[i] *= fac;
   return 1;
}

void NormalizeEVect(float *vec, int len, float norm)
{
   float fac, sqlen = 0.0;

   for (int i = 0; i < len; i++) {
     sqlen += fabs(vec[i]);
   }
   //   fac = norm / sqrt(sqlen);
   fac=sqlen/len;
   for (int i = 0; i < len; i++)
     vec[i] = vec[i]/fac *100.0;
}
   


void computeSift(DARY *img, FeatureDescriptor *ds, int do_normalization){
  
  DARY * grad = new DARY(PATCH_SIZE,PATCH_SIZE);
  DARY * gori = new DARY(PATCH_SIZE,PATCH_SIZE);    
  
  gradAngle(img,grad,gori);
  ds->allocVec(SiftSize);
  float *vec = ds->getVec();
  KeySample(vec, grad, gori, 0.0, OriSize, LocSize);

  float norm_factor; 

  if(do_normalization) {

    NormalizeVect(vec, SiftSize); 
    int changed = FALSE;
    for (int i = 0; i < SiftSize; i++)
      if (vec[i] > MaxIndexVal) { 
        vec[i] = MaxIndexVal;
        changed = TRUE;
      }
    if (changed)
      NormalizeVect(vec, SiftSize); 

    /* Convert float vector to integer. Assume largest value in normalized
       vector is likely to be less than 0.5. */
    
    norm_factor = 512.0; 

  } else {
    /* no normalization, but then we need some magic factor to convert to byte */
    norm_factor = 1.5 * 21.0 * 21.0 / (PATCH_SIZE * PATCH_SIZE);
  } 

  for (int i = 0; i < SiftSize; i++) {
    int intval = (int) (norm_factor * vec[i]);
    vec[i] = (255 < intval) ? 255 : intval;
  }
  
  delete grad;delete gori;
  //int sift_pca_size=128;
   //ds->pca(sift_pca_size, sift_pca_avg, sift_pca_base);	
}


int computeSift(DARY *dx, DARY *dy ,DARY *grad, DARY *gori , FeatureDescriptor *ds, int do_normalization){
  ds->allocVec(SiftSize);
  float *vec = ds->getVec();
  float angle=ds->getAngle();
  if(angle<-M_PI || angle>M_PI)angle=0;

  KeySample(vec, grad, gori, angle, OriSize, LocSize);

  float norm_factor; 

  if(do_normalization) {

    if(!NormalizeVect(vec, SiftSize)) return 0;
 
    int changed = FALSE;
    for (int i = 0; i < SiftSize; i++)
      if (vec[i] > MaxIndexVal) { 
        vec[i] = MaxIndexVal;
        changed = TRUE;
      }
    if (changed)
      NormalizeVect(vec, SiftSize); 

    /* Convert float vector to integer. Assume largest value in normalized
       vector is likely to be less than 0.5. */
   
    norm_factor = 512.0; 

  } else {
    /* no normalization, but then we need some magic factor to convert to byte */
    norm_factor = 1.5 * 21.0 * 21.0 / (PATCH_SIZE * PATCH_SIZE);
  } 
  
  //  printf("nf=%g\n",norm_factor);

  for (int i = 0; i < SiftSize; i++) {
    int intval = (int) (norm_factor * vec[i]);
    vec[i] = (255 < intval) ? 255 : intval;
  }
  
 //int sift_pca_size=128;
   //ds->pca(sift_pca_size, sift_pca_avg, sift_pca_base);	
  return 1;
}
 
void computeESift(DARY *img, FeatureDescriptor *ds){
  
  DARY * grad = new DARY(PATCH_SIZE,PATCH_SIZE);
  DARY * gori = new DARY(PATCH_SIZE,PATCH_SIZE);    
  
  gradAngle(img,grad,gori);
  ds->allocVec(ESiftSize);//EIndexSize*EIndexSize*OriSize;
  float *vec = ds->getVec();

  KeySample(vec, grad, gori, 0.0, EOriSize, ELocSize);
  //grad->writePNG("grad.png");cout << "patch 1" << endl;getchar();
  delete grad;delete gori;
  NormalizeEVect(vec, ESiftSize, 1.0); 
  
  int esift_pca_size=128; 
  ds->pca(esift_pca_size, esift_pca_avg, esift_pca_base);

}


void computeESift(DARY *dx, DARY *dy ,DARY *grad, DARY *gori , FeatureDescriptor *ds){
  ds->allocVec(ESiftSize);//EIndexSize*EIndexSize*OriSize;
  float *vec = ds->getVec();
  float angle=ds->getAngle();
  if(angle<-M_PI || angle>M_PI)angle=0;
  KeySample(vec, grad, gori, angle, EOriSize, ELocSize);
  //grad->writePNG("grad.png");cout << "patch 2" << endl;getchar();
  NormalizeEVect(vec, ESiftSize, 1.0); 
  int esift_pca_size=128; 
  ds->pca(esift_pca_size, esift_pca_avg, esift_pca_base);	
}


/******************SHAPE CONTEXT******************************************/


void computeShape(DARY *img, FeatureDescriptor *ds){

  DARY * grad = new DARY(PATCH_SIZE,PATCH_SIZE);
  DARY * ori = new DARY(PATCH_SIZE,PATCH_SIZE);    
  DARY *dx = new DARY(img->y(),img->x());
  DARY *dy = new DARY(img->y(),img->x());
  DARY *edge = new DARY(img->y(),img->x());
  //float angle=ds->getAngle();
  //canny(img,grad,ori,5,15);
  dX2(img,dx);
  dY2(img,dy);
  for(uint j=0;j<grad->y();j++){
    for(uint i=0;i<grad->x();i++){
      grad->fel[j][i]=sqrt(dx->fel[j][i]*dx->fel[j][i]+dy->fel[j][i]*dy->fel[j][i]); 
      ori->fel[j][i]=atan2(dy->fel[j][i],dx->fel[j][i]);
    }
  } 
  cannyEdges(dx, dy, grad, edge, 5, 15);
  //grad->write("edge.pgm");cout << "OK"<< endl;getchar();
  delete dx; delete dy;delete grad;


  ds->allocVec(ShapeSize);
  float *vec = ds->getVec();
  //KeyLogPolSample(vec, edge, ori,  angle, SOriSize, SrSize, ScSize);
  //getchar();
   KeySample(vec, edge, ori,0.0, SOriSize, SLocSize);
  delete edge;delete ori;
  NormalizeVect(vec, ShapeSize); 
  int intval, changed = FALSE;
  for (int i = 0; i < ShapeSize; i++)
    if (vec[i] > MaxIndexVal) { 
      vec[i] = MaxIndexVal;
       changed = TRUE;
     }
  if (changed)
    NormalizeVect(vec, ShapeSize); 
  /* Convert float vector to integer. 
     Assume largest value in normalized
     vector is likely to be less than 0.5. */
  for (int i = 0; i < ShapeSize; i++) {
    intval = (int) (512.0 * vec[i]);
    vec[i] = (255 < intval) ? 255 : intval;
  }
   
  int sc_pca_size=36;
  ds->pca(sc_pca_size, sc_pca_avg, sc_pca_base);	

} 
  
void computeShape(DARY *dx, DARY *dy ,DARY *grad, DARY *ori, FeatureDescriptor *ds){

  DARY *edge = new DARY(dx->y(),dx->x());

  cannyEdges(dx, dy, grad, edge,   5, 15);
  //  edge->write("edge.pgm");cout << "OK"<< endl;getchar();

  float angle=ds->getAngle();
  if(angle<-M_PI || angle>M_PI)angle=0;

  ds->allocVec(ShapeSize);
  float *vec = ds->getVec();
  KeyLogPolSample(vec, edge, ori,  angle, SOriSize, SrSize, ScSize);
 //getchar();

  //KeySample(vec, edge, ori, angle, SOriSize, SLocSize);
  delete edge;
  NormalizeVect(vec, ShapeSize); 
  int intval, changed = FALSE;
  for (int i = 0; i < ShapeSize; i++)
    if (vec[i] > MaxIndexVal) { 
      vec[i] = MaxIndexVal;
       changed = TRUE;
     }
  if (changed)
    NormalizeVect(vec, ShapeSize); 
  /* Convert float vector to integer. 
     Assume largest value in normalized
     vector is likely to be less than 0.5. */
  for (int i = 0; i < ShapeSize; i++) {
    intval = (int) (512.0 * vec[i]);
    vec[i] = (255 < intval) ? 255 : intval;
  }
 
  //int sc_pca_size=36;
  //ds->pca(sc_pca_size, sc_pca_avg, sc_pca_base);	

} 


void computePcaSift(DARY *img, FeatureDescriptor *ds){
  

  float *tvec = new float[GPLEN];
  uint count=0;
  for(int j=1;j<PATCH_SIZE-1;j++){
    for(int i=1;i<PATCH_SIZE-1;i++){      
      tvec[count++]=img->fel[j][i+1]-img->fel[j][i-1];
      tvec[count++]=img->fel[j+1][i]-img->fel[j-1][i];
    }
  }

  NormalizeEVect(tvec, GPLEN, 1.0); 


  for (int i = 0; i < GPLEN; i++) {
    tvec[i] -= avgs[i];
  }

  ds->allocVec(PCALEN);//EIndexSize*EIndexSize*OriSize;
  float *vec = ds->getVec();

  for (int ldi = 0; ldi < PCALEN; ldi++) {
    for (int x = 0; x < GPLEN; x++)
      vec[ldi] += eigs[ldi][x] * tvec[x];
  }
  delete tvec;
}

void computeNON(DARY *img, FeatureDescriptor *ds){}
