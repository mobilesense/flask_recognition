#include "feature.h"
#include "../gauss_iir/gauss_iir.h"
#include "../util/util.h"
#include <sys/time.h> 

#include <algorithm>

#define ALPHA 0.05  
#define SCALE_LAP_THRES 5
#define MIN_SCALE 2 
#define MAX_SCALE 7
#define BORDER  3

/**
**/

void findAffineRegion(vector<DARY *> image,vector<DARY *> laplacian, vector<float> scale ,
		      vector<FeatureDescriptor*> &cor, int lmax);

int check(int x, int y, int width, int height, int border){
  if(x>border && y>border && x<width-border && y<height-border)return 1;
  else return 0;

}


void thresholdFeatures(vector<FeatureDescriptor*> &cor){
  
  vector<FeatureDescriptor*> cor_tmp;
  for(int i=0;i<(int)cor.size();i++){
    if(cor[i]->getFeatureness()>2000 && fabs(cor[i]->getLap())>15 && cor[i]->getExtr())    
      cor_tmp.push_back(cor[i]);        
  }  
  cor.clear();
  cor=cor_tmp;
}



/******************MULTI HESSIAN DETECTOR**************/

/******************MULTI HESSIAN DETECTOR**************/

/******************MULTI HESSIAN DETECTOR**************/

/******************MULTI HESSIAN DETECTOR**************/

/******************MULTI HESSIAN DETECTOR**************/

/******************MULTI HESSIAN DETECTOR**************/

/******************MULTI HESSIAN DETECTOR**************/

void laplacian(DARY *sm, DARY * lap){
  
  DARY *fxx = new DARY(sm->y(), sm->x());
  DARY *fyy = new DARY(sm->y(), sm->x());
  
  dXX9(sm,fxx);
  dYY9(sm,fyy);
  for(uint x=0;x<fxx->x();x++){
    for(uint y=0;y<fxx->y();y++){
      lap->fel[y][x]=fxx->fel[y][x]+fyy->fel[y][x];
    }
  }
  delete fxx;delete fyy;
  //lap->writePNG("lap.png");cout << "ci"<< endl;getchar();
}



/****** HESSIAN DETECTOR***********/
void hessian(DARY *image_in, DARY *hes, DARY *lap){
    
    
   int col_nb, row_nb;
   float A, B, AB, determinant;
   float C, C2;
   row_nb = image_in->y();
   col_nb = image_in->x(); 
    
   //printf("- Starting to calculate gradients \n");
   DARY  *fxx= new DARY(row_nb,col_nb);	
   DARY  *fyy= new DARY(row_nb,col_nb);	
   DARY  *fxy= new DARY(row_nb,col_nb);	
   dXX9(image_in, fxx);//fxx->writePNG("fxx.png");
   dYY9(image_in, fyy);//fyy->writePNG("fyy.png");
   dXY7(image_in, fxy);//fxy->writePNG("fxy.png");

   int row, col;
 
   for ( row = 0; row < row_nb; row++)
     for ( col = 0; col < col_nb; col++)
       {
	 /*        A = B = C = determinant = trace = 0.0;*/
        A = fxx->fel[row][col];
        B = fyy->fel[row][col];
        C = fxy->fel[row][col];
	C2=110*(C*C);// scaling factor to make equal amplitudes of fxx and fxy
	AB=(A*B);
        determinant = AB - (C2);
	lap->fel[row][col]=A+B;
        hes->fel[row][col] =(determinant);	
       }
   //hes->writePNG("har.png");cout << "ci"<< endl;getchar();
   delete fxx; delete fyy; delete fxy; 
}

void harris(DARY *img,DARY *har,DARY *har11,DARY *har12,DARY *har22){    
  
   int col_nb, row_nb;
   float A, B, C, determinant, trace, t1,t2;
   
   row_nb = img->y();
   col_nb = img->x();
   
   //   printf("- Starting to calculate gradients %d %d\n",col_nb, row_nb);
   DARY  *fx= new DARY(row_nb,col_nb);	
   dX6(img, fx);
   DARY  *fy= new DARY(row_nb,col_nb);	
   dY6(img, fy);
   DARY *fxy  = new DARY(row_nb, col_nb);
   int row, col;
   for ( row = 0; row < row_nb; row++)
       for ( col = 0; col < col_nb; col++){
	   t1 = fx->fel[row][col];
	   t2 = fy->fel[row][col];
	   fx->fel[row][col] = t1*t1;
	   fy->fel[row][col] = t2*t2;
	   fxy->fel[row][col] = t1*t2;
       }

   smoothSqrt(fx,har11);  delete fx;
   smoothSqrt(fy,har22);  delete fy;
   smoothSqrt(fxy,har12); delete fxy;

   for ( row = 0; row < row_nb; row++)
     for ( col = 0; col < col_nb; col++)
       {
	 /*        A = B = C = determinant = trace = 0.0;*/
        A = har11->fel[row][col];
        B = har22->fel[row][col];
        C = har12->fel[row][col];
        determinant = A * B - (C*C);
        trace = A + B;
        har->fel[row][col] = (determinant - ALPHA * (trace*trace));
       }
}


void harris(DARY *dx,DARY *dy,DARY *har,DARY *har11,DARY *har12,DARY *har22){    
  
   int col_nb, row_nb;
   float A, B, C, determinant, trace, t1,t2;
   
   row_nb = dx->y();
   col_nb = dx->x();
   
   //   printf("- Starting to calculate gradients %d %d\n",col_nb, row_nb);
   DARY  *fx= new DARY(row_nb,col_nb);	
   DARY  *fy= new DARY(row_nb,col_nb);	
   DARY *fxy  = new DARY(row_nb, col_nb);
   int row, col;
   for ( row = 0; row < row_nb; row++)
       for ( col = 0; col < col_nb; col++){
	   t1 = dx->fel[row][col];
	   t2 = dy->fel[row][col];
	   fx->fel[row][col] = t1*t1;
	   fy->fel[row][col] = t2*t2;
	   fxy->fel[row][col] = t1*t2;
       }

   smoothSqrt(fx,har11);  delete fx;
   smoothSqrt(fy,har22);  delete fy;
   smoothSqrt(fxy,har12); delete fxy;

   for ( row = 0; row < row_nb; row++)
     for ( col = 0; col < col_nb; col++)
       {
	 /*        A = B = C = determinant = trace = 0.0;*/
        A = har11->fel[row][col];
        B = har22->fel[row][col];
        C = har12->fel[row][col];
        determinant = A * B - (C*C);
        trace = A + B;
        har->fel[row][col] = (determinant - ALPHA * (trace*trace));
       }
}




float inline interpScale(float a, float b, float c, int i, vector<float> scale){
  float ds=interpPeak(a,b,c),sc;
  if(ds>0){
    sc=scale[i]+ds*(scale[i+1]-scale[i]);	
  }else {	
    sc=scale[i-1]+(1-ds)*(scale[i]-scale[i-1]);	
  }
  return sc;
}

int getLapMax(vector<DARY *> lap, vector<float> scale, int level, int minlev, int maxlev, float x, float y,float &sc, int &extr){

  vector<float> llap;
  float fx,fy;
  int flag=1;
  for(int i=0; i<(int)lap.size() && i<(level+maxlev) && flag;i++){
    fx=x/scale[i];
    fy=y/scale[i];
    if(fx>BORDER && fx<(lap[i]->x()-BORDER) && fy>BORDER && fy<(lap[i]->y()-BORDER)){
      llap.push_back(lap[i]->getValue(fx,fy));
      //cout << llap[llap.size()-1]<< " ";
    } else flag=0; 
  } 

  //cout << endl;
  if(level-minlev<2)minlev=level-2;
  //cout << level-minlev<< " max " << level+maxlev<< endl;
  for(int i=level-minlev;i<(level+maxlev) && i<(int)(llap.size()-2);i++){
    //cout << " li-1 "<< llap[i-1] << " "<< llap[i]<< " "<< llap[i+1]<< endl;
    //for local maximum or minimum of laplacian
    if(llap[i]>SCALE_LAP_THRES && llap[i-1]> llap[i-2] && llap[i]> llap[i-1] && llap[i]>llap[i+1] && llap[i+1]>llap[i+2]){
      sc=interpScale(llap[i-1],llap[i],llap[i+1],i,scale);
      extr=1;
      llap.clear();      
      return i;            
    }else if(llap[i]<-SCALE_LAP_THRES && llap[i-1]< llap[i-2] && llap[i]< llap[i-1] && llap[i]<llap[i+1] && llap[i+1]<llap[i+2]){
      sc=interpScale(llap[i-1],llap[i],llap[i+1],i,scale);
      extr=-1;
     llap.clear();      
      return i;            
    }  
  } 
  llap.clear();
  return -1;
}


void harris_lap(vector<DARY *> har,vector<DARY *> har11,vector<DARY *> har12,vector<DARY *> har22, vector<DARY *> lap, vector<float> scale, vector<FeatureDescriptor*>&features, float threshold, int type){
  
  FeatureDescriptor *cor;
  float act_pixel,fx,fy;
  int level;
  int extr;//(int)rint(GAUSS_CUTOFF*scale+2);
  float a,b,c,l1,l2,ea,fcol,frow,sc,dc=0,dr=0,thres=threshold;  
  int cbad=0,cgood=0;
  float aprox=0.0;
  if(type==1)aprox=1.0;
  for(uint i=1;i<har.size()-1;i++){
    for (uint row = BORDER ; row+BORDER  < har[i]->y(); row++){
      for (uint col =  BORDER; col+BORDER  < har[i]->x(); col++){
	act_pixel=har[i]->fel[row][col];
	if(act_pixel > har[i]->fel[row-1][col-1] &&
	   act_pixel > har[i]->fel[row-1][col] &&
	   act_pixel > har[i]->fel[row-1][col+1] &&
	   act_pixel > har[i]->fel[row][col-1] &&
	   act_pixel > har[i]->fel[row][col+1] &&
	   act_pixel > har[i]->fel[row+1][col-1] &&
	   act_pixel > har[i]->fel[row+1][col] &&
	   act_pixel > har[i]->fel[row+1][col+1] &&
	   har[i]->fel[row-1][col]>thres &&
	   har[i]->fel[row][col-1]>thres &&
	   har[i]->fel[row][col+1]>thres &&
	   har[i]->fel[row+1][col]>thres){
	  dc=aprox*interpPeak(har[i]->fel[row][col-1],act_pixel,har[i]->fel[row][col+1]);
	  dr=aprox*interpPeak(har[i]->fel[row-1][col],act_pixel,har[i]->fel[row+1][col]);
	  //cout << har[i]->fel[row][col-1] << " " << act_pixel << " " << har[i]->fel[row][col+1] << endl;
	  //cout <<"dc "<<  dc << endl;
	  //cout << har[i]->fel[row-1][col] << " " << act_pixel << " " << har[i]->fel[row+1][col] << endl;
	  //cout <<"dr "<<  dr << endl;
	  fcol=(float)col+dc; 
	  frow=(float)row+dr;
	  //cout<< col << "  "<< row << "  "<< scale[i]*col << " " << scale[i]*row<<  endl;
	  //cout << col << " " << row<< "  "<< scale[i]*fcol << " " << scale[i]*frow<<  endl;getchar();
	  fx=scale[i]*fcol;//get coordinates at scale level 0
	  fy=scale[i]*frow;
	  level=getLapMax(lap, scale, i, MIN_SCALE, MAX_SCALE, fx, fy, sc, extr);
	  if(level>0){
	    cgood++;
	    //cout << sc << endl; 
	    cor=new FeatureDescriptor(fx,fy, sc, act_pixel);
	    cor->setExtr(extr);
	    cor->setDerLev(i);
	    cor->setIntLev(level);
	    cor->setType(type); 
	    if(type==1){
	      a=har11[i]->fel[row][col];
	      b=har12[i]->fel[row][col];
	      c=har22[i]->fel[row][col];
	      //cout << a << " "<< b << " " << c << endl;
	      invSqrRoot(a,b,c,l2,l1,ea);		      
	      cor->setMi(a,b,b,c,l2,l2/l1,ea);
	    }else {
	      cor->setMi(1,0,0,1);
	    }
	    features.push_back(cor);
	  }else cbad++;
	}
      }
    }   
  }
  //cout << "cgood "<< cgood << " cbad "<< cbad << " all " << (cgood+cbad) << endl;
}

void setMask(int x, int y, int r, DARY *mask){
  int xs=((x-r)>=0)?(x-r):0;
  int ys=((y-r)>=0)?(y-r):0;
  int xe=((x+r)<(int)mask->x())?(x+r):(mask->x()-1);
  int ye=((y+r)<(int)mask->y())?(y+r):(mask->y()-1);

  //cout << xs << " " << ys << " " << xe << " " << ye << endl;

  for(int j=ys;j<=ye;j++){
    for(int i=xs;i<=xe;i++){
      mask->fel[j][i]=0;
    }
  }
  
}

void edge_lap(vector<DARY *> edge,vector<DARY *> mask,vector<DARY *> har11,vector<DARY *> har12,vector<DARY *> har22, vector<DARY *> lap, vector<float> scale, vector<FeatureDescriptor*>&features, float threshold, int radius, int type){
  
  FeatureDescriptor *cor;
  float act_pixel,fx,fy;
  int level,sx=0,sy=0;
  int extr;//(int)rint(GAUSS_CUTOFF*scale+2);
  float a,b,c,l1,l2,ea,el,sc,thres=threshold/1.0;  
  int cbad=0,cgood=0;
  for(uint i=1;i<edge.size()-1;i++){
    for (uint row = BORDER ; row+BORDER  < edge[i]->y(); row++){
      for (uint col =  BORDER; col+BORDER  < edge[i]->x(); col++){
	act_pixel=edge[i]->fel[row][col];
	if(act_pixel > thres && mask[i]->fel[row][col]>0){
	  fx=scale[i]*col;//get coordinates at scale level 0
	  fy=scale[i]*row;
	  level=getLapMax(lap, scale, i, MIN_SCALE, MAX_SCALE, fx, fy, sc, extr);
	  if(level>0){	    
	    sx=(int)(0.5+col/scale[level-i]);
	    sy=(int)(0.5+row/scale[level-i]); 
	    if(mask[i]->fel[row][col]>0){cgood++;
	      //cout << sc << endl; 
	      cor=new FeatureDescriptor(fx,fy, sc, act_pixel);
	      cor->setExtr(extr);
	      cor->setDerLev(i);
	      cor->setIntLev(level);
	      cor->setType(type);
	      setMask(col,row,radius,mask[i]);
	      if(type==1){
		a=har11[i]->fel[row][col];
		b=har12[i]->fel[row][col];
		c=har22[i]->fel[row][col];
		invSqrRoot(a,b,c,l2,l1,ea);		      
		el=l2/l1;
		//cout << a << " "<< b << " " << c << " " << l1 << " " << l2 << endl;getchar();
		cor->setMi(a,b,b,c,l1,el,ea);
	      }else {
		cor->setMi(1,0,0,1);
	      }
	      features.push_back(cor);
	    }
	  }else cbad++;
	}
      }
    }  
    //mask[i]->writePNG("mask1.png");cout << "mask"<< endl;getchar(); 
 
  }
  //cout << "cgood "<< cgood << " cbad "<< cbad << " all " << (cgood+cbad) << endl;
	    
}


void sedge_lap(vector<DARY *> edge,vector<DARY *> mask,vector<DARY *> har11,vector<DARY *> har12,vector<DARY *> har22, vector<DARY *> lap, vector<float> scale, vector<FeatureDescriptor*>&features, float threshold, int radius, float eratio, int type){
  
  FeatureDescriptor *cor;
  float act_pixel,fx,fy,fcol,frow;
  int level,sx=0,sy=0;
  int extr;//(int)rint(GAUSS_CUTOFF*scale+2);
  float a,b,c,l1,l2,ea,el=1,sc,thres=threshold;  
  int cbad=0,cgood=0;
  for(uint i=1;i<edge.size()-1;i++){ 
    for (uint row = BORDER ; row+BORDER  < edge[i]->y(); row++){
      for (uint col =  BORDER; col+BORDER  < edge[i]->x(); col++){
	act_pixel=edge[i]->fel[row][col];
	if(act_pixel > thres &&
	   act_pixel > edge[i]->fel[row-1][col-1] &&
	   act_pixel > edge[i]->fel[row-1][col] &&
	   act_pixel > edge[i]->fel[row-1][col+1] &&
	   act_pixel > edge[i]->fel[row][col-1] &&
	   act_pixel > edge[i]->fel[row][col+1] &&
	   act_pixel > edge[i]->fel[row+1][col-1] &&
	   act_pixel > edge[i]->fel[row+1][col] &&
	   act_pixel > edge[i]->fel[row+1][col+1] &&
	   edge[i]->fel[row][col+1] > thres &&
	   edge[i]->fel[row+1][col] > thres &&
	   edge[i]->fel[row][col-1] > thres &&
	   edge[i]->fel[row-1][col] > thres &&
	   mask[i]->fel[row][col]>0){
	  fcol=(float)col+interpPeak(edge[i]->fel[row][col-1],act_pixel,edge[i]->fel[row][col+1]);
	  frow=(float)row+interpPeak(edge[i]->fel[row-1][col],act_pixel,edge[i]->fel[row+1][col]);
	  fx=scale[i]*fcol;//get coordinates at scale level 0
	  fy=scale[i]*frow; 
	  level=getLapMax(lap, scale, i, MIN_SCALE, MAX_SCALE, fx, fy, sc, extr);
	  //cout << i << " " << level <<  endl;
	  if(level>0){
	    sx=(int)(0.5+fcol/scale[level-i]);
	    sy=(int)(0.5+frow/scale[level-i]); 
	    //cout << act_pixel<< " sx "<< fx<< " sy "<< fy << " i "<< i << " level "<< level << " col " << col << " row " << row << endl;
	    if(mask[i]->fel[row][col]>0){cgood++;
	      //cout << sc << endl; 
	      cor=new FeatureDescriptor(fx,fy, sc, act_pixel);
	      cor->setExtr(extr);
	      cor->setDerLev(i);
	      cor->setIntLev(level);
	      cor->setType(type);
	      setMask(col,row,radius,mask[i]);//mask[i]->writePNG("mask.png");cout << "mask"<< endl;getchar();
	      if(type==1){
		a=har11[i]->fel[row][col];
		b=har12[i]->fel[row][col];
		c=har22[i]->fel[row][col];
		invSqrRoot(a,b,c,l2,l1,ea);
		el=l2/l1;
		//cout << a << ", "<< b << ", "<< b << ", " << c << " " << l1 << " " << l2 << endl;getchar();
		cor->setMi(a,b,b,c,l1,el,ea);
	      }else {
		cor->setMi(1,0,0,1);
	      }
	      if(el>eratio)
		features.push_back(cor);
	      else delete cor;
	    }
	  }else cbad++;
	}
      }
    }
    //mask[i]->writePNG("mask1.png");cout << "mask"<< endl;getchar(); 
  }
  //cout << "cgood "<< cgood << " cbad "<< cbad << " all " << (cgood+cbad) << endl;
	    
}

void buildScaleSpace(DARY *img, vector<DARY *> &sm, vector<float> &sc, float step){
  
  int sx=img->x(),sy=img->y();
  sc.push_back(1);
  sm.push_back(new DARY(sy,sx));
  while(sx>12 && sy>12){
    sx=(int)(0.5+sx/step);sy=(int)(0.5+sy/step);    
    sm.push_back(new DARY(sy,sx));
    sc.push_back(step*sc[sc.size()-1]);
  }
  // sm[0]->set(img);
  if(sm.size()>1) {
    smooth5(img,sm[0]);
    sm[1]->scale(sm[0],step,step);
    float sc2=step*step;
    for(uint i=2;i<sm.size()-1;i+=2){
      sm[i]->scale(sm[i-2],sc2,sc2); 
      sm[i+1]->scale(sm[i-1],sc2,sc2); 
    }  
  }
}

void buildScaleSpaceSqrt2(DARY *img, vector<DARY *> &sm, vector<float> &sc, int max_size){
  int sx=img->x(),sy=img->y();
  sc.push_back(1);
  sm.push_back(new DARY(sy,sx));

  float sqrt2=sqrt(2);

  smooth5(img,sm[0]); 

  for(int i=1;i<max_size;i++) {
    sc.push_back(pow(sqrt2,i));    
    if(i%2==1) {
      sx=(int)(sx/sqrt2);
      sy=(int)(sy/sqrt2);
    } else {
      sx=(int)(sm[i-2]->x()/2);
      sy=(int)(sm[i-2]->y()/2);
    }
    if(!(sx>12 && sy>12)) break;
    
    sm.push_back(new DARY(sy,sx));
    if(i%2==1) 
      sm[i]->scale(sm[i-1],sqrt2,sqrt2);
    else
      sm[i]->scaleHalf(sm[i-2]);
  }

  

}

/*************************************************/
void computeDescriptors( vector<DARY *> sm, vector<float> sc, vector<FeatureDescriptor*>&features, int DESC, int noangle){

  int level=0;
  float x,y;
  float c_scale=PATCH_SIZE/2;
  float angle;
  float scal; 
  vector<FeatureDescriptor*> desc;
  initPatchMask(PATCH_SIZE);
  DARY * patch = new DARY(PATCH_SIZE,PATCH_SIZE);   
  float scale_factor=(PATCH_SIZE/2.88);
  for(uint i=0;i<features.size();i++){
    level=features[i]->getIntLev()-2;    
    angle=features[i]->getAngle();
    if(noangle){
      angle=0;
      features[i]->setAngle(0.0);
    }
    x=(features[i]->getX()/sc[level]);
    y=(features[i]->getY()/sc[level]);

    features[i]->setScale(features[i]->getScale()*scale_factor);//(patch_size/2)*scale/(step*step)

    DARY *imgbs= normalizeAffine(sm[level], x, y, c_scale, 
				 angle, features[i]->getMi11(), features[i]->getMi12(), 
				 features[i]->getMi21(),features[i]->getMi22(),
				 scal,patch, 1.0);
    //patch->writePNG("patch1.png");cout <<i <<  " OK1 "<<features[i]->getX() << " " << features[i]->getY() <<  endl;getchar();
    //if(!(i%100))cout << "\raff_descriptor " << i << " of " << features.size()<< "    "<< flush;
    computeAffineDescriptor(imgbs, patch, scal,DESC,features[i],desc);
    delete imgbs;
  }
  //  cout << desc.size()<< endl;
  delete patch;
  deleteDescriptors(features);
  features=desc;  
} 



/*************************************************/
void computeDescriptors( vector<DARY *> dx,vector<DARY *> dy, 
			 vector<float> sc, vector<FeatureDescriptor*>&features, int DESC, int noangle){

  int level=0;
  int x,y;
  vector<FeatureDescriptor*> desc;
  initPatchMask(PATCH_SIZE);
  DARY * dxpatch = new DARY(PATCH_SIZE,PATCH_SIZE,0.0);   
  DARY * dypatch = new DARY(PATCH_SIZE,PATCH_SIZE,0.0);   

  if(noangle){
    for(uint i=0;i<features.size();i++){
      features[i]->setAngle(0.0);
    }
  }
  float scale_factor=(PATCH_SIZE/2.88);

  for(uint i=0;i<features.size();i++){
    level=features[i]->getIntLev()-2;    //scale/(step*step)=1.44
    x=(int)(0.5+features[i]->getX()/sc[level]);
    y=(int)(0.5+features[i]->getY()/sc[level]);


    features[i]->setScale(features[i]->getScale()*scale_factor);//(patch_size/2)*scale/(step*step)
    dxpatch->crop(dx[level],x,y);
    dypatch->crop(dy[level],x,y);
    //patch_mask->fel[0][0]=-0.00001;patch_mask->writePNG("mask.png");
    //if(features[i]->getScale()>7){dxpatch->writePNG("patch2.png");cout <<i <<  " OK 2 "<<features[i]->getX() << " " << features[i]->getY() <<" " << features[i]->getScale() <<  endl;getchar();}
    //if(!(i%100))cout << "\rdescriptor " << i << " of " << features.size()<< "    "<< flush;

    computeDescriptor( dxpatch, dypatch, DESC, features[i], desc);
  }
  //cout << desc.size()<< endl;
  delete dxpatch; delete dypatch;
  deleteDescriptors(features);
  features=desc;  
} 


void fastEdge(DARY* img, vector<FeatureDescriptor*>&features, 
		      float threshold,float step, int DESC, int noangle){

  vector<DARY *> sm; 
  vector<DARY *> lap; 
  vector<DARY *> edge; 
  vector<DARY *> dx; 
  vector<DARY *> dy; 
  vector<DARY *> grad; 
  vector<DARY *> ori; 
  vector<DARY *> mask; 
  vector<float> sc; 

  buildScaleSpace(img, sm, sc,step);
  for(uint i=0;i<sm.size();i++){
    lap.push_back(new DARY(sm[i]->y(),sm[i]->x()));
    laplacian(sm[i],lap[i]);

    if(i>0 && i<sm.size()-1){
      edge.push_back(new DARY(sm[i]->y(),sm[i]->x(),0.0));//
      dx.push_back(new DARY(sm[i]->y(),sm[i]->x()));//
      dy.push_back(new DARY(sm[i]->y(),sm[i]->x()));//
      grad.push_back(new DARY(sm[i]->y(),sm[i]->x()));//
      ori.push_back(new DARY(sm[i]->y(),sm[i]->x()));//
      mask.push_back(new DARY(sm[i]->y(),sm[i]->x(),255.0));//
      gradAngle(sm[i],dx[i],dy[i],grad[i],ori[i]);    
      cannyEdges(dx[i],dy[i], grad[i], edge[i], 15, 30);
    //char nm[512];sprintf(nm,"sm%d.png",i);sm[i]->writePNG(nm);sprintf(nm,"grad%d.png",i);grad[i]->writePNG(nm);
    //sprintf(nm,"edge%d.png",i);edge[i]->writePNG(nm);sprintf(nm,"lap%d.png",i);lap[i]->writePNG(nm);
    }else{
      edge.push_back(new DARY(1,1,0.0));//
      dx.push_back(new DARY(1,1));//
      dy.push_back(new DARY(1,1));//
      grad.push_back(new DARY(1,1));//
      ori.push_back(new DARY(1,1));//      
      mask.push_back(new DARY(1,1));//      
    }
  }
  edge_lap(edge,mask,edge,edge,edge,lap,sc,features,threshold,2,0);
 
  if(DESC==JLA || DESC==CC || DESC==CF || DESC==KOEN || DESC==SPIN)
    computeDescriptors( sm, sc, features, DESC, noangle);
  else computeDescriptors( dx, dy, sc, features, DESC, noangle);


  for(uint i=0;i<sm.size();i++){
    delete sm[i];delete lap[i];delete edge[i];delete dx[i];delete dy[i];delete grad[i];delete ori[i];delete mask[i];
  }
  sc.clear();
  sm.clear();
  lap.clear();  
  edge.clear(); 
  dx.clear();
  dy.clear();
  grad.clear();
  ori.clear();
  mask.clear();
  
}


void harEdge(DARY* edge,DARY *har){

  for (uint row = 0 ; row  < edge->y(); row++){
    for (uint col =  0; col  < edge->x(); col++){
      if(edge->fel[row][col]<=0 || har->fel[row][col]<=100)
	har->fel[row][col]=0;
    }
  }
}


void fastSEdge(DARY* img, vector<FeatureDescriptor*>&features, 
		      float threshold,float step, int DESC, int noangle){

  vector<DARY *> sm; 
  vector<DARY *> lap; 
  vector<DARY *> edge; 
  vector<DARY *> hes; 
  vector<DARY *> har; 
  vector<DARY *> har11; 
  vector<DARY *> har12; 
  vector<DARY *> har22; 
  vector<DARY *> dx; 
  vector<DARY *> dy; 
  vector<DARY *> grad; 
  vector<DARY *> ori; 
  vector<DARY *> mask; 
  vector<float> sc; 
 
  //cout << "start "<< endl;

  buildScaleSpace(img, sm, sc,step);

  //laplacian(sm[0],lap[0]);
  for(uint i=0;i<sm.size();i++){
    lap.push_back(new DARY(sm[i]->y(),sm[i]->x()));
    hes.push_back(new DARY(sm[i]->y(),sm[i]->x()));//
    hessian(sm[i],hes[i],lap[i]);//

    if(i>0 && i<sm.size()-1){
      edge.push_back(new DARY(sm[i]->y(),sm[i]->x(),0.0));//
      dx.push_back(new DARY(sm[i]->y(),sm[i]->x()));//
      dy.push_back(new DARY(sm[i]->y(),sm[i]->x()));//
      grad.push_back(new DARY(sm[i]->y(),sm[i]->x()));//
      ori.push_back(new DARY(sm[i]->y(),sm[i]->x()));//
      gradAngle(sm[i],dx[i],dy[i],grad[i],ori[i]);    
      cannyEdges(dx[i],dy[i], grad[i], edge[i], 15, 30);
      har.push_back(new DARY(sm[i]->y(),sm[i]->x()));//
      har11.push_back(new DARY(sm[i]->y(),sm[i]->x()));//
      har12.push_back(new DARY(sm[i]->y(),sm[i]->x()));//
      har22.push_back(new DARY(sm[i]->y(),sm[i]->x()));//
      mask.push_back(new DARY(sm[i]->y(),sm[i]->x(),255.0));//
      harris(dx[i],dy[i],har[i],har11[i],har12[i],har22[i]);
      //harEdge(edge[i],har[i]);
      //char nm[512];
      //sprintf(nm,"sm%d.png",i);sm[i]->writePNG(nm);sprintf(nm,"grad%d.png",i);grad[i]->writePNG(nm);
      //sprintf(nm,"edge%d.png",i);edge[i]->writePNG(nm);//sprintf(nm,"lap%d.png",i);lap[i]->writePNG(nm);
      //sprintf(nm,"har%d.png",i);har[i]->writePNG(nm);
      //sprintf(nm,"hes%d.png",i);hes[i]->writePNG(nm);
      //sprintf(nm,"har2%d.png",i);har[i]->writePNG(nm);
    }else{
      edge.push_back(new DARY(1,1,0.0));//
      dx.push_back(new DARY(1,1));//
      dy.push_back(new DARY(1,1));//
      grad.push_back(new DARY(1,1));//
      ori.push_back(new DARY(1,1));//      
      har.push_back(new DARY(1,1));//      
      har11.push_back(new DARY(1,1));//      
      har12.push_back(new DARY(1,1));//      
      har22.push_back(new DARY(1,1));//      
      mask.push_back(new DARY(1,1));//      
    }
  }
   
  sedge_lap(hes,mask,har11,har12,har22,lap,sc,features,threshold,1,0.1,1);//cout<< "hes " << features.size()<< endl;
  sedge_lap(har,mask,har11,har12,har22,lap,sc,features,threshold,1,0.1,1);//cout<< "haredge " << features.size()<< endl;
  //sedge_lap(edge,mask,har11,har12,har22,lap,sc,features,threshold,3,0,1);//cout<< "maxedge " << features.size()<< endl;
  threshold=2;
  edge_lap(edge,mask,har11,har12,har22,lap,sc,features,threshold,2,1);//cout<< "edge2 " << features.size()<< endl;
 
 
  //cout << "done "<< endl;
  // for(uint i=0;i<features.size();i++){
  //  img->fel[(int)features[i]->getY()][(int)features[i]->getX()]=127;
    
  //}
  //img->writePNG("points.png");

  if(DESC==JLA || DESC==CC || DESC==CF || DESC==KOEN || DESC==SPIN)
    computeDescriptors( sm, sc, features, DESC, noangle);
  else computeDescriptors( dx, dy, sc, features, DESC, noangle);
  //cout << "done 2"<< endl;


  for(uint i=0;i<sm.size();i++){
    delete sm[i];delete lap[i];delete edge[i];delete dx[i];delete dy[i];delete grad[i];delete ori[i];
    delete har[i];delete har11[i];delete har12[i];delete har22[i];delete mask[i];delete hes[i];
  }
  sc.clear();
  sm.clear();
  lap.clear();  
  edge.clear(); 
  dx.clear();
  dy.clear();
  grad.clear();
  ori.clear();
  
  hes.clear();  
  har.clear();  
  har11.clear();  
  har12.clear();  
  har22.clear(); 
  mask.clear();
}

 void fastHarrisHessian(DARY *img, vector<FeatureDescriptor*>&corners, float threshold, float step, int DESC, int aff, int noangle){
  
  vector<DARY *> sm; 
  vector<DARY *> lap; 
  vector<DARY *> hes; 
  vector<float> sc;
  vector<DARY *> dx; 
  vector<DARY *> dy; 
  vector<DARY *> har; 
  vector<DARY *> har11; 
  vector<DARY *> har12; 
  vector<DARY *> har22; 

  buildScaleSpace(img, sm, sc,step);
  for(uint i=0;i<sm.size();i++){
    lap.push_back(new DARY(sm[i]->y(),sm[i]->x()));//
    hes.push_back(new DARY(sm[i]->y(),sm[i]->x()));//
    hessian(sm[i],hes[i],lap[i]);//
    dx.push_back(new DARY(sm[i]->y(),sm[i]->x()));//
    dy.push_back(new DARY(sm[i]->y(),sm[i]->x()));//
    //  dX2(sm[i], dx[i]);dY2(sm[i], dy[i]);
    dX6(sm[i], dx[i]);dY6(sm[i], dy[i]);
    har.push_back(new DARY(sm[i]->y(),sm[i]->x()));
    har11.push_back(new DARY(sm[i]->y(),sm[i]->x()));
    har12.push_back(new DARY(sm[i]->y(),sm[i]->x()));
    har22.push_back(new DARY(sm[i]->y(),sm[i]->x()));
    harris(dx[i],dy[i],har[i],har11[i],har12[i],har22[i]);
    //char name[512];sprintf(name,"hes%d.png",i+2);sm[i]->writePNG(name);
    //cout <<i << " "<< sc[i] << endl; 
  }
  
  harris_lap(hes,har11,har12,har22,lap, sc, corners,threshold,1);  
  harris_lap(har,har11,har12,har22,lap,sc,corners,threshold,1);

  if(aff>1)findAffineRegion(sm,lap,sc,corners,aff);


  if(DESC==JLA || DESC==CC || DESC==CF || DESC==KOEN || DESC==SPIN || aff)
    computeDescriptors( sm, sc, corners, DESC, noangle);
  else computeDescriptors( dx, dy, sc, corners, DESC, noangle);



  for(uint i=0;i<sm.size();i++){
    delete sm[i];delete lap[i];delete hes[i];delete har[i];delete har11[i];delete har12[i];delete har22[i];
    delete dx[i];delete dy[i];
  }
  sc.clear();
  sm.clear();
  lap.clear();  
  har.clear();  
  hes.clear();  
  har11.clear();  
  har12.clear();  
  har22.clear();  
  hes.clear();  
  dx.clear();
  dy.clear();
}


struct CmpCor {
  bool operator()(FeatureDescriptor*f1,FeatureDescriptor*f2) {
    return f1->featureness > f2->featureness;
  }
}; 

/*double getmillisecs() {
  struct timeval tv;
  gettimeofday(&tv,NULL);
  return tv.tv_sec*1e3 +tv.tv_usec*1e-3;
}
*/

void fastHessian(DARY *img, vector<FeatureDescriptor*>&features, float threshold, float step, int DESC, int aff, int noangle, int max_desc ){

  //double t0,t1,t11,t2,t3,t4;
  
  vector<DARY *> sm; 
  vector<DARY *> lap; 
  vector<DARY *> hes; 
  vector<float> sc;
  vector<DARY *> dx; 
  vector<DARY *> dy; 

  //t0=getmillisecs();


  buildScaleSpace(img, sm, sc,step);
  //t11 = getmillisecs();
  for(uint i=0;i<sm.size();i++){
    lap.push_back(new DARY(sm[i]->y(),sm[i]->x()));//
    hes.push_back(new DARY(sm[i]->y(),sm[i]->x()));//
    hessian(sm[i],hes[i],lap[i]);//
    dx.push_back(new DARY(sm[i]->y(),sm[i]->x()));//
    dy.push_back(new DARY(sm[i]->y(),sm[i]->x()));//
    dX2(sm[i], dx[i]);dY2(sm[i], dy[i]);
    //char name[512];sprintf(name,"hes%d.png",i+2);hes[i]->writePNG(name);
    //cout <<i << " "<< sc[i] << endl; 
  }
  //t1=getmillisecs();
  
  harris_lap(hes,hes,hes,hes,lap, sc, features,threshold,2);  

  //t2=getmillisecs();

  sort(features.begin(),features.end(),CmpCor());
  if(max_desc && (max_desc<(int)features.size()))
    features.erase(features.begin() + max_desc, features.end()) ;

  //printf( "%lf %lf\n", features[0]->featureness, features[1]->featureness);
  
  if(aff>1)findAffineRegion(sm,lap,sc,features,aff);

   cout << " interest points "<< features.size()<< endl;

  //t3=getmillisecs();


  if(DESC==JLA || DESC==CC || DESC==CF || DESC==KOEN || DESC==SPIN || aff)
    computeDescriptors( sm, sc,features , DESC, noangle);
  else computeDescriptors( dx, dy, sc, features, DESC, noangle);

  for(uint i=0;i<sm.size();i++){
    delete sm[i];delete lap[i];delete hes[i];delete dx[i];delete dy[i];
  }

  //t4=getmillisecs();

  //printf("times %.3f %.3f %.3f %.3f %.3f\n",t11-t0, t1-t11, t2-t1, t3-t2, t4-t3);

  sc.clear();
  sm.clear();
  lap.clear();  
  hes.clear();  
  dx.clear();
  dy.clear();
}


void fastDense(DARY *img, vector<FeatureDescriptor*>&features, 
               int DESC, int xstep, int ystep, int max_pyrlevels ){

  
  vector<DARY *> sm; 
  vector<float> sc;
  vector<DARY *> dx; 
  vector<DARY *> dy; 

  buildScaleSpaceSqrt2(img, sm, sc, max_pyrlevels);

  for(uint i=0;i<sm.size();i++){
    dx.push_back(new DARY(sm[i]->y(),sm[i]->x()));//
    dy.push_back(new DARY(sm[i]->y(),sm[i]->x()));//
    dX2(sm[i], dx[i]);dY2(sm[i], dy[i]);
  }

  vector<FeatureDescriptor*> desc;

  initPatchMask(PATCH_SIZE);
  DARY * dxpatch = new DARY(PATCH_SIZE,PATCH_SIZE,0.0);   
  DARY * dypatch = new DARY(PATCH_SIZE,PATCH_SIZE,0.0);   

  DARY *gpatch = new DARY(PATCH_SIZE,PATCH_SIZE);
  DARY *opatch = new DARY(PATCH_SIZE,PATCH_SIZE);      
  
  for(int level=0;level<sm.size();level++) {
    int xd=dx[level]->x();
    int yd=dx[level]->y();
    float scale=sc[level];

    for(int y=0;y<yd;y+=ystep) {
      if(y-PATCH_SIZE/2<0 || y+PATCH_SIZE/2>yd) continue;

      for(int x=0;x<xd;x+=ystep) {
        if(x-PATCH_SIZE/2<0 || x+PATCH_SIZE/2>xd) continue;
        
        FeatureDescriptor *ds=new FeatureDescriptor();
        ds->init();
        ds->x=x*scale;
        ds->y=y*scale;
        ds->c_scale=scale;
        ds->angle=0;
        ds->featureness=0;

        /*                 computeDescriptor( dxpatch, dypatch, DESC, features[i], desc);*/
        
                
        dxpatch->crop(dx[level],x,y);
        dypatch->crop(dy[level],x,y);       
        gradAngle(dxpatch,dypatch,gpatch,opatch);  /* maybe not translation invariant */
        int ok=compDesc(dxpatch,dypatch,gpatch,opatch,ds,DESC);

        if(!ok) {
          delete ds;
          continue;
        }

        features.push_back(ds);
/*
        printf("desc=[");
        for(int j=0; j<ds->size; j++) printf("%g ",ds->vec[j]);
        printf("]\n");
  */      
      }
    }
  }
  //cout << desc.size()<< endl;
  delete gpatch;delete opatch;
  delete dxpatch; delete dypatch;


  for(uint i=0;i<sm.size();i++){
    delete sm[i];delete dx[i];delete dy[i];
  }

  sc.clear();
  sm.clear();
}



void fastHarris(DARY *img, vector<FeatureDescriptor*>&features, float threshold, float step, int DESC, int aff, int noangle){

  vector<DARY *> sm; 
  vector<DARY *> lap; 
  vector<DARY *> har; 
  vector<DARY *> har11; 
  vector<DARY *> har12; 
  vector<DARY *> har22; 
  vector<DARY *> dx; 
  vector<DARY *> dy; 
  vector<float> sc;
  
  buildScaleSpace(img, sm, sc,step);

  for(uint i=0;i<sm.size();i++){
      lap.push_back(new DARY(sm[i]->y(),sm[i]->x()));
      dx.push_back(new DARY(sm[i]->y(),sm[i]->x()));
      dy.push_back(new DARY(sm[i]->y(),sm[i]->x()));
      //dX2(sm[i], dx[i]);dY2(sm[i], dy[i]);
      dX6(sm[i], dx[i]);dY6(sm[i], dy[i]);

      if(i>0 && i< sm.size()-1){
	har.push_back(new DARY(sm[i]->y(),sm[i]->x()));
	har11.push_back(new DARY(sm[i]->y(),sm[i]->x()));
	har12.push_back(new DARY(sm[i]->y(),sm[i]->x()));
	har22.push_back(new DARY(sm[i]->y(),sm[i]->x()));
	harris(dx[i],dy[i],har[i],har11[i],har12[i],har22[i]);
      }else {
	har.push_back(new DARY(1,1));
	har11.push_back(new DARY(1,1));
	har12.push_back(new DARY(1,1));
	har22.push_back(new DARY(1,1));
      }
      laplacian(sm[i],lap[i]);
      //char nm[512];sprintf(nm,"sm%d.png",i);sm[i]->writePNG(nm);
      //sprintf(nm,"har%d.png",i);har[i]->writePNG(nm);sprintf(nm,"lap%d.png",i);lap[i]->writePNG(nm);
  }

  harris_lap(har,har11,har12,har22,lap,sc,features,threshold,1);

  if(aff>1)findAffineRegion(sm,lap,sc,features,aff);

  cout << " interest points "<< features.size()<< endl;

  if(DESC==JLA || DESC== SHAPE || DESC==CC || DESC==CF || DESC==KOEN || DESC==SPIN || aff)
    computeDescriptors( sm, sc, features, DESC, noangle);
  else computeDescriptors( dx, dy, sc, features, DESC, noangle);


  for(uint i=0;i<sm.size();i++){
    delete sm[i];delete lap[i];delete har[i];delete har11[i];delete har12[i];delete har22[i];
    delete dx[i];delete dy[i];
  }
  sc.clear();
  sm.clear();
  lap.clear();  
  har.clear();  
  har11.clear();  
  har12.clear();  
  har22.clear();  
  dx.clear();  
  dy.clear();  
}





/**********************************AFFINE**********************************************/


void getMi(DARY *img, float &a, float &b, float &c, float x, float y){
    int row_nb= img->y();
    int col_nb= img->y();    
    float  t1,t2;
    DARY  *fx= new DARY(row_nb, col_nb);	
    DARY  *fy= new DARY(row_nb, col_nb);	
    DARY  *fxy  = new DARY(row_nb, col_nb);
    dX6(img, fx);
    dY6(img, fy);
    for (int row = 0; row < row_nb; row++)
	for (int col = 0; col < col_nb; col++){
	    t1 = fx->fel[row][col];
	    t2 = fy->fel[row][col];
	    fx->fel[row][col] = t1*t1;
	    fy->fel[row][col] = t2*t2;
	    fxy->fel[row][col] = t1*t2;
	}          
    a= smoothf((int)x,(int)y, fx, 3);
    c= smoothf((int)x,(int)y, fy, 3);
    b= smoothf((int)x,(int)y, fxy, 3);
    
    delete fx;delete fy;delete fxy;
}


int fastfindAffineRegion(vector<DARY *> image,vector<DARY *> laplacian,vector<float>scale, FeatureDescriptor * cor, int lmax){

  int level = cor->getDerLev();
  float pointx=cor->getX()/scale[level];
  float pointy=cor->getY()/scale[level];
  int sizex=19;//2*(1.44*GAUSS_CUTOFF+3)+1;
  int l,go_on=1;
  float l1=1,l2=1,ea=0;
  float eigen_ratio_act=0.1,eigen_ratio_bef=0.1,deigen=0.02;
  //static float nba=0; 
  DARY *img=new DARY(sizex,sizex);
  DARY *cimg=image[level];
  float a,b,c,u11=cor->getMi11(),u12=cor->getMi12(),u21=cor->getMi12(),u22=cor->getMi22(),u11t,u12t;
  getEigen(u11,u12,u21,u22,l1,l2);
  eigen_ratio_act=1-l2/l1;
 
  //cout << pointx << " " << pointy << " level " << level << " " << eigen_ratio_act <<endl;getchar();
  for(l=0;l<lmax && go_on;l++){ //estimate affine structure
    img->interpolate(cimg,pointx,pointy,u11,u12,u21,u22);
    //img->writePNG("img.png");
    getMi(img, a,b,c, sizex>>1, sizex>>1);
    //cout <<l1 <<"  " << l2 <<" a " << a << " b " <<b << "  c " << c << endl;
    invSqrRoot(a,b,c,l2,l1,ea);
    
    eigen_ratio_bef=eigen_ratio_act;
    eigen_ratio_act=1-l2/l1;
    
    u11t=u11;u12t=u12;
    u11=a*u11t+b*u21;
    u12=a*u12t+b*u22;
    u21=b*u11t+c*u21;
    u22=b*u12t+c*u22;
    
    //cout << u11 << " "<< u12 << endl;
    //cout << u21 << " "<< u22 << endl;
    //cout << endl << l1 << " "<< l2 << endl;
    getEigen(u11,u12,u21,u22,l1,l2); 
    
    if(l>15 || (l1/l2>6) || (l2/l1>6)) {
      delete img;
      return 0;
    }
    
    if(eigen_ratio_act<deigen && eigen_ratio_bef<deigen)go_on=0;       
  }
  delete img;
  //cout <<"eigen_ratio_act "<<eigen_ratio_act<<"  "<< l1<< " "<< l2 << " "<< l<< endl;getchar();
  
  cor->setMi(u11,u12,u21,u22,l1,l2,ea);  
  return l;    
}


void findAffineRegion(vector<DARY *> image,vector<DARY *> laplacian, vector<float>scale,vector<FeatureDescriptor*> &cor, int lmax){
  for(int i=0;i<(int)cor.size();i++){        
    int l=fastfindAffineRegion(image,laplacian,scale,cor[i],lmax);    
    if(l!=0){
      //cout<<"\r  cor  "<<i<<" of "<< size << "  "<<cor[i]->getDerLev()<< "  " << cor[i]->getX() << "  " << cor[i]->getY()<<"  yes  "<< flush;
    }else { 
      //cout<<"\r  cor  "<<i<<" of "<< size << "  "<<cor[i]->getDerLev()<< "  " << cor[i]->getX() << "  " << cor[i]->getY()<<"  no  "<< flush;
      cor.erase((std::vector<FeatureDescriptor*>::iterator)&cor[i]);i--;}    
  }    
}





void extractFeatures(DARY *image, vector<FD*> &features, int detector, int descriptor, int affine, int noangle, float thres, int max_desc ){
  if(detector>0){
    if(detector==HARRIS){
      fastHarris(image, features, thres, 1.2, descriptor, affine, noangle);      
    }else if(detector==HESSIAN){
      fastHessian(image, features, thres, 1.2, descriptor, affine, noangle, max_desc); 
    }else if(detector==EDGE){     
      fastEdge(image, features, thres, 1.2, descriptor, noangle);  
    }else if(detector==HARHES){
      fastHarrisHessian(image, features, thres, 1.2, descriptor, affine, noangle);                  
    }else if(detector==SEDGE){ 
      fastSEdge(image, features, thres, 1.2, descriptor, noangle);                  
    }
  }
}


/**********************DETECT FEATURES*****************************/
/**********************DETECT FEATURES*****************************/


void getBoundingBox(DARY *mask,int &xc, int &yc,int &width, int &height){

  uint xs=mask->x(),ys=mask->y(),xe=0,ye=0;
  for(uint j=0;j<mask->y();j++){
    for(uint i=0;i<mask->x();i++){
      if(mask->bel[j][i]>0){
	if(xs>i)xs=i;
	if(xe<i)xe=i;
	if(ys>j)ys=j;
	if(ye<j)ye=j;
      }
    }
  }
  //xc=(xs+xe)/2;
  //yc=(ys+ye)/2;
  width=(xe-xs)/2;
  height=(ye-ys)/2;
}


void getObjectCenter(DARY *mask,int &x, int &y){
  float fx=0,fy=0,sum=0;
  for(uint j=0;j<mask->y();j++){
    for(uint i=0;i<mask->x();i++){
      if(mask->bel[j][i]>0){
	fx+=i;fy+=j;sum++;
      }
    }
  }
  x=(int)(fx/sum);
  y=(int)(fy/sum);

}
  
int checkOverlap(DARY *mask, FD *cor){/*MASK is in bytes!!!*/
  
  int rad=(int)(0.5+cor->getScale());
  int x=(int)cor->getX();
  int y=(int)cor->getY();
  float fig=0,score;
  //cout << x << " "<< y << " "<< cor->getScale() << " " << rad << endl;getchar();
  for(int j=y-rad;j<=y+rad;j++){
    for(int i=x-rad;i<=x+rad;i++){
      if(i>=0 && j>=0 && i<(int)mask->x() && j<(int)mask->y()){
	if(mask->bel[j][i]>0)
	  fig++;
      }
    }
  } 
  score=fig/square(2*rad+1);
  cor->setArea(score);
  if(score>0.3)return 1;
  else return 0;        
}

void selectWithMask(DARY *mask, vector<FD *> &desc, vector<pair<int,int> > &dims){
  int x,y;
  getObjectCenter(mask,x, y);
  int width,height;
  getBoundingBox(mask,x, y, width, height);
  dims.push_back(make_pair(width,height));
  for(uint i=0;i<desc.size();i++){
    if(checkOverlap(mask, desc[i])){
      desc[i]->setX_Y(x-desc[i]->getX(),y-desc[i]->getY());
    }else{
      delete desc[i];
      desc.erase((std::vector<FD*>::iterator)&desc[i]);
      i--; 
    }    
  }  
  delete mask; 
}   

void detectFeatures(const char *name, DARY *img, DARY *mask, vector<FD *> &features,  vector<pair<int,int> > &dims, int detector, int DECSRIPTOR, float image_scale, int image_shift, float det_threshold, int max_desc = 0){
  
  vector<FD *> desc;
  if(image_shift==0 && image_scale!=1){
    extractFeatures(img, desc, detector, DECSRIPTOR, 0, 0, det_threshold, max_desc);
  }else{
    DARY *img2 = new DARY((int)(img->y()*image_scale),(int)(img->x()*image_scale+2*image_shift));
    img2->scale(img,1.0/image_scale,1.0/image_scale);
    extractFeatures(img2, desc, detector, DECSRIPTOR, 0, 0, det_threshold, max_desc);
    //img2->writePNG("test.png");getchar();
    delete img2;
  }
  for(uint d=0;d<desc.size();d++){ 
    desc[d]->setImageName(name);         
  }
  cout << "tot nb "<< desc.size()<< endl; 
  if(mask!=NULL)selectWithMask(mask,desc, dims);
  //displayFeatures(img, desc, "feature.png", 255,1);getchar();//remove object center
  delete img;
  //cout << desc.size()<< endl;
  features.insert(features.begin(),desc.begin(),desc.end());
  //cout << "done"<< endl;getchar();
  desc.clear();
}


void detectFeatures(const char *filein, vector<FD *> &features, int mask, vector<pair<int,int> > &dims, int detector, int DECSRIPTOR, float image_scale, int image_shift, float det_threshold, int max_desc = 0){

  vector<FD *> desc;
  DARY *img = new DARY(filein); 
  img->toGRAY();
  img->char2float();
  DARY *imask=NULL;
  if(mask){
    char tname[512];
    char name[512];

    strcpy(tname,filein);
    char *dir=strrchr(tname,'/');
    if(!dir)dir=tname;
    else dir++;
    strcpy(name,dir);
    char *suf=strrchr(name,'.');
    sprintf(suf,"-map.png");
    sprintf(dir,"maps/%s",name);
    //cout <<imgname<< " "<<  tname<< endl;
    
    imask = new DARY(tname);
    imask->toGRAY();
  }
  detectFeatures(filein, img, imask, features, dims, detector, DECSRIPTOR, image_scale, image_shift, det_threshold);
}
