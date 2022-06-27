#include "feature.h"
#include "../matrix/matrix.h"

void displayCircle(DARY *im, FeatureDescriptor *cor, float colori){
  int x=(int)cor->getX();
  int y=(int)cor->getY();
  //cout << x <<" display " << y << endl;
  for(int i=-1;i<=1;i++){
    if(y>0 && y<im->y() && y+i>0 && y+i<im->y() && x>0 && x<im->x() && x+i>0 && x+i< im->x())
      im->fel[y][x+i]=im->fel[y+i][x]=colori;
  }
    return;
  float sizef=1*(cor->getScale());//(patch_size/2)*scale/(2*step)
  float thick=1;
  //cout << "OK1 " << x << "  " << y << "  " << sizef << endl;
  float color=0;
  for(int nb=0;nb<2;nb++){ 
  for(float sizefi = sizef-thick;sizefi<=sizef+thick;sizefi+=0.2){
    
    uint size=(uint)rint(sizefi);
    for(uint i=0;i<=size;i++){
      int yi=(int)rint(sqrt((sizefi*sizefi-i*i)));
      if(x+i>0 && x+i<im->x()){
	if(y+yi<im->y()&& y+yi>0)
	  im->fel[y+yi][x+i]=color;
	if(y-yi>0 && y-yi<im->y()) 
	  im->fel[y-yi][x+i]=color;
      }
      if(x-i>0 && x-i<im->x()){
	if(y+yi<im->y()&& y+yi>0)
	  im->fel[y+yi][x-i]=color;
	if(y-yi>0 && y-yi<im->y()) 
	  im->fel[y-yi][x-i]=color;
      }	  
      if(y+i<im->y() && y+i>0){
	if(x+yi<im->x()&& x+yi>0)
	  im->fel[y+i][x+yi]=color;
	if(x-yi>0 && x-yi<im->x())
	  im->fel[y+i][x-yi]=color;
      } 
      if(y-i<im->y() && y-i>0){
	if(x+yi<im->x()&& x+yi>0)
	  im->fel[y-i][x+yi]=color;
	if(x-yi>0 && x-yi<im->x())
	  im->fel[y-i][x-yi]=color;
      }
    } 
  }  
  thick=0;color=colori;
  }

  //cout << "OK2 " << endl;

  //imn->write("imn.pgm");cout << "wrote" << endl;getchar();

}

void displayEllipse(DARY *im, FeatureDescriptor *cor, float colori){
  
  //cout << "OKe "<< endl;
  Matrix M(2,2),U,V,D,rot;
  M(1,1)=cor->getMi11();
  M(1,2)=cor->getMi12();
  M(2,1)=cor->getMi21();
  M(2,2)=cor->getMi22();
  M.svd(U,D,V);
  rot=V;
  //cout << "OKe2"<< endl;
  float scalex = cor->getScale();
  //if(scalex<8)return;
  float scale=1;
  float color=0;
  //  for(float sf=0.98;sf<1.01;sf+=0.005){   
  float sf=1;//(patch_size/2)*scale/(2*step)//2 for - 2 level in harris lap
  float v=sf*scale*D(1,1)*scalex;
  float n=sf*scale*D(2,2)*scalex;
  //cout << scalex << " scale " << endl;getchar();
  uint x=(int)cor->getX();
  uint y=(int)cor->getY();

  for(int i=-1;i<=1;i++){
    if(y>0 && y<im->y() && y+i>0 && y+i<im->y() && x>0 && x<im->x() && x+i>0 && x+i< im->x())
      im->fel[y][x+i]=im->fel[y+i][x]=colori;
  }
  //return;
  float thick=0;
  for(int nb=0;nb<2;nb++){
    for(float a=v-thick;a<=v+thick;a+=0.4)  
      for(float b=n-thick;b<=n+thick;b+=0.4){  
      
  //cout << a << " " << b << endl;
  float b2=b*b,a2=a*a;
  float b2_a2=b2/a2;
  float a2_b2=a2/b2;
  int xe = (int)cor->getX();
  int ye = (int)cor->getY();
  
  float xt, yt;
  int px,py;
  
  int rnb=im->y();
  int cnb=im->x();

  for(float i=0;i<a;i++){               
    yt=sqrt((float)(b2-b2_a2*i*i));
    px=xe+(int)rint(rot(1,1)*i+rot(1,2)*yt);
    py=ye+(int)rint(rot(2,1)*i+rot(2,2)*yt);
    if(px>=0 && px<cnb && py>=0 && py<rnb)
      im->fel[py][px]=color;
    
    px=xe+(int)rint(rot(1,1)*i-rot(1,2)*yt);
    py=ye+(int)rint(rot(2,1)*i-rot(2,2)*yt);
    if(px>=0 && px<cnb && py>=0 && py<rnb)
      im->fel[py][px]=color;
    
    px=xe+(int)rint(-rot(1,1)*i+rot(1,2)*yt);
    py=ye+(int)rint(-rot(2,1)*i+rot(2,2)*yt);
    if(px>=0 && px<cnb && py>=0 && py<rnb)
      im->fel[py][px]=color;
    
    px=xe+(int)rint(-rot(1,1)*i-rot(1,2)*yt);
    py=ye+(int)rint(-rot(2,1)*i-rot(2,2)*yt);
    if(px>=0 && px<cnb && py>=0 && py<rnb)
      im->fel[py][px]=color;
  }
  for(float j=0;j<b;j++){               
    xt=(int)rint(sqrt((float)(a2-a2_b2*j*j)));
    px=xe+(int)rint(rot(1,1)*xt+rot(1,2)*j);
    py=ye+(int)rint(rot(2,1)*xt+rot(2,2)*j);
    if(px>=0 && px<cnb && py>=0 && py<rnb)
      im->fel[py][px]=color;
    
    px=xe+(int)rint(rot(1,1)*xt-rot(1,2)*j);
    py=ye+(int)rint(rot(2,1)*xt-rot(2,2)*j);
    if(px>=0 && px<cnb && py>=0 && py<rnb)
      im->fel[py][px]=color;
    
    px=xe+(int)rint(-rot(1,1)*xt+rot(1,2)*j);
    py=ye+(int)rint(-rot(2,1)*xt+rot(2,2)*j);
    if(px>=0 && px<cnb && py>=0 && py<rnb)
      im->fel[py][px]=color;
    
    px=xe+(int)rint(-rot(1,1)*xt-rot(1,2)*j);
    py=ye+(int)rint(-rot(2,1)*xt-rot(2,2)*j);
    if(px>=0 && px<cnb && py>=0 && py<rnb)
      im->fel[py][px]=color;
  }
  }
  thick=0;
  color=colori;
  }
}

void displayEllipse(DARY *im, int xe, int ye, float scalex, float m11, float m12, float m21, float m22, float color){
  
  Matrix M(2,2),U,V,D,rot;
  M(1,1)=m11;
  M(1,2)=m12;
  M(2,1)=m21;
  M(2,2)=m22;
  M.svd(U,D,V);
  rot=V;
  // <<" DDD "<< D << endl;
  // cout << " l1 "<< ((m11+m22)+sqrt((m11+m22)*(m11+m22)-4*(m11*m22-m12*m21)))/2.0<< endl;
  //cout << " l2 "<< ((m11+m22)-sqrt((m11+m22)*(m11+m22)-4*(m11*m22-m12*m21)))/2.0<< endl;

  float a,b;
  // a=1.0/sqrt(D(1,1));b=1.0/sqrt(D(2,2));
  a=1*D(1,1)*scalex;b=1*D(2,2)*scalex;
  // cout << a << " " << b << endl;
  float b2=b*b,a2=a*a;
  float b2_a2=b2/a2;
  float a2_b2=a2/b2;
  
  float xt, yt;
  int px,py;
  int rnb=im->y();
  int cnb=im->x();

  //if(xe-2>=0 && xe+2<cnb && ye-2>=0 && ye+2<im->y()){
  for(int j=xe-2;j<=xe+2;j++)
    for(int i=0;i<=0;i++){
      if(j>=0 && j<cnb && ye+i>=0 && ye+i<rnb)
	im->fel[ye+i][j]=color;
    }  
  for(int j=ye-2;j<=ye+2;j++)
    for(int i=0;i<=0;i++){
      if(xe+i>=0 && xe+i<cnb && j>=0 && j<rnb)
	im->fel[j][xe+i]=color;
    }
  //}

  for(float i=0;i<a;i++){               
    yt=sqrt((float)(b2-b2_a2*i*i));
    px=xe+(int)rint(rot(1,1)*i+rot(1,2)*yt);
    py=ye+(int)rint(rot(2,1)*i+rot(2,2)*yt);
    if(px>=0 && px<cnb && py>=0 && py<rnb)
      im->fel[py][px]=color;
    
    px=xe+(int)rint(rot(1,1)*i-rot(1,2)*yt);
    py=ye+(int)rint(rot(2,1)*i-rot(2,2)*yt);
    if(px>=0 && px<cnb && py>=0 && py<rnb)
      im->fel[py][px]=color;
    
    px=xe+(int)rint(-rot(1,1)*i+rot(1,2)*yt);
    py=ye+(int)rint(-rot(2,1)*i+rot(2,2)*yt);
    if(px>=0 && px<cnb && py>=0 && py<rnb)
      im->fel[py][px]=color;
    
    px=xe+(int)rint(-rot(1,1)*i-rot(1,2)*yt);
    py=ye+(int)rint(-rot(2,1)*i-rot(2,2)*yt);
    if(px>=0 && px<cnb && py>=0 && py<rnb)
      im->fel[py][px]=color;
  }
  for(float j=1;j<b;j++){               
    xt=(int)rint(sqrt((float)(a2-a2_b2*j*j)));
    px=xe+(int)rint(rot(1,1)*xt+rot(1,2)*j);
    py=ye+(int)rint(rot(2,1)*xt+rot(2,2)*j);
    if(px>=0 && px<cnb && py>=0 && py<rnb)
      im->fel[py][px]=color;
    
    px=xe+(int)rint(rot(1,1)*xt-rot(1,2)*j);
    py=ye+(int)rint(rot(2,1)*xt-rot(2,2)*j);
    if(px>=0 && px<cnb && py>=0 && py<rnb)
      im->fel[py][px]=color;
    
    px=xe+(int)rint(-rot(1,1)*xt+rot(1,2)*j);
    py=ye+(int)rint(-rot(2,1)*xt+rot(2,2)*j);
    if(px>=0 && px<cnb && py>=0 && py<rnb)
      im->fel[py][px]=color;
    
    px=xe+(int)rint(-rot(1,1)*xt-rot(1,2)*j);
    py=ye+(int)rint(-rot(2,1)*xt-rot(2,2)*j);
    if(px>=0 && px<cnb && py>=0 && py<rnb)
      im->fel[py][px]=color;
  }
}




void displayEllipse(DARY *im, Matrix *M, int x, int y, float scalex){
  Matrix Mi;
  Mi =(*M);// M->inverse(); 
  int  size=(int)(3*rint(scalex)),xi,yi;
  float color=255;
  Matrix X(2,1),Xe;
  for(int i=-2;i<=2;i++){
    im->fel[y+i][x]=im->fel[y][x+i]=255;
  }
  int rnb=im->y();
  int cnb=im->x();

  for(int i=-size;i<=size;i++){
      X(1,1)=i;
      X(2,1)=(int)rint(sqrt((float)(size*size-i*i)));
      Xe=Mi*X;
      xi=(int)rint(Xe(1,1));
      yi=(int)rint(Xe(2,1));
      if(y+yi>0 && y+yi<rnb && x+xi>0 && x+xi<cnb)
	im->fel[y+yi][x+xi]=color;
      if(y-yi>0 && y-yi<rnb&& x-xi>0 && x-xi<cnb)
	im->fel[y-yi][x-xi]=color;
      X(1,1)=X(2,1);
      X(2,1)=i;
      Xe=Mi*X;
      xi=(int)rint(Xe(1,1));
      yi=(int)rint(Xe(2,1));
      if(y+yi>0 && y+yi<rnb && x+xi>0 && x+xi<cnb)
	im->fel[y+yi][x+xi]=color;
      if( y-yi>0 && y-yi<rnb&& x-xi>0 && x-xi<cnb)
	im->fel[y-yi][x-xi]=color;
  }
}

void displayFeatures(const char *filein, vector<FeatureDescriptor*> cor, char *filename, float color, int circle){
    DARY *img = new DARY(filein);
    img->toGRAY();
    img->char2float();
    displayFeatures(img, cor, filename, color, circle );
} 


void displayFeatures(DARY *image, vector<FeatureDescriptor*> cor, char *filename, float color, int circle){
    DARY *im = new DARY(image);
    color=400;
    for(uint i=0;i<cor.size();i++){   
	  if(circle)displayCircle(im, cor[i], color);
	  else  displayEllipse(im, cor[i], color);
    } 
    //im->write(filename);
    
    DARY *cim = new DARY(im->y(),im->x(),"3uchar");

    for(uint j=0;j<cim->y();j++){
      for(uint i=0;i<cim->x();i++){
	if(im->fel[j][i]==400){
	  cim->belr[j][i]=255;
	  cim->belg[j][i]=255;
	  cim->belb[j][i]=0;
	}else{
	  cim->belr[j][i]=(uchar)im->fel[j][i];
	  cim->belg[j][i]=(uchar)im->fel[j][i];
	  cim->belb[j][i]=(uchar)im->fel[j][i];
 
	} 
      }
    }    

    cim->writePNG(filename);
    delete im;delete cim;
}
void displayCircle(DARY *im,int x, int y, float scalex){
  int  size=(int)(rint(scalex)),yi;
  float color=255;
  for(int i=-2;i<=2;i++)im->fel[y+i][x]=im->fel[y][x+i]=color;
  for(int i=-size;i<=size;i++){
    yi=(int)rint(sqrt((float)(size*size-i*i)));	
    im->fel[y+yi][x+i]=im->fel[y-yi][x+i]=color;
    im->fel[y+i][x+yi]=im->fel[y+i][x-yi]=color;  
  }
}
