
#include "feature.h"
#include "../matrix/matrix.h"
#include "../util/util.h"
#include "../gauss_iir/gauss_iir.h"
#define ORI_THRESHOLD 0.80

void computeHistAngle(DARY *grad_im, DARY *angle_im,vector<float> &angles);

FeatureDescriptor::~FeatureDescriptor(){
  if(size>0)delete[] vec;
  if(cluster_dist!=NULL)delete cluster_dist;
  if(imagename!=NULL)delete[] imagename; 
  for(uint i=0;i<occurences.size();i++)delete occurences[i];
  occurences.clear();
}
 
FeatureDescriptor::FeatureDescriptor(float xin, float yin, float scale_in, float featureness_in){
  init();
  x=xin;y=yin;c_scale=scale_in;featureness=featureness_in;
}

void deleteDescriptors(vector<FeatureDescriptor *> &desc){
  for(uint i=0;i<desc.size();i++){
    if(desc[i]!=NULL){
      if(desc[i]->features.size()>0)deleteDescriptors(desc[i]->features);
      delete desc[i];
    }
  }
  desc.clear(); 
}
void FeatureDescriptor::Cout(int dim){
  cout << "x " << x << " y " << y << " s " << c_scale<< " nbf " << nbf << " var " << var;
  if(size>0){
    cout<< " vec "; 
    for(int i=0;i<dim;i++)cout<< vec[i]<< " ";
  }
  cout << endl;
}

void  FeatureDescriptor::findNearestNeighbor(vector<FD*> features, int dim, int &ind, float &sim){
  sim=MAX_FLOAT;
  float d;
  ind=-1;
  for(uint i=0;i<features.size();i++){
    if(features[i]!=NULL){
      d=var+features[i]->getVar()+distEuc(vec,features[i]->getVec(),dim);
      if(d<sim){
	sim=d;
	ind=i;
      }      
    }
  }  
}


/**********************************************/
void FeatureDescriptor::copy(FeatureDescriptor* ds){
  x=ds->getX(); 
  y=ds->getY();
  type=ds->getType();
  featureness=ds->getFeatureness();
  angle=ds->getAngle();
  eangle=ds->getEAngle();
  c_scale=ds->getScale();
  der_sig=ds->getDerSig();
  int_sig=ds->getIntSig();
  extr=ds->getExtr();
  l1=ds->getL1();
  l2=ds->getL2();
  mi11=ds->getMi11();
  mi12=ds->getMi12();
  mi21=ds->getMi21();
  mi22=ds->getMi22();
  nbf=ds->getNbf();
  var=ds->getVar();
  obj=ds->getObj();
  weight=ds->getWeight();
  radius=ds->radius;
  lap=ds->getLap();
  area=ds->getArea();
  tree_lev=ds->getTreeLevel();
  size=ds->getSize();
  imagename=NULL;
  cluster_dist=NULL; 
  if(size>0){
    allocVec(ds->getSize());
    float *vec2=ds->getVec();
    memcpy(vec, vec2, sizeof(float)*size);
  }
  //for(int i=0;i<size;i++)vec[i]=ds->getV(i);
}

void FeatureDescriptor::init(){
	x=0; 
	y=0;
	type=0;
	obj=0;
	lap=0;
	featureness=1;
	angle=1000;
	c_scale=1;
	der_sig=0;
	int_sig=0;
	extr=0;
	l1=1;
	l2=0;
	eangle=0;
	mi11=1;
	mi12=0;
	mi21=0;     
	mi22=1;
	nbf=1;
	radius=0;
	weight=0;
	tree_lev=0;
	size=0;
	var=0;
	area=0;	
	mean_dist=0;
	vec=NULL;
	imagename=NULL;
	cluster_dist=NULL; 
}

/**********************************************/
void FeatureDescriptor::allocVec(int size_in){
    if(size_in >0){
	if(vec!=NULL)delete [] vec;
	size=size_in;
	vec = new float[size];
	for(int i=0;i<size;i++)vec[i]=0;
    }
}
  
/**********************************************/
void  FeatureDescriptor::setImageName(const char* name) { 
  if(imagename==NULL)imagename=new char[512];
  strcpy(imagename,name);
}


/**********************************************/
void FeatureDescriptor::read( ifstream &input, int size_in){
  input >> x;
  input >> y; 
  input >> featureness; 
  input >> c_scale;
  input >> angle;
  input >> obj;
  input >> type;
  input >> lap;
  input >> extr;
  input >> mi11;
  input >> mi12;
  input >> mi21;
  input >> mi22;
  d_scale=c_scale;
  if(size_in>0){
      //input >> d_scale;
    allocVec(size_in);
    for(int j=0;j<size;j++){
      input >> vec[j];
    }
  } 
} 

/**********************************************/
void FeatureDescriptor::readCommon( ifstream &input, int size_in){
  double a,b,c;
  Matrix U(2,2,0.0),D,Vi,V;
  input >> x;
  if(input.eof()) return; // will be deleted
  input >> y;
  input >> a;
  input >> b;
  input >> c;
  U(1,1)=a;
  U(1,2)=b;
  U(2,1)=b;
  U(2,2)=c;
  U.svd(Vi,D,V);
  D(1,1)=(1.0/sqrt(D(1,1)));
  D(2,2)=(1.0/sqrt(D(2,2)));
  a=sqrt(D(2,2)*D(1,1));
  //cout << D(1,1)<< " " <<  D(2,2)<< "  " << tmp1<< endl;
  D.tabMat[2][2]/=a;
  D.tabMat[1][1]/=a;
  U=V*D*V.transpose();
  c_scale=a;
  //cout << x << " " << y << " " << a << " "<< b << " "<< c << " " << size_in<<endl;getchar();
  mi11=U(1,1);
  mi12=U(1,2);
  mi21=U(2,1);
  mi22=U(2,2);
  
  if(size_in>0){
    allocVec(size_in);
    for(int j=0;j<size;j++){
      input >> vec[j];
    }
  } 
} 



/**********************************************/
void FeatureDescriptor::write(ofstream &output){
  output << x << " " << y << " " << featureness << " " << c_scale << " ";
  output << angle << " " << obj << " "<< type << " " << lap << " " <<extr;
  output << " " << mi11 << " " << mi12 << " " << mi21 << " " << mi22<< " ";
  if(size>0){
      //output  << "  " << d_scale << " ";
    for(int j=0; j<size;j++){
      output << vec[j] << " ";
    }
  }
  output << endl;
}
/**********************************************/
void FeatureDescriptor::writeBin(FILE *fid, float *buf){
  int pos=0;
  buf[pos++]=x;//1
  buf[pos++]=y; //2
  buf[pos++]=featureness; 
  buf[pos++]=c_scale;
  buf[pos++]=angle;
  buf[pos++]=area;
  buf[pos++]=type;
  buf[pos++]=lap;
  buf[pos++]=extr;
  buf[pos++]=weight;
  buf[pos++]=obj;
  buf[pos++]=tree_lev;
  buf[pos++]=radius;
  buf[pos++]=mi11;
  buf[pos++]=mi12;
  buf[pos++]=mi21;
  buf[pos++]=mi22;
  buf[pos++]=l2;
  buf[pos++]=eangle;//19
  buf[pos++]=mean_dist;//20
  for(int v=0;v<size;v++)
    buf[pos++]=vec[v];
  fwrite(buf, sizeof(float),size+20, fid);
}

 /**********************************************/
void FeatureDescriptor::readBin(FILE *fid, int size, float *buf){

  allocVec(size);
  fread(buf, sizeof(float), size+20 , fid);
  
  int pos=0;
  x=buf[pos++];
  y=buf[pos++]; 
  featureness=buf[pos++]; 
  c_scale=buf[pos++];
  angle=buf[pos++];
  area=buf[pos++];
  type=(int)buf[pos++];
  lap=buf[pos++];
  extr=(int)buf[pos++];
  weight=buf[pos++];
  obj=(int)buf[pos++];
  tree_lev=(int)buf[pos++];
  radius=buf[pos++];
  mi11=buf[pos++];
  mi12=buf[pos++];
  mi21=buf[pos++];
  mi22=buf[pos++];
  l2=buf[pos++];
  eangle=buf[pos++];
  mean_dist=buf[pos++];
  for(int v=0;v<size;v++)
    vec[v]=buf[pos++];
}


/**********************************************/
void FeatureDescriptor::writeCommon(ofstream &output){

    Matrix U(2,2,0.0),D,Vi,V;
    U(1,1)=mi11;
    U(1,2)=mi12;
    U(2,1)=mi21;
    U(2,2)=mi22;
    U.svd(Vi,D,V);
    D=D*(c_scale);
    D(1,1)=1.0/(D(1,1)*D(1,1));
    D(2,2)=1.0/(D(2,2)*D(2,2));
    U=V*D*V.transpose();
    
    output << x << " " << y << " " << U(1,1)<< " "<< U(1,2)<< " "<< U(2,2);
    if(size>0){
      for(int j=0; j<size;j++){
	output  << " " << vec[j];
      }
    }
    output << endl;
} 

/**********************************************/
void FeatureDescriptor::changeBase(float *mat){
  for(int v=0;v<size;v++){
    vec[v] = vec[v]*mat[v];   
  }
}

/**********************************************/
void FeatureDescriptor::pca(int dim, float *avg, float *base){
  float *outvec = new float[dim];
  for(int v=0;v<dim;v++)outvec[v]=0;
  
  for(int v=0;v<size;v++)vec[v]-=avg[v];  
  uint cnt=0;
  for(int i=0;i<dim;i++){
    for(int v=0;v<size;v++,cnt++){
      outvec[i] += vec[v]*base[cnt];   
    }
    //cout << outvec[i]<< endl;
  }
  for(int i=0;i<dim;i++){
    vec[i]=outvec[i];
  } 
  //getchar();
  delete []outvec;
  size=dim;
  //vec = outvec;
}

/**********************************************/
void loadFeatures( const char* points1, vector<FeatureDescriptor*> &cor1, int format){
    FeatureDescriptor* cor;
    ifstream input1(points1);
    if(!input1)return;
    int cor_nb1,size; 
    float tmp;
    input1 >> tmp;
    input1 >> cor_nb1;   //cout << cor_nb1 << endl;
    if(tmp<=1.0)size=0;
    else size=(int)tmp;
    if(cor_nb1==0)return;
    while(input1.good()){
      cor = new FeatureDescriptor();
      if(format==0)
	cor->read(input1,size);
      else
	cor->readCommon(input1,size);
      cor1.push_back(cor); //cout << "read " <<cor1.size()<<endl;     
    }
    cor1.erase((std::vector<FeatureDescriptor*>::iterator) &cor1[(int)cor1.size()-1]);
    if(cor_nb1!=(int)cor1.size()){
      cout << "warning:"<< endl<<"in header: "<<cor_nb1 << ", in file: "<< cor1.size()<< endl; 
    }

}

/**********************************************/
void writeFeatures(vector<FeatureDescriptor*> cor, const char* points_out, int format){
  if (cor.size()==0){
      cout << " descriptors nb " << cor.size() << endl;
      return; 
  }
    ofstream output(points_out);
    if(!output)cout << "error opening " << points_out<< endl;

    output << cor[0]->getSize()<< endl;
    output << cor.size()<< endl;
    for(unsigned int i=0; i<cor.size();i++){
      if(format==0)
	cor[i]->write(output);
      else  cor[i]->writeCommon(output);
    }  
    output.close();  
} 


/*struct CmpCor {
  bool operator()(FeatureDescriptor*f1,FeatureDescriptor*f2) {
    return f1->featureness > f2->featureness;
  }
};
*/
void writeFeaturesHerve(vector<FeatureDescriptor*> cor, const char * points_out) {
  /*  sort(cor.begin(),cor.end(),CmpCor());
  int nmax = 600;
  int n = min(nmax, int(cor.size()));*/
  int n = int(cor.size());

  FILE * fo = fopen(points_out,"wb");
  if(!fo) fprintf(stderr, "big error :cannot open file");

  for(int i = 0; i < n; i++) { 
    FeatureDescriptor* fd = cor[i];
    int j;
    int d = fd->getSize();
    double norm = 0.0;
    for (j = 0 ; j < d ; j++)
      norm += ( fd->vec[j] ) * ( fd->vec[j] ) ;
    norm = sqrt( norm ) ;
    
    for (j = 0 ; j < d ; j++)
      fd->vec[j] /= norm;
    
    fwrite (&d, sizeof (int), 1, fo);
    fwrite (fd->vec, sizeof (float), d, fo);
  }

  fclose(fo);

}

void writeCoordinates(vector<FeatureDescriptor*> cor, const char * coord_out) {
	
  int n = int(cor.size());

  FILE * fo = fopen(coord_out,"wb");
  if(!fo) fprintf(stderr, "big error :cannot open file for coordinates");

  for(int i = 0; i < n; i++) { 
    FeatureDescriptor* fd = cor[i];
       
    fprintf(fo,"%g %g %g   %g %g %g %g   %g\n",
            fd->x,fd->y,fd->c_scale,
            fd->mi11,fd->mi12,fd->mi21,fd->mi22,
            fd->angle);
  }

  fclose(fo);
 
}

void writeFeaturesSameh(const vector<FeatureDescriptor*> & cor, const char * points_out) {
  int n = int(cor.size());
  int looseprec=0;
  FILE * fo = fopen(points_out,"wb");
  if(!fo) fprintf(stderr, "big error :cannot open file for coordinates");

  for(int i = 0; i < n; i++) { 
    FeatureDescriptor* fd = cor[i];
    int m=fd->getSize();
    unsigned char buf[m];
    
    struct {
      float x, y, scale, angle, mi11, mi12, mi21, mi22, cornerness;
      int dim;
    } geom={
      fd->x,fd->y,fd->c_scale,fd->angle,
      fd->mi11,fd->mi12,fd->mi21,fd->mi22,
      fd->featureness,m
    };
    
    fwrite(&geom,sizeof(geom),1,fo);
    
    for(int i=0;i<m;i++) {
      float v=fd->vec[i];
      if(v!=floor(v)) 
        looseprec++;
      buf[i]=(int)v;
    }
    fwrite(buf,1,m,fo);    
    
  }

  if(looseprec) 
    fprintf(stderr,
            "warning: rounding in %s lost precision in %d components!\n",
            points_out,looseprec);

  fclose(fo);
}




/**********************************************/
void saveBinFeatures(vector<FeatureDescriptor*> cor, const char* points_out){
  cout << "Saving " << cor.size() << " features in " << points_out << "... "<< flush;
  if (cor.size()==0){
    cout << "error saving features  " << endl;
    return; 
  }
  FILE *fid;
  fid=fopen(points_out, "wb");

  if(!fid)cout << "error opening " << points_out<< endl;

  int size=cor[0]->getSize();
  fwrite(&size, sizeof(int), 1, fid);
  uint nb=cor.size();
  fwrite(&nb, sizeof(uint), 1, fid);
  
  int pos=0;
  uint totsize=(size+16)*nb;
  float *write_vec = new float[totsize];
  for(unsigned int i=0; i<cor.size();i++){
    write_vec[pos++]=cor[i]->getX();
    write_vec[pos++]=cor[i]->getY(); 
    write_vec[pos++]=cor[i]->getFeatureness(); 
    write_vec[pos++]=cor[i]->getScale();
    write_vec[pos++]=cor[i]->getAngle();
    write_vec[pos++]=cor[i]->getObj();
    write_vec[pos++]=cor[i]->getArea();
    write_vec[pos++]=cor[i]->getType();
    write_vec[pos++]=cor[i]->getLap();
    write_vec[pos++]=cor[i]->getExtr();
    write_vec[pos++]=cor[i]->getMi11();
    write_vec[pos++]=cor[i]->getMi12();
    write_vec[pos++]=cor[i]->getMi21();
    write_vec[pos++]=cor[i]->getMi22();
    write_vec[pos++]=cor[i]->getL2();
    write_vec[pos++]=cor[i]->getEAngle();
    for(int v=0;v<size;v++)
      write_vec[pos++]=cor[i]->getV(v);
  }  
  fwrite(write_vec, sizeof(float),totsize , fid);
  
  fclose(fid); 
  delete [] write_vec;
  cout << "done"<< endl;
} 
/**********************************************/
void loadBinFeatures(const char* points_out, vector<FeatureDescriptor*> &cor){
  cout << "Loading features from " << points_out << "... "<<  flush;
    
  FILE *fid;
  fid=fopen(points_out, "rb");

  if(!fid)cout << "error opening " << points_out<< endl;

  int size;  

  fread(&size, sizeof(int), 1, fid);
  uint nb;
  fread(&nb, sizeof(uint), 1, fid);
  
  int pos=0;
  uint totsize=(size+16)*nb;
  float *write_vec = new float[totsize];

  fread(write_vec, sizeof(float), totsize , fid);
  float *vec;
  for(unsigned int i=0; i<nb;i++){
    cor.push_back(new FeatureDescriptor());
    cor[i]->setX_Y(write_vec[pos],write_vec[pos+1]);pos+=2;
    cor[i]->setFeatureness(write_vec[pos++]);
    cor[i]->setScale(write_vec[pos++]);
    cor[i]->setAngle(write_vec[pos++]);
    cor[i]->setObj((int)write_vec[pos++]);
    cor[i]->setArea(write_vec[pos++]);
    cor[i]->setType((int)write_vec[pos++]);
    cor[i]->setLap(write_vec[pos++]);
    cor[i]->setExtr((int)write_vec[pos++]);
    cor[i]->setMi(write_vec[pos],write_vec[pos+1],write_vec[pos+2],write_vec[pos+3]);
    pos+=4;
    cor[i]->setL2(write_vec[pos++]);
    cor[i]->setEAngle(write_vec[pos++]);
    cor[i]->allocVec(size);

    vec=cor[i]->getVec();
    for(int v=0;v<size;v++)
      vec[v]=write_vec[pos++];
  }  
  
  fclose(fid); 
  delete [] write_vec;
  cout <<cor.size()<<  " done"<< endl;

} 


/************************************************************************/
/************************************************************************/
/************************************************************************/
/************************************************************************/




/**********************************************/
void smoothHistogram(float *hist, int bins)
{
	int i;
	float prev, temp;
   
	prev = hist[bins - 1];
	for (i = 0; i < bins; i++) {
		temp = hist[i];
		hist[i] = (prev + hist[i] + hist[(i + 1 == bins) ? 0 : i + 1]) / 3.0;
		prev = temp;
	}
}



/******************Lowe method*********************/
void computeHistAngle(DARY *grad, DARY *ori,vector<float> &angles){
	if(ori==NULL ||grad==NULL ){return;}
	int i, r, c, rows, cols, bin, prev, next,OriBins=36;
	float hist[OriBins],  gval, langle,dbin,fbin,  interp, maxval = 0.0;/*OriBins=36*/
   
	rows = grad->y();
	cols = grad->x();
   
	for (i = 0; i < OriBins; i++)
		hist[i] = 0.0;

	for (r = 1; r < rows - 1; r++)
	  for (c = 1; c <= cols - 1; c++){
	    gval = grad->fel[r][c];   
	    if (gval > 1.0  &&  patch_mask->fel[r][c]>0) {
	      /* Ori is in range of -PI to PI. */
	      langle = ori->fel[r][c] + M_PI;
	      fbin =  (OriBins * langle / (M_2PI));
	      bin=(int)floor(fbin);
	      dbin=fbin-bin;
	      //assert(bin >= 0 && bin <= OriBins);
	      bin = (bin < OriBins)?bin:(0);
	      hist[bin] +=  (1-dbin)*gval * patch_mask->fel[r][c];
	      bin = (bin+1 < OriBins)?bin+1:(0);
	      hist[bin] +=  (dbin) * gval * patch_mask->fel[r][c];
	    } 
	  }
	
	/* Apply smoothing 6 times for accurate Gaussian approximation. */
	for (i = 0; i < 6; i++)
	  smoothHistogram(hist, OriBins);
	
	
	/* Find maximum value in histogram. */
	int maxi=-1; 
	for (i = 0; i < OriBins; i++){
	  if (hist[i] > maxval){
	    maxval = hist[i];
	    maxi=i;
	  }
	}
	//cout << endl;
	for (i = 0; i < OriBins; i++){
	  prev = (i == 0 ? OriBins - 1 : i - 1);
	  next = (i == OriBins - 1 ? 0 : i + 1);
	  //cout << i  << " "<< hist[i];
	  if (hist[i] > hist[prev] && hist[i] > hist[next] && hist[i]/maxval>ORI_THRESHOLD){
	    interp = interpPeak(hist[prev], hist[i], hist[next]);
	    angles.push_back(M_2PI * (i + 0.5 + interp) / OriBins - M_PI); 
	    //angles.push_back(M_2PI * (i + interp) / OriBins - M_PI); 
	    //cout <<i <<" " <<  angles[angles.size()-1];      
	  }
	  // cout << endl;
	}   
}

int normalize(DARY * img, int x, int y, float radius){
	float sum=0;
	float gsum=0; 

	for(uint j=0;j<img->y();j++){ 
		for(uint i=0;i<img->x();i++){ 
			if(patch_mask->fel[j][i]>0){
				sum+=img->fel[j][i]; 
				gsum++;
			}
		} 
	}    
	sum=sum/gsum;       
        
	float var=0;
	for(uint j=0;j<img->y();j++){ 
		for(uint i=0;i<img->x();i++){ 
			if(patch_mask->fel[j][i]>0){	
				var+=(sum-img->fel[j][i])*(sum-img->fel[j][i]);	
			}
		}
	}     
        if(var==0) 
          return 0;

	var=sqrt(var/gsum);            

    //  cout << "mean "<<sum<< " " <<img->fel[y][x] << " var " << var << endl;
	float fac=50.0/var;
	float max=0,min=1000;
	for(uint j=0;j<img->y();j++){ 
		for(uint i=0;i<img->x();i++){ 
			img->fel[j][i]=128+fac*(img->fel[j][i]-sum);
			if(max<img->fel[j][i])max=img->fel[j][i];
			if(min>img->fel[j][i])min=img->fel[j][i];
			if(img->fel[j][i]>255)img->fel[j][i]=255;
			if(img->fel[j][i]<0)img->fel[j][i]=0;
		}
	}   
    // cout << "max " << max << " min "<< min <<endl;

        return 1;
}



/************NORMALIZATION PATCH****************/

DARY *patch_mask = new DARY(PATCH_SIZE,PATCH_SIZE,1.0);
//float PATCH_SUM;
void initPatchMask(int size){ 
	int center=size>>1;
	float radius = center*center;
	float sigma=0.9*radius;
	float disq;
	for(int i=0;i<size;i++)
	  for(int j=0;j<size;j++){
	    disq=(i-center)*(i-center)+(j-center)*(j-center);
	    if(disq < radius){
	      patch_mask->fel[j][i]= exp(- disq / sigma);
	      //mask->fel[j][i]= 255*exp(- disq / sigma);   
	      //cout << patch_mask->fel[j][i]<< endl; 
	      //PATCH_SUM+=patch_mask->fel[j][i];
	    }else { 
	      patch_mask->fel[j][i]=0;
	    }		
	  } 
	
	//patch_mask->normalize(0,1);patch_mask->write("mask.pgm");cout << "mask "<< endl;getchar();
} 


void  compDesc(DARY *patch, FD *ds, int DESC){
  

    switch ( DESC )
    {
      case JLA:
	computeJLA(patch, ds);
        break;
      case CC:
 	computeCC(patch, ds);
        break;
      case CF:
 	computeCF(patch, ds);
        break;
      case KOEN:
 	computeKoen(patch, ds);
        break;
      case SPIN:
 	computeSpin(patch, ds);
        break;
      case SHAPE:
 	computeShape(patch, ds);
       break;
        case MOM:
	computeMoments(patch, ds);
        break;
      case PCA:
 	computePcaSift(patch, ds);
        break;
      case SIFT:
 	computeSift(patch, ds, 1);
        break;
      case SIFTNONORM:
 	computeSift(patch, ds, 0);
        break;
      case GLOH:
	computeESift(patch, ds);
        break;
    default:
        ;
    }


}


int  compDesc(DARY *dx, DARY *dy, DARY *grad, DARY *ori, FD *ds, int DESC){
  

    switch ( DESC )
    {
      case SHAPE:
 	computeShape(dx,dy,grad,ori, ds);
       break;
        case MOM:
	  //computeMoments(patch, ds);
        break;
      case PCA:
 	//computePcaSift(patch, ds);
        break;
      case SIFT:
 	return computeSift(dx,dy,grad,ori, ds, 1);
        break;
      case SIFTNONORM:
 	return computeSift(dx,dy,grad,ori, ds, 0);
        break;
      case GLOH:
	computeESift(dx,dy,grad,ori, ds);
        break;
    default:
      printf("compDesc(DARY *dx, DARY *dy, DARY *grad, DARY *ori, FD *ds, int DESC) is nop !\n");
        ;
    }
    return 1;
}


void computeDescriptor(DARY *dx, DARY *dy, int DESC,
			FeatureDescriptor *cor, vector<FeatureDescriptor *> &desc){
  
  vector<float> angles;
  DARY *gpatch = new DARY(dx->y(),dx->x());
  DARY *opatch = new DARY(dx->y(),dx->x());
  gradAngle(dx,dy,gpatch,opatch);  
  computeHistAngle(gpatch,opatch,angles);
  
  for(uint i=0;i<angles.size();i++){
    FeatureDescriptor *ds= new FeatureDescriptor();
    ds->copy(cor);
    ds->setAngle(angles[i]);
    compDesc(dx,dy,gpatch,opatch,ds,DESC);
    desc.push_back(ds);
    //{gpatch->writePNG("norm.png");cout << angles[i]<< " " << i << " angle  of " << angles.size() << endl;getchar();}
  }      
  delete gpatch;delete opatch;
  angles.clear();
}



/**********************************************/

/**********************COMPUTE DESCRIPTORS FROM EXTERNAL DETECTORS************************/
void computeAffineDescriptor( DARY *imgbs, DARY *patch, float scal, int DESC,
			FeatureDescriptor *cor, vector<FeatureDescriptor *> &desc){

  if(cor->getAngle()<M_PI && cor->getAngle()>-M_PI){  
    int patch_ok=normalize(patch,patch->x()>>1,patch->y()>>1,patch->x()>>1);
    if(!patch_ok) return;
    FeatureDescriptor *ds= new FeatureDescriptor();
    ds->copy(cor);
    compDesc(patch,ds, DESC);
    desc.push_back(ds);
    return;
  }   

  DARY * grad = new DARY(patch->x(),patch->x());   
  DARY * ori = new DARY(patch->x(),patch->x());   
  vector<float> angles;
  gradAngle(patch, grad, ori);
  computeHistAngle(grad,ori,angles);
  delete grad;delete ori;
  for(uint i=0;i<angles.size();i++){
    patch->interpolate(imgbs,imgbs->x()>>1,imgbs->y()>>1,1.0/scal,1.0/scal,-180*angles[i]/M_PI); // normalization done
    int patch_ok=normalize(patch,patch->x()>>1,patch->y()>>1,patch->x()>>1);
    if(!patch_ok) continue;    
    FeatureDescriptor *ds = new FeatureDescriptor();
    ds->copy(cor);
    ds->setAngle(angles[i]);
    compDesc(patch,ds, DESC);
    desc.push_back(ds);
    //if(scal>0){patch->writePNG("norm.png");cout << angles[i]<< " " << i << " angle  of " << angles.size() << endl;getchar();}
  }      
  angles.clear();
}

DARY* normalizeAffine(DARY *image, float x, float y, float c_scale, 
		      float angle, float mi11, float mi12,float mi21, float mi22, float &scal,
		      DARY *patch, float DESC_SCALE){
  
  // int imb_size=2*(int)(1.414*GAUSS_CUTOFF*c_scale)+1;//1.414 margin for rotation, CUTOFF for min gaussian kernel size 
  int imb_size=2*(int)(1.414*DESC_SCALE*c_scale)+1;//9 pixels margin for blur and rotation
  //  float scal=(2*(int)(GAUSS_CUTOFF*c_scale)+1)/((float)patch->x());//scale for sampling without margin
  scal=(2.0*DESC_SCALE*c_scale+1)/((float)patch->x());//scale for sampling without margin
  float lecos=1;
  float lesin=0; 

  if(angle<M_PI && angle>-M_PI){// if angle is already  estimated 
    lecos=cos(-angle);
    lesin=sin(-angle);
  } 
     
  float m11=(mi11*lecos-mi12*lesin);
  float m12=(mi11*lesin+mi12*lecos);
  float m21=(mi21*lecos-mi22*lesin);
  float m22=(mi21*lesin+mi22*lecos); 
  
  //  smooth before sampling
  //cout <<" x "<< x << " y "<< y << " cs "<< c_scale << " sc " << scal<< " imb " <<  imb_size<<" " <<  1.5*scal << " "<< PATCH_SIZE*scal<<  endl;
    //cout << angle << " m " <<mi11<< " " << mi12 << " " << mi21 << " " << mi22 << endl;  
  DARY * imgb=new DARY(imb_size,imb_size);
  imgb->interpolate(image,x,y,m11,m12,m21,m22);//imgb->writePNG("bsnorm.png");//normalize affine with scale=1 
  if(scal>1){
    DARY * imgbs=new DARY(imb_size,imb_size);
    smooth(imgb,imgbs,(scal)/1.3);//imgbs->writePNG("smnorm.png");//smooth 
    imgb->set(imgbs);
    //cout << c_scale << " " << 1.5*scal << endl; 
    delete imgbs;
  }
 
  patch->interpolate(imgb,imb_size>>1,imb_size>>1,1.0/scal,1.0/scal,0.0);//not rotated
  //patch->writePNG("patch.png"); getchar();
  return imgb;
}




/**********************COMPUTE DESCRIPTORS FROM EXTERNAL DETECTORS************************/
void computeAffineDescriptor( DARY *image, DARY *patch, int DESC, FeatureDescriptor *cor, vector<FeatureDescriptor *> &desc, float DESC_SCALE){
     
  float angle=cor->getAngle();
  float scal; 

  // DARY *imgbs= normalizeAffine(image, cor, patch, scal);
  //cor->Cout();
  DARY *imgbs = normalizeAffine(image, cor->getX(), cor->getY(), cor->getScale(), 
			       angle, cor->getMi11(), cor->getMi12(), cor->getMi21(),cor->getMi22(),
				scal,patch,DESC_SCALE);  
  computeAffineDescriptor(imgbs, patch, scal,DESC,cor,desc);
  delete imgbs;
  //estimate orientation angle
 
} 
 



/**********************COMPUTE DESCRIPTORS FROM EXTERNAL DETECTORS************************/
void computeAffineDescriptors(DARY *image,  vector<FeatureDescriptor *> &desc, int DESC, float DESC_SCALE){
  if(DESC!=MSIFT)initPatchMask(PATCH_SIZE);
  else DESC=SIFT;
    DARY * patch = new DARY(PATCH_SIZE,PATCH_SIZE);   
    vector<FeatureDescriptor *> tmpdesc;
    for(unsigned int c=0;c<desc.size();c++){
      computeAffineDescriptor( image, patch, DESC, desc[c], tmpdesc,DESC_SCALE);
	if(!(c%100))cout << "\rdescriptor "<< c<< " of "<< desc.size() << "    " << flush;
    }
    delete patch;
    desc.clear();
    desc=tmpdesc;
    cout << desc.size() << endl;
    //tmpdesc.clear();
}
 

void pca(vector<FD*> &features, float *mvec, float *base, int newdim){
  cout<< endl;
  for(uint i=0;i<features.size();i++){
    features[i]->pca(newdim,mvec,base);
    //for(int j=0;j<newdim;j++)cout << features[i]->getV(j)<< endl;
    //getchar();
  }
}

void pcaBase(vector<FD*> &features,float *mvec, float *base, int newdim){
  cout << "Estimating pca base for "<<  features.size()<< " features..." << endl;
  if(features.size()<100){
    cout << "Too few features for pca"<< endl;
    return; 
  }
  uint dim=features[0]->getSize();
  
  // compute mean  

  for(uint d=0;d<dim;d++){
    mvec[d]=0;
  }    
  float *vec;
  for(uint i=0;i<features.size();i++){
    vec=features[i]->getVec();
    for(uint d=0;d<dim;d++){
      mvec[d]+=vec[d];
    }    
  }
  float fnb=(float)features.size();
  for(uint d=0;d<dim;d++){
    mvec[d]=mvec[d]/fnb;
  }    
  for(uint i=0;i<features.size();i++){
    vec=features[i]->getVec();
    for(uint d=0;d<dim;d++){
      vec[d]-=mvec[d];
    }    
  }
  


  Matrix cov(dim,dim,0.0);
  float *f1;
  for(uint i=0;i<features.size();i++){
    f1=features[i]->getVec();
    for(uint d1=0;d1<dim;d1++){
      for(uint d2=d1;d2<dim;d2++){
	cov.tabMat[d1+1][d2+1]+=((f1[d1])*(f1[d2]));
      }
    }
  }
  
  for(uint d1=1;d1<=dim;d1++){
    for(uint d2=d1+1;d2<=dim;d2++){
      cov.tabMat[d2][d1]=cov.tabMat[d1][d2];
    }
  }

  
  Matrix V(dim,dim);
  Vector d(dim);
  //cov.write("cov.mat",1);
  cov.jacobi(d, V);
  //d.write("d.mat",1);
  //V.write("V.mat",1);

  int cnt=0;
  for(int d1=1;d1<=newdim;d1++){
    for(uint d2=1;d2<=dim;d2++){
      base[cnt++]=V.tabMat[d2][d1];
    }
  }
    

  float *outvec = new float[newdim];
  for(uint f=0;f<features.size();f++){
    vec=features[f]->getVec();
    for(int v=0;v<newdim;v++)outvec[v]=0;
    uint cnt=0;
    for(int i=0;i<newdim;i++){
      for(uint v=0;v<dim;v++,cnt++){
	outvec[i] += vec[v]*base[cnt];   
      }    
    }
    for(int i=0;i<newdim;i++){
      vec[i]=outvec[i];
    } 
    features[f]->setSize(newdim);
  }

  /*
  writeFeatures(features, "features.max.pca",1);

  Matrix cov2(newdim,newdim,0.0);
  for(uint i=0;i<features.size();i++){
    f1=features[i]->getVec();
    for(uint d1=0;d1<newdim;d1++){
      for(uint d2=d1;d2<newdim;d2++){
	cov2.tabMat[d1+1][d2+1]+=((f1[d1])*(f1[d2]));
      }
    }
  }
  
  for(uint d1=1;d1<=newdim;d1++){
    for(uint d2=d1+1;d2<=newdim;d2++){
      cov2.tabMat[d2][d1]=cov2.tabMat[d1][d2];
    }
  }

  cov2.write("cov2.mat",1);
  */


  delete []outvec;
}

int loadPCAVectors(const char *filename, float *&vmean, uint &mean_dim, float *&vbase, uint &base_dim){
  cout << "Loading pca vectors from " <<filename << endl;
  ifstream input(filename);
  if(!input){
    cout << "No pca vectors in " <<filename << endl;
    return 0;
  }
  input >> mean_dim;  
  vmean=new float[mean_dim];
  float val;
  for(uint i=0;i<mean_dim;i++){
    input >> val;
    vmean[i]=val;
  }
  input >> base_dim;  
  vbase=new float[base_dim];  
  for(uint i=0;i<base_dim;i++){
    input >> val;
    vbase[i]=val;
  }
  input.close();
  return 1;
}
void savePCAVectors(const char *filename, float *vmean, uint size1 ,float *vbase, uint size2){
  cout << "Saving pca vectors in " <<filename << endl;

  ofstream output(filename);
  if(!output){
    cout<< "Saving error "<< endl;
    exit(0);
  }
  output <<  size1 << endl;    

  for(uint i=0;i<size1;i++){
    output  <<  vmean[i] << endl;
  }
  output << size2 << endl;   
  for(uint i=0;i<size2;i++){
    output << vbase[i] << endl; 
  }
  output.close();
}






