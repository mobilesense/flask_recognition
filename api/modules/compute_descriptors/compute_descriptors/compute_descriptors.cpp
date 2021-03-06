#include "../descriptor/feature.h"
#include "../util/util.h"
#include <string.h>

int main(int argc, char **argv){
     
  if(argc<2){     
    cout << "Interest point detectors/descriptors implemented by Krystian.Mikolajczyk@inrialpes.fr\n";
    cout << "at INRIA Rhone-Alpes.[ref. www.inrialpes.fr/movi/people/Mikolajczyk/Affine]" <<endl;
    cout << "Options:"<< endl;
    cout << "     -harlap - harris-laplace detector"  << endl;
    cout << "     -heslap - hessian-laplace detector"  << endl;
    cout << "     -haraff - harris-affine detector"  << endl;
    cout << "     -hesaff - hessian-affine detector"  << endl;
    cout << "     -harhes - harris-hessian-laplace detector"  << endl;
    cout << "     -dense xstep ystep nbscales - dense detector"  << endl;
    cout << "     -sedgelap - edge-laplace detector"  << endl;
    cout << "     -jla  - steerable filters,  similarity= " << endl; 
    cout << "     -sift - sift [D. Lowe],  similarity=" << endl;
    cout << "     -siftnonorm - sift [D. Lowe],  without normalization" << endl;
    cout << "     -msift - Mahalanobis sift, similarity= " << endl;
    cout << "     -gloh - extended sift,  similarity= " << endl;
    cout << "     -mom  - moments,  similarity= " << endl;
    cout << "     -koen - differential invariants,  similarity= " << endl;
    cout << "     -cf   - complex filters [F. Schaffalitzky],  similarity=" << endl;
    cout << "     -sc   - shape context,  similarity=45000 " << endl;
    cout << "     -spin - spin,  similarity= "  << endl;
    cout << "     -gpca - gradient pca [Y. Ke],  similarity="  << endl;
    cout << "     -cc - cross correlation,  similarity="  << endl;
    cout << "     -i image.pgm  - input image pgm, ppm, png" << endl;
    cout << "     -i2 image.jpg - input of any format of ImageMagick (WARN: uses only green channel)" << endl;
    cout << "     -pca input.basis - projects the descriptors with pca basis"  << endl;
    cout << "     -p1 image.pgm.points - input regions format 1" << endl; 
    cout << "     -p2 image.pgm.points - input regions format 2" << endl; 
    cout << "     -o1 out.desc - saves descriptors in out.desc output format1" << endl; 
    cout << "     -o2 out.desc - saves descriptors in out.desc output format2" << endl; 
    cout << "     -o3 out.siftbin - saves descriptors in binary format" << endl;
    cout << "     -o4 out.siftgeo - binary descriptor format used by Bigimbaz" << endl;
    cout << "     -coord out.coord - saves coordinates in binary format" << endl;
    cout << "     -noangle - computes rotation variant descriptors (no rotation esimation)" << endl; 
    cout << "     -DC - draws regions as circles in out.desc.png" << endl; 
    cout << "     -DR - draws regions as ellipses in out.desc.png" << endl; 
    cout << "     -c 255 - draws points in grayvalue [0,...,255]" << endl; 
    cout << "     -thres - threshod value" << endl;
    cout << "     -max - maximum for the number of computed descriptors in HESSAFF" << endl;
    cout <<"example:\n     "<< argv[0]<< " -jla -i image.png -p1 image.png.points -DR " << endl;
    cout <<"               "<< argv[0]<< "-harlap -gloh -i image.png  -DR " << endl << endl;
    cout <<" file format 1:"<< endl;
    cout <<"vector_dimension" << endl;
    cout << "nb_of_descriptors" << endl;
    cout << "x y a b c desc_1 desc_2 ......desc_vector_dimension"<< endl;
    cout <<"--------------------" << endl << endl;
    cout <<"where a(x-u)(x-u)+2b(x-u)(y-v)+c(y-v)(y-v)=1 " <<  endl <<  endl;
    cout <<" file format 2:"<< endl;
    cout <<"vector_dimension"  << endl;
    cout <<"nb_of_descriptors" << endl;
    cout <<"x y cornerness scale angle object_index  point_type laplacian extremum_type mi11 mi12 mi21 mi22 desc_1 ...... desc_vector_dimension"<< endl;
    cout <<"--------------------" << endl << endl;
    cout <<"distance=(descA_1-descB_1)^2+...+(descA_vector_dimension-descB_vector_dimension)^2"<< endl<< endl;
    cout <<" input.basis format:"<< endl;
    cout <<"nb_of_dimensions"<< endl;
    cout <<"mean_v1"<< endl;
    cout <<"mean_v2"<< endl;
    cout <<"."<< endl;
    cout <<"."<< endl;
    cout <<"mean_v_nb_of_dim"<< endl;
    cout <<"nb_of_dimensions*nb_of_pca_vectors"<< endl;
    cout <<"pca_vector_v1"<< endl;
    cout <<"pca_vector_v2"<< endl;
    cout <<"."<< endl;
    cout <<"."<< endl;
    cout <<"--------------------" << endl;    
    exit(-1);     
  }         
  
  char input[512];
  char pcafile[512];
  char input_desc[512];
  char output[512];  
  char coord_file[512]; 
  int in=0,p=2,out=0,coord=0,dr=0,color=255, fpca=0, filter=0, detector=0,descriptor=0,noangle=0, max_desc=0;
  float cluster_thres=-1000000000, thres=500;

  int dense_xstep=4, dense_ystep=4, dense_pyrlevels=5;

  for(int i=1;i<argc;i++){ 
    if(!strcmp(argv[i],"-i")){
      in = 1; sprintf(input,argv[i+1]);
    }else if(!strcmp(argv[i],"-i2")){
      in = 2; sprintf(input, argv[i+1]);
    }else if(!strcmp(argv[i],"-p1")){
      p=1;
      sprintf(input_desc,argv[i+1]);
    }else if(!strcmp(argv[i],"-p2")){
      p=0;
      sprintf(input_desc,argv[i+1]);      
    }else if(!strcmp(argv[i],"-o1")){
      out=1;sprintf(output,argv[i+1]);
    }else if(!strcmp(argv[i],"-o2")){
      out=2;sprintf(output,argv[i+1]);
    }else if(!strcmp(argv[i],"-o3")){
      out=3;sprintf(output,argv[i+1]);
    }else if(!strcmp(argv[i],"-o4")){
      out=4;sprintf(output,argv[i+1]);
    }else if(!strcmp(argv[i],"-coord")){
      coord=1;sprintf(coord_file,argv[i+1]);
    }else if(!strcmp(argv[i],"-filter")){
      filter=1;
      cluster_thres=atof(argv[i+1]);     
    }else if(!strcmp(argv[i],"-thres")){
      thres=atof(argv[i+1]);     
    }else if(!strcmp(argv[i],"-noangle")){
      noangle=1;
    }else if(!strcmp(argv[i],"-DR")){
      dr=1; 
    }else if(!strcmp(argv[i],"-DC")){
      dr=2; 
    }else if(!strcmp(argv[i],"-c")){
      color=atoi(argv[i+1]);
    }else if(!strcmp(argv[i],"-max")){
      max_desc=atoi(argv[i+1]);
    }else if(!strcmp(argv[i],"-pca")){
      fpca=1;
      sprintf(pcafile,argv[i+1]);
    }
  }  
 
  DARY *image;
  if(in){
    cout << "\ncomputing descriptors in image " << input <<endl;

    if(in==1) image = new ImageContent(input);
    else {
      image = new ImageContent(input,in);
    }

    if(image->x()<=12 || image->y()<=12) {
      fprintf(stderr,"image too small, aborting\n");
      exit(2);
    }
    image->toGRAY();   
    image->char2float();      
  }else {
    cout << "No input image "<< endl;
    exit(-1);
  }
  
  int aff=0;
  cout << "detector: ";
  for(int i=1;i<argc;i++){ 
    if(!strcmp(argv[i],"-harlap")){
      detector=HARRIS;if(!out)sprintf(input,"%s.harlap",input);cout << "harlap ";
    }else if(!strcmp(argv[i],"-heslap")){
      detector=HESSIAN;if(!out)sprintf(input,"%s.heslap",input);cout << "heslap ";
    }else if(!strcmp(argv[i],"-edgelap")){
      detector=EDGE;if(!out)sprintf(input,"%s.edgelap",input);cout << "edgelap ";
    }else if(!strcmp(argv[i],"-haraff")){
      detector=HARRIS;aff=16;if(!out)sprintf(input,"%s.haraff",input);cout << "haraff ";
    }else if(!strcmp(argv[i],"-hesaff")){
      detector=HESSIAN;aff=16;if(!out)sprintf(input,"%s.hesaff",input);cout << "hesaff ";
    }else if(!strcmp(argv[i],"-harhes")){
      detector=HARHES;if(!out)sprintf(input,"%s.harhes",input);cout << "harheslap ";
    }else if(!strcmp(argv[i],"-sedgelap")){
      detector=SEDGE;if(!out)sprintf(input,"%s.sedgelap",input);cout << "sedgelap ";
    } else if(!strcmp(argv[i],"-dense") && i+3<argc &&
              sscanf(argv[i+1],"%d",&dense_xstep)==1 && 
              sscanf(argv[i+2],"%d",&dense_ystep)==1 && 
              sscanf(argv[i+3],"%d",&dense_pyrlevels)==1) {
      detector=DENSE; 
      printf("dense %d %d %d",dense_xstep, dense_ystep, dense_pyrlevels);     
    }
  }   
  if(detector==0)cout << "non";
  cout << endl; 

 
  cout << "descriptor: ";
  for(int i=1;i<argc;i++){ 
    if(!strcmp(argv[i],"-jla")){
      descriptor=JLA;if(!out)sprintf(output,"%s.jla",input);cout << "jla";
    }else if(!strcmp(argv[i],"-koen")){
      descriptor=KOEN;if(!out)sprintf(output,"%s.koen",input);cout << "koen ";
    }else if(!strcmp(argv[i],"-cf")){
      descriptor=CF;if(!out)sprintf(output,"%s.cf",input);cout << "cf ";
    }else if(!strcmp(argv[i],"-sc")){
      descriptor=SHAPE;if(!out)sprintf(output,"%s.sc",input);cout << "sc ";
    }else if(!strcmp(argv[i],"-sift")){
      descriptor=SIFT;if(!out)sprintf(output,"%s.sift",input);cout << "sift ";
    }else if(!strcmp(argv[i],"-siftnonorm")){
      descriptor=SIFTNONORM;if(!out)sprintf(output,"%s.siftnonorm",input);cout << "siftnonorm ";
    }else if(!strcmp(argv[i],"-msift")){
      descriptor=MSIFT;if(!out)sprintf(output,"%s.msift",input);cout << "msift ";
    }else if(!strcmp(argv[i],"-gloh")){
      descriptor=GLOH;if(!out)sprintf(output,"%s.gloh",input);cout << "gloh ";
    }else if(!strcmp(argv[i],"-gpca")){
      descriptor=PCA;if(!out)sprintf(output,"%s.pca",input);cout << "pca "; 
    }else if(!strcmp(argv[i],"-cc")){
      descriptor=CC;if(!out)sprintf(output,"%s.cc",input);cout << "cc ";
    }else if(!strcmp(argv[i],"-mom")){
      descriptor=MOM;if(!out)sprintf(output,"%s.mom",input);cout << "mom ";
    }else if(!strcmp(argv[i],"-spin")){
      descriptor=SPIN;if(!out)sprintf(output,"%s.spin",input);cout << "spin ";
    }else if(!strcmp(argv[i],"-pca")){
      fpca=1;
    }
  }
  if(descriptor==0)cout << "non";
  cout << endl;
  cout <<"threshold: "<< thres <<   endl;

  vector<FeatureDescriptor*> desc;   

  if(p==2 && detector==0){
    cout << "neihter detector nor input feature file was selected"<< endl;
    exit(-1);
  }else if(p!=2 && detector>0){
    cout << "eihter detector or input feature file can be selected"<< endl;
    exit(-1);
  }else if((descriptor==0)){
    cout << "WARNING! no descriptor is selected."<< endl;
    // exit(-1);
  }else if(p!=2){   
    cout << "loading points from "<< input_desc<<"... ";  
    loadFeatures(input_desc,desc,p); 
    cout << desc.size()<< endl;
    if(noangle){
      for(uint c=0;c<desc.size();c++){	 
	desc[c]->setAngle(0.0);
      }
    }  
  }
 
  Timer timer;
  if(detector>0){
    if(detector!=DENSE) 
      extractFeatures(image, desc, detector, descriptor, aff, noangle, thres, max_desc);
    else 
      fastDense(image, desc, descriptor, dense_xstep, dense_ystep, dense_pyrlevels);
  }else {
    computeAffineDescriptors(image,desc,descriptor, 1);
  }

  
  if(filter){ 
    cout << "filtering features with threshold= "<< cluster_thres<<" ..."<< endl;
    //filterFeatures(desc, cluster_thres);
  }

  if(fpca){
    uint mean_dim; 
    uint base_dim;
    float *mvec;// = new float[mean_dim];
    float *base;// = new float[base_dim];
    loadPCAVectors(pcafile, mvec,mean_dim, base, base_dim);
    pca(desc,mvec,base, mean_dim);
    delete []mvec;
    delete []base;
  }
 
  timer.stop();
  cout << "saving " << desc.size() << " feaures in output file: " << output<< endl;
  /*
  vector<FD *> aggoOld;
  vector<ClStep *> vClusterTrace;
  int K=desc.size()/500;
  if(K<2)K=2;
  kmeansAgglo(desc, aggoOld, vClusterTrace, K , CLUSTER_SIMILARITY);
  timer.stop();
  */
 if(!out)
    writeFeatures(desc, output,1); 
  if(out==2)
    writeFeatures(desc, output,0);  
  if(out==1)
    writeFeatures(desc, output,1);
  if(out==3)
    writeFeaturesHerve(desc, output);
  if(out==4)
    writeFeaturesSameh(desc, output);
  
  if(coord==1)
    writeCoordinates(desc, coord_file);

  if(dr){
    char draw[512];
    sprintf(draw,"%s.png",output);  
    cout << "drawing points in: " << draw<< endl;
    if(dr==1)
      displayFeatures(image, desc, draw,color);
    if(dr==2)
      displayFeatures(image, desc, draw,color,1);
  }
    
  delete image;
  desc.clear();
  return 0;
}
  
