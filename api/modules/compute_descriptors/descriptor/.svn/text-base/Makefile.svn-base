 LIB =    libdescriptor.a
FILES   =  FeatureDescriptor  SiftDescriptor \
	   KoenDescriptor JetLocalAffine \
	   FeatureDetector  \
           MomentDescriptor ComplexGauss \
	   SpinDescriptor DisplayFeature   \
	   EdgeDetector  \
	   CCDescriptor	
	   
include ../makefile.inc

CCDescriptor.o: CCDescriptor.cpp feature.h ../ImageContent/imageContent.h \
  cc_base.h cc_pca_base.h ../gauss_iir/gauss_iir.h \
  ../gauss_iir/../ImageContent/imageContent.h
ComplexGauss.o: ComplexGauss.cpp feature.h ../ImageContent/imageContent.h \
  cf_base.h cf_pca_base.h complexFilters.h
DisplayFeature.o: DisplayFeature.cpp feature.h \
  ../ImageContent/imageContent.h ../matrix/matrix.h ../matrix/matrix.hh \
  ../matrix/./matrix.icc ../matrix/vector.hh ../matrix/./vector.icc \
  ../matrix/points.hh ../matrix/quatern.hh
EdgeDetector.o: EdgeDetector.cpp feature.h ../ImageContent/imageContent.h \
  ../gauss_iir/gauss_iir.h ../gauss_iir/../ImageContent/imageContent.h
FeatureDescriptor.o: FeatureDescriptor.cpp feature.h \
  ../ImageContent/imageContent.h ../matrix/matrix.h ../matrix/matrix.hh \
  ../matrix/./matrix.icc ../matrix/vector.hh ../matrix/./vector.icc \
  ../matrix/points.hh ../matrix/quatern.hh ../util/util.h \
  ../gauss_iir/gauss_iir.h ../gauss_iir/../ImageContent/imageContent.h
FeatureDetector.o: FeatureDetector.cpp feature.h \
  ../ImageContent/imageContent.h ../gauss_iir/gauss_iir.h \
  ../gauss_iir/../ImageContent/imageContent.h ../util/util.h
JetLocalAffine.o: JetLocalAffine.cpp feature.h \
  ../ImageContent/imageContent.h gaussFilters.h jla_base.h jla_pca_base.h
KoenDescriptor.o: KoenDescriptor.cpp feature.h \
  ../ImageContent/imageContent.h koen_base.h koen_pca_base.h
MomentDescriptor.o: MomentDescriptor.cpp feature.h \
  ../ImageContent/imageContent.h mom_base.h mom_pca_base.h \
  ../gauss_iir/gauss_iir.h ../gauss_iir/../ImageContent/imageContent.h
SiftDescriptor.o: SiftDescriptor.cpp feature.h \
  ../ImageContent/imageContent.h ../gauss_iir/gauss_iir.h \
  ../gauss_iir/../ImageContent/imageContent.h sift_base.h sift_pca_base.h \
  esift_base.h esift_pca_base.h sc_base.h sc_pca_base.h pca_base.h \
  pca_eig.h
SpinDescriptor.o: SpinDescriptor.cpp feature.h \
  ../ImageContent/imageContent.h spin_base.h spin_pca_base.h
