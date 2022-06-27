.SUFFIXES: .c .o .cu .cpp

include ../compile/makefile.inc



FLAGS=$(CFLAGS) -I../ 
LIBS=$(LDFLAGS) $(LAPACKLDFLAGS) $(THREADLDFLAGS) $(YAELLDFLAGS) "-lm"


OBJS=siftgeo.o siftgeo_binarize.o siftgeo_and_nn.o siftgeo_agg.o 
CLUSTOBJS=../utils/hkm.o ../utils/annvw.o 


ifeq "$(USEMPI)" "yes"
  EXTRATARGET=cluster_mpi
endif


all: ann_vw cluster  siftgeo_disp _siftgeo.so $(EXTRATARGET)



############################################################
# generic rules


.c.o:
	$(CC) -c $< -o $@ $(FLAGS) $(EXTRACFLAGS) -lm

############################################################
# specific rules


siftgeo_binarize.o: EXTRACFLAGS=$(LAPACKCFLAGS)

ann_vw.o: EXTRACFLAGS=$(THREADCFLAGS)

cluster: cluster.o siftgeo.o siftgeo_binarize.o  $(CLUSTOBJS) siftgeo_and_nn.o 
	$(CC) -o $@ $^  $(LIBS)


siftgeo_disp: siftgeo_disp.o siftgeo.o 
	$(CC) -o $@ $^  $(LIBS)

ann_vw: ann_vw.o siftgeo.o  $(CLUSTOBJS) siftgeo_and_nn.o
	$(CC) -o $@ $^  $(LIBS) 

cluster_mpi: cluster_mpi.o siftgeo.o $(CLUSTOBJS) siftgeo_and_nn.o siftgeo_binarize.o 
	$(CC) -o $@ $^  $(LIBS) $(MPILDFLAGS)

cluster_mpi.o: EXTRACFLAGS=$(MPICFLAGS)

############################################################
# wrapper


siftgeo_wrap.c : siftgeo.swg siftgeo.h siftgeo_binarize.h siftgeo_agg.h ../yael/common.swg siftgeo_and_nn.h siftgeo_agg.h 
	$(SWIG) -I.. $<

_siftgeo.so: siftgeo_wrap.o $(OBJS) 
	gcc -o $@  $(WRAPLDFLAGS) $^ $(LIBS)

siftgeo_wrap.o: EXTRACFLAGS=$(PYTHONCFLAGS)


############################################################
# indent


indent: nn.c ../utils/nn.h siftgeo.c siftgeo.h 
	$(INDENT) $^

############################################################
# dependencies



# generate dependencies with
# for i in *.c; do cpp -I.. -MM $i; done
ann_vw.o: ann_vw.c ../yael/clustering.h ../yael/vector.h ../yael/nn.h \
  ../yael/machinedeps.h ../yael/sorting.h ../utils/annvw.h \
  ../siftgeo/siftgeo.h ../siftgeo/siftgeo_and_nn.h ../siftgeo/siftgeo.h \
  siftgeo_binarize.h
cluster.o: cluster.c ../yael/machinedeps.h ../yael/vector.h \
  ../yael/nn.h ../yael/clustering.h ../utils/hkm.h ../siftgeo/siftgeo.h \
  ../siftgeo/siftgeo_and_nn.h ../siftgeo/siftgeo.h siftgeo_binarize.h
cluster_mpi.o: cluster_mpi.c ../yael/nn.h ../yael/vector.h \
  ../siftgeo/siftgeo.h ../siftgeo/siftgeo_binarize.h
siftgeo.o: siftgeo.c ../siftgeo/siftgeo.h
siftgeo_agg.o: siftgeo_agg.c ../yael/vector.h ../yael/sorting.h \
  siftgeo_agg.h siftgeo.h
siftgeo_and_nn.o: siftgeo_and_nn.c siftgeo_and_nn.h siftgeo.h \
  ../yael/vector.h
siftgeo_binarize.o: siftgeo_binarize.c ../yael/vector.h ../yael/matrix.h \
  ../yael/nn.h ../yael/sorting.h ../siftgeo/siftgeo.h \
  ../siftgeo/siftgeo_binarize.h
siftgeo_disp.o: siftgeo_disp.c ../siftgeo/siftgeo.h


clean: 
	rm -f *.o cluster  siftgeo.py _siftgeo.so *_wrap.c \
                  siftgeo_disp ann_vw *.pyc cluster_mpi


