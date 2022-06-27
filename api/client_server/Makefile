
include ../makefile.inc



all: invfile.py _invfile.so compute_descriptors



INCLUDES=-I../invfile -I../siftgeo -I../hash $(ITCFLAGS)


invfile_wrap.c invfile.py: invfile.swg \
	../siftgeo/ann.h ../invfile/invfile.h ../hash/hash.h\
	../siftgeo/siftgeo.h
	swig -python $(LAPACKCFLAGS) $(INCLUDES) invfile.swg


invfile_wrap.o: invfile_wrap.c 
	gcc -g -c -fpic $< -o $@ \
             -I$(PYTHONINCLUDE) $(INCLUDES)

OBJS=../invfile/invfile.o ../hash/hash.o ../siftgeo/siftgeo.o \
../siftgeo/ann.o ../lattice/lattice.o  ../siftgeo/clustering.o\
 $(addprefix ../siftgeo/,$(CUDAOBJ))


_invfile.so: invfile_wrap.o $(OBJS)
	gcc -o $@  -shared $^ $(ITLDFLAGS) $(LAPACKLDFLAGS) $(CUDALDFLAGS)

compute_descriptors: $(MIKOCOMPUTEDESC)
	ln -s $< $@


clean:
	rm -f invfile_wrap.c invfile.py invfile_wrap.o _invfile.so compute_descriptors *.o


