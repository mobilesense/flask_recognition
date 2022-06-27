
include ../compile/makefile.inc

.SUFFIXES: .c .o .h

FLAGS=$(CFLAGS) -I.. $(YAELCFLAGS)
LIBS=$(LDFLAGS) -lm  $(YAELLDFLAGS)


OBJS=ivfgeo.o ../siftgeo/siftgeo.o
PQOBJS=ivfgeo_pq.o ../utils/pq.o ../utils/ivfpq.o ../siftgeo/siftgeo_and_nn.o 


all: ivfgeo_display ivfgeo.py _ivfgeo.so ivfgeo_merge 


ivfgeo_display: ivfgeo_display.o  $(OBJS) 
	$(CC) $(FLAGS) -o $@ $^  $(LIBS)

ivfgeo_summary: ivfgeo_summary.o  $(OBJS) 
	$(CC) $(FLAGS) -o $@ $^  $(LIBS)

ivfgeo_merge: ivfgeo_merge.o  $(OBJS) 
	$(CC) $(FLAGS) -o $@ $^  $(LIBS)

.c.o:
	$(CC) -c $< -o $@ $(FLAGS) $(EXTRACFLAGS)

indent: ivfgeo.h ivfgeo.c
	$(INDENT) $^

######################################################
# wrapper
 
ivfgeo_wrap.c ivfgeo.py: ivfgeo.swg ivfgeo.h ../yael/common.swg ivfgeo_pq.h
	$(SWIG) -I.. $< 

_ivfgeo.so: ivfgeo_wrap.o $(OBJS)  $(PQOBJS)
	$(CC) -o $@ $(WRAPLDFLAGS) $^ $(LIBS) $(LAPACKLDFLAGS)

ivfgeo_wrap.o: EXTRACFLAGS=$(PYTHONCFLAGS)



# generate dependencies with
# for i in *.c; do cpp -I.. -MM $i; done
ivfgeo.o: ivfgeo.c ../siftgeo/siftgeo.h ../yael/sorting.h \
  ../yael/machinedeps.h ivfgeo.h
ivfgeo_display.o: ivfgeo_display.c ivfgeo.h ../siftgeo/siftgeo.h
ivfgeo_merge.o: ivfgeo_merge.c ../ivfgeo/ivfgeo.h ../siftgeo/siftgeo.h
ivfgeo_pq.o: ivfgeo_pq.c ../yael/vector.h ../yael/sorting.h \
  ../siftgeo/siftgeo_and_nn.h ../siftgeo/siftgeo.h ivfgeo_pq.h \
  ../utils/ivfpq.h ../utils/pq.h ../siftgeo/siftgeo.h ivfgeo.h
ivfgeo_summary.o: ivfgeo_summary.c ivfgeo.h ../siftgeo/siftgeo.h

clean: 
	rm -f *.o ivfgeo_query ivfgeo_create ivfgeo_display ivfgeo_merge _ivfgeo.so *.pyc ivfgeo.py ivfgeo_wrap.c

