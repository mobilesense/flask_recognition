
CFLAGS=-fPIC -Wall -g -O3 -lm

PYTHONCFLAGS=-I/usr/include/python2.7

YAELCFLAGS=-I/home/amine/Bureau/flask_recognition/api/modules/compile/..
YAELLDFLAGS=-L/home/amine/Bureau/flask_recognition/api/modules/compile/../yael -Wl,-rpath,/home/amine/Bureau/flask_recognition/api/modules/compile/../yael -lyael


SWIG=swig -python

WRAPLDFLAGS=-shared
LAPACKLDFLAGS=/usr/lib/libblas.so.3gf /usr/lib/liblapack.so.3gf

USETHREADS=yes
THREADCFLAGS=-DHAVE_THREADS
THREADLDFLAGS=-lpthread

USEFF=no
FFLDFLAGS=-lz -lavformat -lavcodec -lavdevice -lswscale -lavutil -lm
FFCFLAGS=-D_ISOC9X_SOURCE

USEMPI=no
MPILDFLAGS=-lmpi
MPICFLAGS=

SHAREDEXT=so
SHAREDFLAGS=-shared

