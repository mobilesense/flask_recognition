#!/bin/bash

# Customize makefile.inc


# Detect architecture
if [ `uname` == Darwin ] ; then
  conf=mac
elif [ `uname -m` == 'x86_64' ] ; then
  conf=linux64
elif  [ `uname -m` == 'i686' ] || [ `uname -m` == 'i386' ]  ; then
  conf=linux32
fi


# Defaults parameters
cflags='-fPIC -Wall -g -O3 -lm'
ldflags='-g -fPIC -lm'

yaelprefix=$PWD/../yael


lapackldflags='-lblas -llapack'
if [ -e /usr/lib/libblas.so.3gf ]; then
   # special names for libs on ubuntu
   lapackldflags="/usr/lib/libblas.so.3gf /usr/lib/liblapack.so.3gf"
elif [ -d /usr/lib64/atlas/ ]; then 
   lapackldflags="/usr/lib64/atlas/libblas.so.3 /usr/lib64/atlas/liblapack.so.3"
elif [ -d /usr/lib/atlas/ ]; then
   lapackldflags="/usr/lib/atlas/libblas.so.3 /usr/lib/atlas/liblapack.so.3"
else
   echo -n "using default locations for blas and lapack"
fi





useff=no
ffldflags='-lz -lavformat -lavcodec -lavdevice -lswscale -lavutil -lm'
ffcflags=-D_ISOC9X_SOURCE


usempi=no
mpicflags=""
mpildflags="-lmpi"

if [ $conf == mac ]; then
  wrapldflags="-Wl,-F. -bundle -undefined dynamic_lookup"
  sharedext=dylib
  sharedflags="-dynamiclib -install_name $yaelprefix/libyael.dylib"
else
  wrapldflags="-shared"
  sharedflags="-shared"
  sharedext=so  
fi

function usage () {
    cat <<EOF 1>&2
usage: $0 
  [--debug] 
  [--yael=yaelprefix] 
  [--swig=swigpath]
  [--lapack=ldflagsforlapack+blas]
EOF
    exit $1
}



# Search latest python version
for pysubver in {7,6,5,4,x} ; do
    if [ -f "/usr/include/python2.${pysubver}/Python.h" ] ; then
	echo "Found python development version 2.$pysubver"
	break
    fi
done


if [ "x$pysubver" == "xx" ] ; then
    echo "# Error: no python directory (python-dev) found"
    exit 1
fi



# Parse command line arguments
while [ $# -gt 0 ] ; do

    a=$1
    shift

    # echo $a

    case $a in 
	-h | --help) usage 0 ;; 

	--debug) cflags="${cflags/ -O3/}" ;;
	
	--yael=*) yaelprefix=${a#--yael=};;
        --swig=*) swig=${a#--swig=} ;;
	--lapack=*) lapackldflags=${a#--lapack=} ;;
        --enable-ffmpeg) useff=yes ;;
	--ffmpeg=*) ffprefix=${a#--ffmpeg=}
	    ffldflags="-L$ffprefix/lib $ffldflags"
	    ffcflags="-I$ffprefix/include $ffcflags" ;;
        --enable-mpi) usempi=yes ;;
	--mpi=*) mpiprefix=${a#--mpi=}
	    mpildflags="-L$mpiprefix/lib $mpildflags"
	    mpicflags="-I$mpiprefix/include $mpicflags" ;;

	*)  echo "unknown option $a" 1>&2; exit 1;;
    esac
done

yaellib=${yaelprefix}
yaelinc=${yaelprefix%/*}


if [ -z "$swig" ]; then
  if which swig ; then
    swig=swig
  else 
    echo "Error: no swig executable found. Provide one with --swig"
    exit 1
  fi
fi


cat <<EOF | tee makefile.inc

CFLAGS=$cflags

PYTHONCFLAGS=-I/usr/include/python2.$pysubver

YAELCFLAGS=-I$yaelinc
YAELLDFLAGS=-L$yaellib -Wl,-rpath,$yaellib -lyael


SWIG=$swig -python

WRAPLDFLAGS=$wrapldflags
LAPACKLDFLAGS=$lapackldflags

USETHREADS=yes
THREADCFLAGS=-DHAVE_THREADS
THREADLDFLAGS=-lpthread

USEFF=$useff
FFLDFLAGS=$ffldflags
FFCFLAGS=$ffcflags

USEMPI=$usempi
MPILDFLAGS=$mpildflags
MPICFLAGS=$mpicflags

SHAREDEXT=$sharedext
SHAREDFLAGS=$sharedflags

EOF

