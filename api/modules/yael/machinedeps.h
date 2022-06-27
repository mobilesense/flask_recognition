#ifndef MACHINEDEPS_H_INCLUDED
#define MACHINEDEPS_H_INCLUDED

#ifdef __APPLE__
#define HAVE_QSORT_R
#endif

#ifdef __linux__
#define HAVE_TLS
#else
#define __thread 
#endif

int count_cpu(void);

#ifndef __APPLE__

double log2(double x);

#endif

#ifdef __linux__
#include <malloc.h>
#else
#include <stdlib.h>

static void *memalign (size_t ignored, size_t nbytes)
{
  return malloc (nbytes);
}
#endif



#endif

