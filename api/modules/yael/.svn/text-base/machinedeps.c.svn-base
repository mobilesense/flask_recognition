#include <math.h>
#include <sys/types.h>


#include "machinedeps.h"


#ifdef __linux__

#define __USE_GNU
#include <sched.h>

int count_cpu (void)
{
  cpu_set_t set;
  sched_getaffinity (0, sizeof (cpu_set_t), &set);
  int i, count = 0;
  for (i = 0; i < CPU_SETSIZE; i++)
    if (CPU_ISSET (i, &set))
      count++;
  return count;
}


#elif defined(__APPLE__)

#include <sys/types.h>
#include <sys/sysctl.h>


int count_cpu (void) {
  int count=-1;
  size_t count_size=sizeof(count);
  sysctlbyname("hw.ncpu",&count,&count_size,NULL,0);
  return count;
}

#else

int count_cpu() {
  return 1;
}


#endif

#ifndef __APPLE__

double log2(double x) {
  return log(x)/M_LN2;
}


#endif

