#include <sys/time.h>
#include <sys/resource.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <errno.h>
#include <stdio.h>
#include <string.h>

int main(int argc, char **argv)
{
  /* Do not optimize this var, even if it's value is never used */
  volatile char total = 0;

  struct stat    stats;
  struct rlimit  limit;

  char *name;
  char *data;
  int   fd;
  long  idx;
  long  page;

  /* Command line check */
  if (argc < 2) {
    fprintf(stderr, "Usage: %s <ivf file>\n", argv[0]);
    return -1;
  }

  name = argv[1];

  /* Get file size */
  if (stat(name, &stats) < 0) {
    fprintf(stderr, "%s: %s\n", name, strerror(errno));
    return -2;
  }

  /* Set max locked memory to unlimited - probably requires to be root */
  limit.rlim_cur = RLIM_INFINITY;
  limit.rlim_max = RLIM_INFINITY;
  if (setrlimit(RLIMIT_MEMLOCK, &limit) < 0) {
    fprintf(stderr, "setrlimit: %s\n", strerror(errno));
    return -2;    
  }

  /* Open and map file into memory */
  fd = open(name, O_RDONLY);
  if (fd < 0) {
    fprintf(stderr, "%s: %s\n", name, strerror(errno));
    return -2;
  }

  data = mmap(NULL, stats.st_size, PROT_READ, MAP_SHARED, fd, 0);

  if (data == MAP_FAILED) {
    fprintf(stderr, "mmap: %s\n", strerror(errno));
    return -2;    
  }

  /* Prevent data to be swapped */
  if (mlock(data, stats.st_size) < 0) {
    fprintf(stderr, "mlock: %s\n", strerror(errno));
    return -2;    
  }

  /* Trigger the load of the data */
  page = sysconf(_SC_PAGESIZE);
  for (idx = 0; idx < stats.st_size; idx += page) total += data[idx];

  /* Infinite sleep as this process is expected to be terminated by a signal only */
  select(0, NULL, NULL, NULL, NULL);

  return 0;
}
