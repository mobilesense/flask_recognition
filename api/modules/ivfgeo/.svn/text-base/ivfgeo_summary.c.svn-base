#include <stdlib.h>
#include <stdio.h>

#include "ivfgeo.h"


int main (int argc, char ** argv)
{
  if (argc < 2) {
    fprintf(stderr, "Usage: %s <ivfgeo file> [<format>]\n", argv[0]);
    return 1;
  }
  char *fivfname = argv[1];
  int format = 0;

  if (argc >= 3) {
    format = atoi(argv[2]);
  }

  /* Create the inverted file, and add the elements to it */
  ivfgeo_t * ivf = ivfgeo_read (fivfname, format);

  if (ivf == NULL)
    return 2;

  ivfgeo_summary (ivf);

  ivfgeo_delete (ivf);
  return 0;
}
