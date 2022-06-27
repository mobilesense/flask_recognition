#include <stdio.h>


#include "ivfgeo.h"


int main (int argc, char ** argv)
{

  if(argc!=2) {
    fprintf(stderr,"usage: %s infilename.ivfgeo\n",argv[0]);
    return 0; 
  }

  char * fivfname = argv[1];

  /* Create the inverted file, and add the elements to it */
  fprintf (stderr, "* Read the Inverted file\n");
  ivfgeo_t * ivf = ivfgeo_read (fivfname, 0);

  if (ivf == NULL)
    return 1;

  ivfgeo_display (ivf);

  ivfgeo_delete (ivf);
  return 0;
}
