#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <siftgeo/siftgeo.h>


/* Display a set of descriptors */
void display_points_custom (const point_t * des, int n)
{
  int i, j;

  for (i = 0; i < n; i++) {
    const geom_t *g = &des[i].geom;

    fprintf (stdout,
             "%.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f "
             "%.3f", g->x, g->y, g->scale, g->angle,
             g->mi11, g->mi12, g->mi21, g->mi22, g->cornerness);

    if (des[i].desc == NULL) {
      fprintf (stdout, " %-7d", des[i].vw);
      fprintf (stdout, " %llx", des[i].binsign);
    } else
      for (j = 0 ; j < des[i].dim ; j++)
	fprintf (stdout, " %d", des[i].desc[j]);

    fprintf (stdout, "\n");
  }
}


int main (int argc, char **argv) 
{
  int i = 0;
  FILE * fi = stdin;

  int mode = 0;
  int display_mode = 0;

  if (argc > 1) {
    if (strcmp (argv[1], "-vwgeo") == 0)
      mode = 1;
    else if (strcmp (argv[1], "-vwsgeo") == 0)
      mode = 2;
  }

  if (argc > 2) {
    if (strcmp (argv[2], "-compact") == 0)
      display_mode = 1;
  }

  while (!feof (fi)) {
    point_t * siftgeo = (point_t *) malloc (sizeof (point_t));

    int ret = read_point_t (fi, siftgeo, mode);
    if (ret != 1)
      break;

    if (display_mode == 0) {
      fprintf (stdout, "des=%-5d ", i);
      display_points (siftgeo, 1);
    }
    else if (display_mode == 1)
      display_points_custom (siftgeo, 1);

    delete_points (siftgeo, 1);
    i++;
  }
  return 0;
}


