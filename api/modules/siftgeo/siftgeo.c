#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include <siftgeo/siftgeo.h>


/*---------------------------------------------------------------------------*/
/* General purpose functions                                                 */
/*---------------------------------------------------------------------------*/

pointset_t *pointset_alloc (int npt)
{
  vwgeoset_t *ret = pointset_new ();
  ret->n = npt;
  ret->pts = calloc (sizeof (point_t), npt);
  return ret;
}

vwgeoset_t *pointset_new ()
{
  vwgeoset_t *ret = malloc (sizeof (pointset_t));
  ret->pts = NULL;
  ret->n = 0;
  return ret;
}


pointset_t *pointset_from_vw (int *vw, int n)
{
  pointset_t *ps = pointset_alloc (n);
  int i;
  for (i = 0; i < n; i++)
    ps->pts[i].vw = vw[i];
  return ps;
}

int *vw_from_pointset(const pointset_t *ps) {
  int *vw=malloc(sizeof(int)*ps->n);
  int i;
  
  for(i=0;i<ps->n;i++) 
    vw[i]=ps->pts[i].vw;

  return vw;
}


void delete_points (point_t * corners, int nc)
{
  int i;
  for (i = 0; i < nc; i++)
    free (corners[i].desc);
  free (corners);
}


void pointset_delete (vwgeoset_t * vgs)
{
  delete_points (vgs->pts, vgs->n);
  free (vgs);
}


/* Compute the norm of a given set of vwgeo */
double vwgeoset_norm1 (const vwgeoset_t * vwg)
{
  return vwg->n;
}


/* Compute the norm of a given set of vwgeo */
double vwgeoset_norm2 (const vwgeoset_t * vwg)
{
  int i;
  double s = 0, stmp = 0;

  if (vwg->n > 0)
    stmp = 1;

  for (i = 1; i < vwg->n; i++) {
    if (vwg->pts[i].vw == vwg->pts[i - 1].vw)
      stmp++;
    else {
      assert (vwg->pts[i].vw > vwg->pts[i - 1].vw);
      s += stmp * stmp;
      stmp = 1;
    }
  }

  s += stmp * stmp;

  return sqrt(s);
}


/* make new set of vwgeo's with geom_t from the descriptors and vw from vws 
 * k=# of vw per point_t */
point_t *siftgeo_to_vwgeo (point_t * des, int n, int *vws, int k)
{
  int i;
  point_t *ret = malloc (sizeof (point_t) * n * k);
  for (i = 0; i < n; i++) {
    int j;
    for (j = 0; j < k; j++) {
      int ii = i * k + j;
      ret[ii].geom = des[i].geom;
      ret[ii].binsign = des[i].binsign;
      ret[ii].vw = vws[ii];
      ret[ii].desc = NULL;
    }
  }
  return ret;
}


void pointset_append (pointset_t * ps, pointset_t * ps2)
{
  int i;
  ps->pts = realloc (ps->pts, sizeof (point_t) * (ps->n + ps2->n));

  for (i = 0; i < ps2->n; i++)
    ps->pts[i + ps->n] = ps2->pts[i];
  ps->n += ps2->n;
  ps2->n = 0;

}


void pointset_affine_transform (pointset_t * ps, double aff[6], double scalef)
{
  int i;
  for (i = 0; i < ps->n; i++) {
    geom_t *g = &ps->pts[i].geom;
    double x = g->x, y = g->y;

    g->x = aff[0] * x + aff[1] * y + aff[2];
    g->y = aff[3] * x + aff[4] * y + aff[5];

    g->scale *= scalef;

  }
}


pointset_t *pointset_dup (pointset_t * ps)
{
  pointset_t *ps_new = pointset_new ();

  ps_new->n = ps->n;
  ps_new->pts = malloc (sizeof (ps_new->pts[0]) * ps->n);
  memcpy (ps_new->pts, ps->pts, sizeof (ps_new->pts[0]) * ps->n);
  int i;
  for (i = 0; i < ps->n; i++) {
    point_t *p = &ps->pts[i];
    if (p->desc) {
      point_t *p_new = &ps_new->pts[i];
      p_new->desc = malloc (sizeof (p->desc[0]) * p->dim);
      memcpy (p_new->desc, p->desc, sizeof (p->desc[0]) * p->dim);
    }
  }

  return ps_new;
}


/*---------------------------------------------------------------------------*/
/* I/O                                                                       */
/*---------------------------------------------------------------------------*/

int read_point_t (FILE * f, point_t * corner, int fmt)
{

#define RETERR(no) { \
                     perror ("read_point_t error "#no);\
                     return -1;\
                   }

  corner->binsign = 0;

  if (fmt==0 || fmt==3) {                      /* siftgeo or siftbin */

    if(fmt==0) {
      int r = fread (&corner->geom, sizeof (geom_t), 1, f);
      if (r != 1) {
        if (feof (f))
          return 0;
        RETERR (3);
      }
    } else {
      int unused;
      int r = fread (&unused, sizeof (int), 1, f);
      if (r == 0) {
        if (feof (f))
          return 0;
        RETERR (3);
      }
      if (fread (&unused, sizeof (int), 1, f)!=1 || 
          fread (&corner->geom.x, sizeof (float), 1, f)!=1 || 
          fread (&corner->geom.y, sizeof (float), 1, f)!=1 || 
          fread (&corner->geom.scale, sizeof (float), 1, f)!=1 || 
          fread (&corner->geom.angle, sizeof (float), 1, f)!=1) {
        RETERR (2);     
      }
      corner->geom.mi11=1;
      corner->geom.mi12=0;
      corner->geom.mi21=0;
      corner->geom.mi22=1;
      corner->geom.cornerness=0;
    }     

    if (fread (&corner->dim, sizeof (int), 1, f) != 1)
      RETERR (4);

    int d = corner->dim;

    if (d <= 0 || d > 2048 * 1024) {
      fprintf (stderr, "read_point_t: weird descriptor dimension %d!\n", d);
      return -1;
    }

    corner->desc = malloc (d);
    assert (corner->desc);

    

    if (fread (corner->desc, 1, d, f) != d) {
      free (corner->desc);
      RETERR (5);
    }

    corner->binsign = 0;
  } else if (fmt==1 || fmt==2) { /* vwgeo and vwsgeo */
    int r = fread (&corner->vw, sizeof (int), 1, f);
    if (r == 0) {
      if (feof (f))
        return 0;
      RETERR (1);
    }

    if (fread (&corner->geom, sizeof (geom_t), 1, f) != 1)
      RETERR (2);

    if (fmt == 2) {
      if (fread (&corner->binsign, sizeof (binsign_t), 1, f) != 1)
        RETERR (6);
    }
    corner->desc = NULL;
  } else assert(0);
#undef RETERR
  return 1;
}



int read_points (const char *fname, point_t ** corners_out, int fmt)
{

  FILE *f = fopen (fname, "r");

  if (!f) {
    fprintf (stderr, "could not open %s\n", fname);
    perror ("");
    return -1;
  }
  int n = read_points_file (f, corners_out, fmt);

  fclose (f);

  return n;
}


/* add some points from a file, according to the outcome of 
 * the check (1=keep, 0=ignore, -1: stop processing)
 */

static int add_some_points_file (FILE * f,
                                 point_t ** corners_io,
                                 int nc,
                                 int fmt,
                                 int (*check) (void *arg_check,
                                               geom_t * geom),
                                 void *arg_check)
{

  int na = nc;                  // does not harm if na is lower than real na
  point_t *corners = *corners_io;

  while (1) {

    point_t corner;

    int ret = read_point_t (f, &corner, fmt);

    if (ret < 0) {
      *corners_io = NULL;
      delete_points (corners, nc);
      return -1;
    }
    if (ret == 0)
      break;

    int keep = check ? (*check) (arg_check, &corner.geom) : 1;

    if (keep < 0)
      break;

    if (keep) {
      if (nc >= na) {
        na = (na + 1) * 3 / 2;
        corners = realloc (corners, sizeof (point_t) * na);
      }

      corners[nc] = corner;

      nc++;

    } else {
      free (corner.desc);
    }
  }

  *corners_io = corners;

  return nc;


}


int read_points_file (FILE * f, point_t ** corners_out, int fmt)
{
  *corners_out = NULL;
  return add_some_points_file (f, corners_out, 0, fmt, NULL, NULL);
}


int pointset_file_size_and_dim (const char *fname, int fmt, int *d_io)
{
  long int file_size;

  {
    struct stat buf;
    int ret = stat (fname, &buf);
    if (ret < 0) {
      perror ("npt_from_file stat");
      return -1;
    }
    file_size = buf.st_size;
  }

  if (file_size == 0)
    return 0;

  int sz_1, d;

  if (fmt==1 || fmt==2) {
    sz_1 = sizeof (int) + sizeof (geom_t);
    if (fmt == 2)
      sz_1 += sizeof (binsign_t);
    d = 0;
  } else if(fmt==0 || fmt==3) {
    int skip=fmt==0 ? sizeof (geom_t) : 2*sizeof(int)+4*sizeof(float);

    if(d_io && *d_io>=0) {
      d=*d_io;
    } else { /* read dimension from file */
      FILE *f = fopen (fname, "r");
      if (!f) {
        perror ("npt_from_file fopen");
        return -1;
      }
      
      if (fseek (f, skip, SEEK_SET) < 0 ||
          fread (&d, sizeof (int), 1, f) != 1) {
        perror ("npt_from_file read dim");
        return -1;
      }
      fclose (f);
      if (d < 1 || d > 50000) {
        fprintf (stderr, "npt_from_file: weird dimension %d !!", d);
        return -1;
      }
    }
    sz_1 = skip + sizeof (int) + d * sizeof (char);
  } else assert(0);    

  long n = file_size / sz_1;

  if (file_size % sz_1 != 0) {
    fprintf (stderr, "npt_from_file: weird file size %ld d=%d fmt=%d\n",
             file_size, d, fmt);
    return -1;
  }
  
  if(d_io) *d_io=d;

  return n;
}

int count_points (const char *fname, int fmt)
{
  return pointset_file_size_and_dim(fname,fmt,NULL);
}

int read_points_string (const char *string, int string_length,
                        point_t ** corners_out, int fmt)
{
  const unsigned char *rp = (const unsigned char *) string;
  int n = 1, i, d;

  if (string_length == 0) {
    *corners_out = NULL;
    return 0;
  }
  int sz_1;

  if (fmt==1 || fmt==2) {
    sz_1 = sizeof (int) + sizeof (geom_t);
    if (fmt == 2)
      sz_1 += sizeof (binsign_t);
    d = 0;
  } else if(fmt==0 || fmt==3) {
    int skip=fmt==0 ? sizeof (geom_t) : 2*sizeof(int)+4*sizeof(float);    
    memcpy (&d, string + skip, sizeof (int));
    if (d < 1 || d > 50000) {
      fprintf (stderr, "read_points_string: weird dimension %d !!", d);
      return -1;
    }
    sz_1 = sizeof (geom_t) + sizeof (int) + d * sizeof (char);
  } else assert(0);


  n = string_length / sz_1;

  if (string_length % sz_1 != 0) {
    fprintf (stderr, "read_points_string: weird string size %d d=%d fmt=%d",
             string_length, d, fmt);
    return -1;
  }

  point_t *corners = (point_t *) malloc (n * sizeof (point_t));

#define R(v) {memcpy(&(v),rp,sizeof(v)); rp+=sizeof(v); }

  for (i = 0; i < n; i++) {
    corners[i].binsign = 0;

    if (fmt==1 || fmt==2) {
      R (corners[i].vw);
      R (corners[i].geom);
      if (fmt == 2) {
        R (corners[i].binsign)
      } else
        corners[i].binsign = 0;
      corners[i].desc = NULL;
    } else if(fmt==0 || fmt==3) {
      if(fmt==0) {
        R (corners[i].geom);
      } else {
        int ignored;
        R(ignored); R(ignored); 
        geom_t *g=&corners[i].geom;
        R(g->x); R(g->y); R(g->scale); R(g->angle);
        g->mi11=1; g->mi12=0;
        g->mi21=0; g->mi22=1;        
      }
      R (corners[i].dim);
      if (corners[i].dim != d) {
        fprintf (stderr, "read_points_string: non homogeneous size !!");
        return -1;
      }
      corners[i].desc = malloc (d);
      memcpy (corners[i].desc, rp, d);
      corners[i].binsign = 0;
      rp += d;
    }
  }
  *corners_out = corners;

  return n;
#undef R
}



int write_point_t (FILE * f, point_t * corner, int fmt)
{
#define RETERR(no) { \
  perror ("write_point_t error "#no);\
  return -1;\
}\

  if (fmt) {
    if (fwrite (&corner->vw, sizeof (int), 1, f) != 1)
      RETERR (1);
    if (fwrite (&corner->geom, sizeof (geom_t), 1, f) != 1)
      RETERR (2);
    if (fmt == 2) {
      if (fwrite (&corner->binsign, sizeof (binsign_t), 1, f) != 1)
        RETERR (6);
    }
  } else {
    if (fwrite (&corner->geom, sizeof (geom_t), 1, f) != 1)
      RETERR (3);
    if (fwrite (&corner->dim, sizeof (int), 1, f) != 1)
      RETERR (4);
    if (fwrite (corner->desc, 1, corner->dim, f) != corner->dim)
      RETERR (5);
  }
#undef RETERR
  return 1;
}


int write_points (const char *fname, point_t * corners, int nc, int fmt)
{
  FILE *f = fopen (fname, "w");

  if (!f) {
    fprintf (stderr, "write_points: could not open %s\n", fname);
    perror ("");
    return -1;
  }

  int ret = write_points_file (f, corners, nc, fmt);
  fclose (f);
  return ret;
}

int write_points_string (char *str, int len, point_t * corners, int nc,
                         int fmt)
{
  point_t *corner;
  char *write;
  int idx;

  write = str;

  switch (fmt) {
  case 0:
    for (idx = 0, corner = corners; idx < nc; idx++, corner++) {
      memcpy (write, &corner->geom, sizeof (corner->geom));
      write += sizeof (corner->geom);
      memcpy (write, &corner->dim, sizeof (corner->dim));
      write += sizeof (corner->dim);
      memcpy (write, corner->desc, corner->dim * sizeof (char));
      write += corner->dim * sizeof (char);
    }
    break;
  case 1:
    for (idx = 0, corner = corners; idx < nc; idx++, corner++) {
      memcpy (write, &corner->vw, sizeof (corner->vw));
      write += sizeof (corner->vw);
      memcpy (write, &corner->geom, sizeof (corner->geom));
      write += sizeof (corner->geom);
    }
    break;
  case 2:
    for (idx = 0, corner = corners; idx < nc; idx++, corner++) {
      memcpy (write, &corner->vw, sizeof (corner->vw));
      write += sizeof (corner->vw);
      memcpy (write, &corner->geom, sizeof (corner->geom));
      write += sizeof (corner->geom);
      memcpy (write, &corner->binsign, sizeof (corner->binsign));
      write += sizeof (corner->binsign);
    }
    break;
  default:
    fprintf (stderr, "write_points_string: invalid format %d\n", fmt);
    return -1;
  }

  return 0;
}

int write_points_file (FILE * f, point_t * corners, int nc, int fmt)
{
  int i;

  for (i = 0; i < nc; i++) {
    if (write_point_t (f, &corners[i], fmt) < 0) {
      return -1;
    }
  }

  return 0;
}


typedef struct {
  int *mask;
  int n;
} cm_t;


static int check_with_mask (void *arg_check, geom_t * geom)
{
  cm_t *cm = arg_check;
  return cm->mask[cm->n++];
}

int read_points_add_with_mask (const char *fname,
                               point_t ** corners_io,
                               int nc, int fmt, int *mask)
{

  FILE *f = fopen (fname, "r");

  if (!f) {
    fprintf (stderr, "read_points_add_with_mask %s", fname);
    perror ("");
    return -1;
  }

  cm_t cm = { mask, 0 };

  nc = add_some_points_file (f, corners_io, nc, fmt, check_with_mask, &cm);

  fclose (f);

  return nc;
}

int read_points_add (const char *fname,
                     point_t ** corners_io, int nc, int fmt)
{

  FILE *f = fopen (fname, "r");

  if (!f) {
    fprintf (stderr, "read_points_add %s", fname);
    perror ("");
    return -1;
  }

  nc = add_some_points_file (f, corners_io, nc, fmt, NULL, NULL);

  fclose (f);

  return nc;
}



/* Display a set of descriptors */
void display_points (const point_t * des, int n)
{
  int i, j;

  for (i = 0; i < n; i++) {
    const geom_t *g = &des[i].geom;

    fprintf (stdout,
             "x=%.3f y=%.3f sca=%.3f ang=%.3f mi11=%.3f mi12=%.3f mi21=%.3f mi22=%.3f "
             "cornss=%.3f", g->x, g->y, g->scale, g->angle,
             g->mi11, g->mi12, g->mi21, g->mi22, g->cornerness);

    if (des[i].desc == NULL) {
      fprintf (stdout, " vw=%-7d", des[i].vw);
      fprintf (stdout, " binsign=%llx", des[i].binsign);
    } else
      for (j = 0; j < des[i].dim; j++)
        fprintf (stdout, " %d", des[i].desc[j]);

    fprintf (stdout, "\n");
  }
}




static int check_cornerness (void *arg_check, geom_t * geom)
{
  double thresh = *(double *) arg_check;
  return geom->cornerness > thresh;
}







pointset_t *pointset_read (char *fname, int fmt)
{
  vwgeoset_t *ret = pointset_new ();
  ret->n = read_points (fname, &ret->pts, fmt);
  return ret;
}


static int check_n_remain(void *arg_check,geom_t * geom) {
  int *remain=arg_check;
  (*remain)--;
  if(*remain<0) return -1;
  return 1;
}

pointset_t *pointset_read_file_max (FILE *f, int n, int fmt) {
  
  vwgeoset_t *ret = pointset_new ();
  ret->n = add_some_points_file(f,&ret->pts,0,fmt,&check_n_remain,&n);
  return ret;
} 

pointset_t *pointset_read_cornerness (char *fname, double min_cornerness,
                                      int fmt)
{
  pointset_t *ret = pointset_new ();

  FILE *f = fopen (fname, "r");

  if (!f) {
    fprintf (stderr, "pointset_read_cornerness %s", fname);
    perror ("");
    return NULL;
  }

  ret->n =
      add_some_points_file (f, &ret->pts, 0, fmt, check_cornerness,
                            &min_cornerness);

  fclose (f);

  return ret;
}

void pointset_write (char *fname, pointset_t * ps, int fmt)
{
  write_points (fname, ps->pts, ps->n, fmt);
}


/*---------------------------------------------------------------------------*/
/* Sort and filtering functions                                              */
/*---------------------------------------------------------------------------*/

static int vwsgeo_comp_cornerness (const void *vwgeo1, const void *vwgeo2)
{
  point_t *vwg1 = (point_t *) vwgeo1;
  point_t *vwg2 = (point_t *) vwgeo2;

  return vwg2->geom.cornerness - vwg1->geom.cornerness;
}

/* Filter the vwgeo set such that only the n vwgeo with highest cornerness are kept */
void vwgeoset_filter_n_cornerness_max (vwgeoset_t * vgs, int n)
{
  qsort (vgs->pts, vgs->n, sizeof (point_t), vwsgeo_comp_cornerness);

  if (n < vgs->n) {
    point_t *pts = vgs->pts;
    pts = (point_t *) realloc (pts, sizeof (point_t) * n);

    assert (pts != NULL);

    vgs->pts = pts;
    vgs->n = n;
  }
}


static int vwgeo_comp_vw (const void *vwgeo1, const void *vwgeo2)
{
  point_t *vwg1 = (point_t *) vwgeo1;
  point_t *vwg2 = (point_t *) vwgeo2;


  if (vwg1->vw < vwg2->vw)
    return -1;

  if (vwg1->vw > vwg2->vw)
    return 1;

  return 0;
}



/* sort all the vwgeo set such that the indexes of the visual word 
   are in increasing order                                                     */
void vwgeoset_sort (vwgeoset_t * vgs)
{
  qsort (vgs->pts, vgs->n, sizeof (point_t), vwgeo_comp_vw);
}


/* Filter the vwgeo set such that only the points assigned to centroids 
   which are at a similar distance (disratio*dismin) from the nearest centroid */
void vwgeoset_filter_ma_nkeep (vwgeoset_t * vgs, int k, float *vw_centroids_dis,
                               double disratio,int *nkeep)
{
  int i, ii = 0;
  int n = vgs->n / k;

  assert (vgs->n % k == 0);

  point_t *des = vgs->pts;

  for (i = 0; i < n; i++) {
    int j;
    float mindis = 1e9;

    /* find the minimum centroid distance */
    for (j = 0; j < k; j++) {
      if (des[i * k + j].vw < 0)
        break;                  /* not assigned by ann MA */
      if (vw_centroids_dis[i * k + j] < mindis)
        mindis = vw_centroids_dis[i * k + j];
    }

    int ii0=ii;

    for (j = 0; j < k; j++) {
      if (des[i * k + j].vw < 0)
        break;
      if (vw_centroids_dis[i * k + j] < mindis * disratio) {
        des[ii] = des[i * k + j];
        ii++;
      }
    }

    if(nkeep) nkeep[i]=ii-ii0;

  }

  /* At this point, ii contains the number of descriptors -> realloc the correct size */
  vgs->pts = realloc (des, sizeof (point_t) * ii);
  //  printf ("\n%d -> %d | %f\n", vgs->n, ii, ii / (double) n);
  vgs->n = ii;
}


void vwgeoset_filter_thresh_nkeep (vwgeoset_t * vgs, int k, float *vw_centroids_dis,
                                   double thresh,int *nkeep) {
  int i, ii = 0;
  int n = vgs->n / k;

  assert (vgs->n % k == 0);

  point_t *des = vgs->pts;
  
  for (i = 0; i < n; i++) {
    int j;
    int ii0=ii;

    for (j = 0; j < k; j++) {
      if (des[i * k + j].vw < 0)
        break;
      if (vw_centroids_dis[i * k + j] < thresh) {
        des[ii] = des[i * k + j];
        ii++;
      }
    }

    if(nkeep) nkeep[i]=ii-ii0;

  }

  /* At this point, ii contains the number of descriptors -> realloc the correct size */
  vgs->pts = realloc (des, sizeof (point_t) * ii);

  vgs->n = ii;

}


void vwgeoset_filter_ma (vwgeoset_t * vgs, int k, float *vw_centroids_dis,
                         double disratio) {
  vwgeoset_filter_ma_nkeep(vgs,k,vw_centroids_dis,disratio,NULL);
}


static int cmp_cornerness (const void *vwgeo1, const void *vwgeo2)
{
  point_t *vwg1 = (point_t *) vwgeo1;
  point_t *vwg2 = (point_t *) vwgeo2;


  if (vwg1->geom.cornerness > vwg2->geom.cornerness)
    return -1;

  if (vwg1->geom.cornerness < vwg2->geom.cornerness)
    return 1;

  return 0;
}


/* sort by cornerness */
void pointset_sort_by_cornerness (pointset_t * vgs)
{
  qsort (vgs->pts, vgs->n, sizeof (point_t), cmp_cornerness);
}


/* sort according to a permutation */
void pointset_sort_by_permutation (pointset_t * vgs, const int *order)
{
  int i;
  int n = vgs->n;
  int *o = malloc (n * sizeof (*o));

  for (i = 0; i < n; i++)
    o[order[i]] = i;

  for (i = 0; i < n; i++)
    while (o[i] != i) {
      int newpos = o[i];

      point_t pt = vgs->pts[i];
      vgs->pts[i] = vgs->pts[newpos];
      vgs->pts[newpos] = pt;

      int j = o[newpos];
      o[newpos] = o[i];
      o[i] = j;
    }

  free (o);
}


void pointset_crop (pointset_t * ps, float xmin, float ymin, float xmax,
                    float ymax, int keep_outside)
{
  int i, wp = 0;

  keep_outside = !!keep_outside;        /* want 0 or 1 */
  for (i = 0; i < ps->n; i++) {
    point_t *pt = &ps->pts[i];
    if ((xmin < pt->geom.x && pt->geom.x < xmax &&
         ymin < pt->geom.y && pt->geom.y < ymax) ^ keep_outside)
      ps->pts[wp++] = ps->pts[i];
    else
      free (pt->desc);
  }
  ps->n = wp;
}



void pointset_crop_polygon (pointset_t * ps, double *coeffs, int n)
{
  int i, wp = 0, j;

  for (i = 0; i < ps->n; i++) {
    point_t *pt = &ps->pts[i];
    for (j = 0; j < n; j++) {
      double a = coeffs[3 * j], b = coeffs[3 * j + 1], c = coeffs[3 * j + 2];
      if (!(pt->geom.x * a + pt->geom.y * b + c > 0))
        break;
    }
    if (j == n)
      ps->pts[wp++] = ps->pts[i];
    else
      free (pt->desc);
  }
  ps->n = wp;
}


void pointset_crop_n (pointset_t * ps, int n)
{
  int i;
  for (i = n; i < ps->n; i++)
    free (ps->pts[i].desc);
  ps->n = n;
}


void vwgeoset_filter_duplicate_vw (vwgeoset_t * vgs)
{
  int i, ii = 1;
  point_t *des = vgs->pts;

  if (vgs->n == 0)
    return;

  vwgeoset_sort (vgs);

  for (i = 1; i < vgs->n; i++) {

    if (des[i].vw == des[i - 1].vw)
      continue;

    des[ii++] = des[i];
  }

  pointset_crop_n (vgs, ii);
}


void pointset_filter_random (pointset_t * ps, float v)
{
  int i, ii = 1;

  if (ps->n == 0)
    return;

  for (i = 1; i < ps->n ; i++) {

    if (drand48 () < v) {       /* Keep this descriptor ? */
      if (i > ii)
        ps->pts[ii] = ps->pts[i];
      ii++;
    }
    else
      free (ps->pts[i].desc);
  }
  ps->n = ii;
}


void pointset_filter_cornerness (pointset_t * ps, float cornerness_min)
{
  int i, j;
  j = 0;
  for (i = 0; i < ps->n; i++) {
    if (ps->pts[i].geom.cornerness >= cornerness_min)
      ps->pts[j++] = ps->pts[i];
    else
      free (ps->pts[i].desc);
  }
  ps->n = j;
}


void vwgeoset_mask_binsign (pointset_t * ps, binsign_t mask) {
  int i;

  for (i = 0; i < ps->n; i++) 
    ps->pts[i].binsign &= mask;

}
