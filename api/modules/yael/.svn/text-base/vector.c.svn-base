#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <math.h>

#include "vector.h"


#ifdef __linux__
#include <malloc.h>
#else
static void *memalign(size_t ignored,size_t nbytes) {
  return malloc(nbytes);
} 
#endif


float *fvec_new (long n)
{
  float *ret = (float *) memalign (16, sizeof (*ret) * n);
  if (!ret) {
    fprintf (stderr, "fvec_new %ld : out of memory\n", n);
    abort();
  }
  return ret;
}

int *ivec_new (long n)
{
  int *ret = (int *) malloc (sizeof (*ret) * n);
  if (!ret) {
    fprintf (stderr, "ivec_new %ld : out of memory\n", n);
    abort();
  }
  return ret;
}

float *fvec_new_0 (long n)
{
  float *ret = (float *) calloc (sizeof (*ret), n);
  if (!ret) {
    fprintf (stderr, "ivec_new_0 %ld : out of memory\n", n);
    abort();
  }
  return ret;
}

float *fvec_new_rand(long n) {
  float *f=fvec_new (n);
  long i;
  for(i=0;i<n;i++) f[i]=drand48();
  return f;
}



int *ivec_new_0 (long n)
{
  int *ret = (int *) calloc (sizeof (*ret), n);
  if (!ret) {
    fprintf (stderr, "ivec_new_0 %ld : out of memory\n", n);
    abort();
  }
  return ret;
}


float *fvec_new_set (long n, float val)
{
  int i;
  float *ret = (float *) calloc (sizeof (*ret), n);
  if (!ret) {
    fprintf (stderr, "fvec_new_set %ld : out of memory\n", n);
    abort();
  }

  for (i = 0 ; i < n ; i++)
    ret[i] = val;

  return ret;
}


int *ivec_new_set (long n, int val)
{
  int i;
  int *ret = (int *) calloc (sizeof (*ret), n);
  if (!ret) {
    fprintf (stderr, "ivec_new_set %ld : out of memory\n", n);
    abort();
  }

  for (i = 0 ; i < n ; i++)
    ret[i] = val;

  return ret;
}


int * ivec_new_range (long a, long b)
{
  int i;
  int *ret = (int *) calloc (sizeof (*ret), b - a);
  if (!ret) {
    fprintf (stderr, "ivec_new_range : out of memory\n");
    abort();
  }

  for (i = a ; i < b ; i++)
    ret[i - a] = i;

  return ret;
}


int * ivec_new_copy (const int * v, long n)
{
  int *ret = (int *) malloc (sizeof (int) * n);
  if (!ret) {
    fprintf (stderr, "ivec_new %ld : out of memory\n", n);
    abort();
  }
  
  memcpy (ret, v, n * sizeof (*ret));
  return ret;
}


float * fvec_new_copy (const float * v, long n) {
  float *ret=fvec_new(n);
  memcpy (ret, v, n * sizeof (*ret));
  return ret;  
}

int *ivec_new_random_perm (long n)
{

  int *perm = ivec_new (n);
  int i;

  for (i = 0; i < n; i++)
    perm[i] = i;

  for (i = 0; i < n - 1; i++) {
    int j = i + random () % (n - i);
    /* swap i and j */
    int p = perm[i];
    perm[i] = perm[j];
    perm[j] = p;
  }

  return perm;
}



int *ivec_new_histogram (int k, int *v, int n)
{
  int i;
  int *h = ivec_new_0 (k);

  for (i = 0; i < n; i++) {
    assert (v[i] >= 0 && v[i] < k);
    h[v[i]]++;
  }

  return h;
}

void fvec_splat_add(const float *a,int n,
                    const int *assign,float *accu) {
  int i;
  for(i=0;i<n;i++) 
    accu[assign[i]] += a[i];
}


void fvec_isplat_add(const float *a,int n,
                     const int *assign,float *accu) {
  int i;
  for(i=0;i<n;i++) 
    accu[i] += a[assign[i]];
  
}

int ivec_hash(const int * v, long n) {
  unsigned int *k=(unsigned int*)v;
  unsigned int b    = 378551;
  unsigned int a    = 63689;
  unsigned int hash = 0;
  int i;
  for(i = 0; i < n; i++) {
    hash = hash * a + k[i];
    a    = a * b;
  }
  return hash;
}


int ivec_count_occurrences(const int * v, long n, int val) {
  int count=0;
  while(n--) if(v[n]==val) count++;
  return count;
}

void ivec_accumulate_slices(const int *v,int *sl,int n) {
  int i;
  int accu=0,j=0;
  for(i=0;i<n;i++) {
    while(j<sl[i]) 
      accu+=v[j++];
    sl[i]=accu;
  }
}



/*---------------------------------------------------------------------------*/
/* Input/Output functions                                                    */
/*---------------------------------------------------------------------------*/


int fvec_fwrite(FILE *fo, const float *v, int d) 
{
  int ret;
  ret = fwrite (&d, sizeof (int), 1, fo);
  if (ret != 1) {
    perror ("fvec_fwrite: write error 1");
    return -1;
  }
  ret = fwrite (v, sizeof (float), d, fo);
  if (ret != d) {
    perror ("fvec_fwrite: write error 2");
    return -1;
  }  
  return 0;
}

int fvec_fwrite_raw(FILE *fo, const float *v, long d) {
  long ret = fwrite (v, sizeof (float), d, fo);
  if (ret != d) {
    perror ("fvec_fwrite: write error 2");
    return -1;
  }  
  return 0;
}


int fvecs_fwrite (FILE *fo, int d, int n, const float *vf)
{
  int i;
  /* write down the vectors as fvecs */
  for (i = 0; i < n; i++) {
    if(fvec_fwrite(fo, vf+i*d, d)<0) 
      return i;
  }
  return n;
}



int fvecs_write (const char *fname, int d, int n, const float *vf)
{
  FILE *fo = fopen (fname, "w");
  if (!fo) {
    perror ("fvecs_write: cannot open file");
    return -1;
  }

  int ret = fvecs_fwrite (fo, d, n, vf);

  fclose (fo);
  return ret;
}




int fvecs_new_read (const char *fname, int *d_out, float **vf_out)
{

  FILE *f = fopen (fname, "r");
  if (!f) {
    fprintf (stderr, "fvecs_new_read: could not open %s\n", fname);
    perror ("");
    return -1;
  }

  int n=fvecs_new_fread_max(f,d_out,vf_out,-1);

  fclose (f);
  
  return n;
}

int fvecs_new_fread_max (FILE *f, int *d_out, float **vf_out, long nmax)
{
  int n, d = -1;

  float *vf=NULL;
  long nalloc=0;
  
  for(n=0; nmax<0 || n<nmax ; n++) {
    int new_d;

    if (fread (&new_d, sizeof (int), 1, f) != 1) {
      if (feof (f))
        break;
      else {
        perror ("fvecs_new_read error 1");
        // TODO free mem
        return -1;
      }
    }

    if (n == 0)
      d = new_d;
    else if (new_d != d) {
      fprintf (stderr, "fvecs_new_read error 2");
      return -1;
    }

    if((n+1)*d>nalloc) {
      nalloc=nalloc==0 ? 1024 : nalloc*3/2;
      vf=realloc(vf,nalloc*sizeof(float));
    }

    if (fread (vf+n*d, sizeof (float), d, f) != d) {
      fprintf (stderr, "fvecs_new_read error 3");
      return -1;
    }

  }

  *vf_out = vf;
  *d_out = d;
  return n;
}

int fvecs_new_read_sparse (const char *fname, int d, float **vf_out) {
  float *vf=NULL;
  long n=0,na=0;
  float *vals=fvec_new(d);
  int *idx=ivec_new(d);
  
  FILE *f = fopen (fname, "r");
#define E(msg) {                                                \
  fprintf (stderr, "fvecs_new_read_sparse %s: " msg , fname);   \
  perror ("");                                                  \
  free(vf); free(vals); free(idx);                              \
  return -1;                                                    \
}
  if (!f) E("");
  
  while(!feof(f)) {
    int nz,ret,nz2;
    ret=fread(&nz,sizeof(int),1,f);
    if(ret!=1) {
      if(feof(f)) break;
      E("err 1");
    }
    if(fread(idx,sizeof(int),nz,f)!=nz) E("err 2");
    if(fread(&nz2,sizeof(int),1,f)!=1) E("err 3");
    if(nz!=nz2) E("err 4");
    if(fread(vals,sizeof(float),nz,f)!=nz) E("err 5");
    
    if(n>=na) {
      na=(na+1)*3/2;
      vf=realloc(vf,na*sizeof(float)*d);
    }
    
    float *dense=spfvec_to_fvec (idx,vals,nz,d);
    memcpy(vf+n*d,dense,sizeof(float)*d);
    free(dense);
    
    n++;       
  }
#undef E
  free(vals);
  free(idx);
  fclose(f);
  *vf_out=vf;
  return n;
}




int fvecs_read (const char *fname, int d, int n, float *a)
{
  FILE *f = fopen (fname, "r");
  if (!f) {
    fprintf (stderr, "fvecs_read: could not open %s\n", fname);
    perror ("");
    return -1;
  }

  int i;
  for (i = 0; i < n; i++) {
    int new_d;

    if (fread (&new_d, sizeof (int), 1, f) != 1) {
      if (feof (f))
        break;
      else {
        perror ("fvecs_read error 1");
        // TODO free mem
        return -1;
      }
    }

    if (new_d != d) {
      fprintf (stderr, "fvecs_read error 2");
      return -1;
    }

    if (fread (a + d * (long) i, sizeof (float), d, f) != d) {
      fprintf (stderr, "fvecs_read error 3");
      return -1;
    }
    n++;
  }
  fclose (f);

  return i;
}



int fvec_read (const char *fname, int d, float *a, int o_f)
{
  FILE *f = fopen (fname, "r");
  if (!f) {
    fprintf (stderr, "fvec_read: could not open %s\n", fname);
    perror ("");
    return -1;
  }

  if (fseek (f, (sizeof (int) + sizeof (float) * d) * o_f, SEEK_SET) != 0) {
    fprintf (stderr, "fvec_read %s fseek", fname);
    perror ("");
    return -1;
  }

  int d2 = 0;
  int ret = fread (&d2, sizeof (int), 1, f);
  if (d != d2) {
    fprintf (stderr, "fvec_read %s bad dim: %d!=%d\n", fname, d, d2);
    return -1;
  }

  ret = fread (a, sizeof (float), d, f);
  fclose (f);
  return ret;
}



int fvec_fread (FILE * f, float * v)
{
  int d;
  int ret = fread (&d, sizeof (int), 1, f);

  if (ret != 1) {
    perror ("# fvec_fread error 1");
    return -1;
  }

  ret = fread (v, sizeof (*v), d, f);
  if (ret != d) {
    perror ("# fvec_fread error 2");
    return -1;
  }

  return d;
}


float *fvec_fread_raw(FILE * f, long d) {
  float *v=fvec_new(d);

  int ret = fread (v, sizeof (*v), d, f);
  if (ret != d) {
    free(v);
    perror ("# fvec_fread error 2");
    return NULL;
  }
  return v;
}

int ivec_fwrite (FILE *f, const int *v, int d)
{
  int ret = fwrite (&d, sizeof (d), 1, f);
  if (ret != 1) {
    perror ("ivec_fwrite: write error 1");
    return -1;
  }

  ret = fwrite (v, sizeof (*v), d, f);
  if (ret != d) {
    perror ("ivec_fwrite: write error 2");
    return -2;
  }
  return 0;
}


int ivecs_fwrite(FILE *f, int d, int n, const int *v)
{
  int i;
  for (i = 0 ; i < n ; i++)
    ivec_fwrite (f, v, d);
  return n;
}


int ivecs_write (const char *fname, int d, int n, const int *v)
{
  int ret = 0;
  FILE *f = fopen (fname, "w");
  if (!f) {
    perror ("ivecs_write");
    return -1;
  }

  ret = ivecs_fwrite (f, d, n, v);

  fclose (f);
  return ret;
}


int ivecs_new_read (const char *fname, int *d_out, int **vi_out)
{

  FILE *f = fopen (fname, "r");
  if (!f) {
    fprintf (stderr, "fvecs_new_read: could not open %s\n", fname);
    perror ("");
    return -1;
  }

  int n, d = -1;

  struct buf_l {
    struct buf_l *next;
    int l[0];
  };
  struct buf_l *bufs = NULL;

  n = 0;

  while (1) {
    int new_d;

    if (fread (&new_d, sizeof (int), 1, f) != 1) {
      if (feof (f))
        break;
      else {
        perror ("load_clusters_fvecs error 1");
        // TODO free mem
        return -1;
      }
    }

    if (n == 0)
      d = new_d;
    else if (new_d != d) {
      fprintf (stderr, "load_clusters_fvecs error 2");
      return -1;
    }

    {
      struct buf_l *new = malloc (sizeof (struct buf_l) + sizeof (int) * d);
      new->next = bufs;
      bufs = new;
    }
    if (fread (bufs->l, sizeof (int), d, f) != d) {
      fprintf (stderr, "load_clusters_fvecs error 3");
      return -1;
    }
    n++;
  }
  fclose (f);

  int * vi = ivec_new (n * d);

  int i = n;
  while (i--) {
    memcpy (vi + i * d, bufs->l, sizeof (int) * d);
    struct buf_l *tmp = bufs;
    bufs = bufs->next;
    free (tmp);
  }
  *vi_out = vi;
  *d_out = d;
  return n;
}


int *ivec_new_read(const char *fname, int *d_out) {
  int d, n;
  int *vi;
  n = ivecs_new_read(fname,&d,&vi);
  if (n<0) 
    return NULL;
  assert (n==1);
  if (d_out) 
    *d_out=d;
  return vi;
}


int ivec_fread (FILE * f, int * v)
{
  int d;
  int ret = fread (&d, sizeof (int), 1, f);

  if (ret != 1) {
    perror ("# ivec_fread error 1");
    return -1;
  }

  ret = fread (v, sizeof (*v), d, f);
  if (ret != d) {
    perror ("# ivec_fread error 2");
    return -1;
  }
  return d;
}



void fvec_print (const float * v, int n)
{
  int i;
  printf ("[");
  for (i = 0 ; i < n ; i++)
    printf ("%f ", v[i]);
  printf ("]\n");
}


void ivec_print (const int * v, int n)
{
  int i;
  printf ("[");
  for (i = 0 ; i < n ; i++)
    printf ("%d ", v[i]);
  printf ("]\n");
}


void fvec_fprintf (FILE * f, const float *v, int n, const char *fmt)
{
  int i;
  if (fmt == NULL)
    fmt = "%f ";

  for (i = 0; i < n; i++)
    fprintf (f, fmt, v[i]);
}


void ivec_fprintf (FILE * f, const int *v, int n, const char *fmt)
{
  int i;
  if (fmt == NULL)
    fmt = "%d ";

  for (i = 0; i < n; i++)
    fprintf (f, fmt, v[i]);
}




/*---------------------------------------------------------------------------*/
/* Vector manipulation                                                       */
/*---------------------------------------------------------------------------*/

void fvec_0(float * v, int n)
{
  memset (v, 0, n * sizeof (*v));
}


void ivec_0(int * v, int n)
{
  memset (v, 0, n * sizeof (*v));
}


void fvec_set (float * v, int n, float val)
{
  int i;
  for (i = 0 ; i < n ; i++)
    v[i] = val;
}


void ivec_set (int * v, int n, int val)
{
  int i;
  for (i = 0 ; i < n ; i++)
    v[i] = val;
}


void fvec_incr (float * v, int n, double scal)
{
  int i = 0;
  for (i = 0 ; i < n ; i++)
    v[i] += scal;
}


void fvec_decr (float * v, int n, double scal)
{
  int i = 0;
  for (i = 0 ; i < n ; i++)
    v[i] -= scal;
}


void fvec_mul_by (float * v, int n, double scal)
{
  int i = 0;
  for (i = 0 ; i < n ; i++)
    v[i] *= scal;
}


void fvec_div_by (float * v, int n, double scal)
{
  int i = 0;
  for (i = 0 ; i < n ; i++)
    v[i] /= scal;
}


void fvec_add (float * v1, const float * v2, int n)
{
  int i = 0;
  for (i = 0 ; i < n ; i++)
    v1[i] += v2[i];
}


void fvec_sub (float * v1, const float * v2, int n)
{
  int i = 0;
  for (i = 0 ; i < n ; i++)
    v1[i] -= v2[i];
}


void fvec_mul (float * v1, const float * v2, int n)
{
  int i = 0;
  for (i = 0 ; i < n ; i++)
    v1[i] *= v2[i];
}


void fvec_div (float * v1, const float * v2, int n)
{
  int i = 0;
  for (i = 0 ; i < n ; i++)
    v1[i] /= v2[i];
}


void fvec_normalize (float * v, int n, double norm)
{
  if(norm==0) return;

  double nr = fvec_norm (v, n, norm);

  /*  if(nr!=0)*/
  fvec_mul_by (v, n, 1 / nr);
  
}

int fvec_purge_nans(float * v, long n, float replace_value) {
  long i,count=0;
  
  for(i=0;i<n;i++) if(isnan(v[i])) {
    count++;
    v[i]=replace_value;
  }
  
  return count;
}


void fvec_sqrt (float * v, int n)
{
  int i;
  for (i = 0 ; i < n ; i++)
    v[i] = sqrt (v[i]);
}


void fvec_sqr (float * v, int n)
{
  int i;
  for (i = 0 ; i < n ; i++)
    v[i] =  v[i] * v[i];
}


/*---------------------------------------------------------------------------*/
/* Vector measures and statistics                                            */
/*---------------------------------------------------------------------------*/

double fvec_sum (const float * v, int n)
{
  int i;
  double s = 0;
  for (i = 0 ; i < n ; i++)
    s += v[i];

  return s;
}


long long ivec_sum (const int * v, int n)
{
  int i;
  long long s = 0;
  for (i = 0 ; i < n ; i++)
    s += v[i];

  return s;
}

long long ivec_sum_2 (const int * v, int n)
{
  int i;
  long long s = 0;
  for (i = 0 ; i < n ; i++)
    s += v[i]*(long long)v[i];

  return s;
}


double fvec_norm (const float * v, int n, double norm)
{

  if(norm==0) return n;

  int i;
  double s = 0;

  if(norm==1) {
    for (i = 0 ; i < n ; i++)
      s += fabs(v[i]);
    return s;
  }

  if(norm==2) {
    for (i = 0 ; i < n ; i++)
      s += v[i]*v[i];
    return sqrt(s);
  }

  for (i = 0 ; i < n ; i++)
    s += pow (v[i], norm);

  return pow (s, 1 / norm);
}


double fvec_norm2sqr (const float * v, int n) {
  double s=0;
  int i;
  for (i = 0 ; i < n ; i++)
    s += v[i]*v[i];
  return s;
}

int fvec_nz (const float * v, int n)
{
  int i, nz = 0;
  for (i = 0 ; i < n ; i++)
    if (v[i] != 0)
      nz++;

  return nz;
}


int ivec_nz (const int * v, int n)
{
  int i, nz = 0;
  for (i = 0 ; i < n ; i++)
    if (v[i] != 0)
      nz++;

  return nz;
}


int fvec_find (const float *v, int n, int ** nzpos_out)
{
  int nz = fvec_nz (v, n);
  int * nzpos = ivec_new (nz);
  int i, ii = 0;

  for (i = 0 ; i < n ; i++) 
    if (v[i] != 0) {
      nzpos[ii] = i;
      ii++;
    }

  *nzpos_out = nzpos;
  return nz;
}


int ivec_find (const int *v, int n, int ** nzpos_out)
{
  int nz = ivec_nz (v, n);
  int * nzpos = ivec_new (nz);
  int i, ii = 0;

  for (i = 0 ; i < n ; i++) 
    if (v[i] != 0) {
      nzpos[ii] = i;
      ii++;
    }

  *nzpos_out = nzpos;
  return nz;
}


void ivec_shuffle (int * v, long n)
{
  int i;

  for (i = 0; i < n - 1; i++) {
    int j = i + random () % (n - i);
    /* swap i and j */
    int p = v[i];
    v[i] = v[j];
    v[j] = p;
  }
}


double entropy (const float *pmf, int n)
{
  int i;
  double minusent = 0;

  for (i = 0 ; i < n ; i++)
    if (pmf[i] > 0)
      minusent += pmf[i] * log (pmf[i]);

  return - minusent / log(2);
}

/*! @brief entropy of a binary variable */
double binary_entropy (double p)
{
  if (p == 0 || p == 1)
    return 0;
  return -(p * log (p) + (1-p) * log (1-p)) / log (2);
}


/*---------------------------------------------------------------------------*/
/* Distances                                                                 */
/*---------------------------------------------------------------------------*/

int ivec_distance_hamming (const int * v1, const int * v2, int n)
{
  int i, dis = 0; 

  for (i = 0 ; i < n ; i++)
    dis += (v1[i] == v2[i] ? 0 : 1);

  return dis;
}


double fvec_distance_L2 (const float * v1, const float * v2, int n)
{
  return sqrt (fvec_distance_L2sqr (v1, v2, n));
}


double fvec_distance_L1 (const float * v1, const float * v2, int n) {
  int i;
  double dis = 0;

  for (i = 0 ; i < n ; i++) {
    dis += fabs(v1[i] - v2[i]);
  }

  return dis;  
}

double fvec_distance_L2sqr (const float * v1, const float * v2, int n)
{
  int i;
  double dis = 0, a;

  for (i = 0 ; i < n ; i++) {
    a = (double) v1[i] - v2[i];
    dis += a * a;
  }

  return dis;  
}


double fvec_inner_product (const float * v1, const float * v2, int n)
{
  int i;
  double res = 0;
  for (i = 0 ; i < n ; i++)
    res += v1[i] * v2[i];
  return res;
}


/*---------------------------------------------------------------------------*/
/* Sparse vector handling                                                    */
/*---------------------------------------------------------------------------*/

int fvec_to_spfvec (float * v, int n, int ** idx_out, float ** v_out)
{
  int i, ii = 0;
  int nz = fvec_nz (v, n);
  int * idx = ivec_new (nz);
  float * val = fvec_new (nz);

  for (i = 0 ; i < n ; i++) 
    if (v[i] != 0) {
      idx[ii] = i;
      val[ii] = v[i];
      ii++;
    }

  *idx_out = idx;
  *v_out = val;
  return nz;
}


int ivec_to_spivec (int * v, int n, int ** idx_out, int ** v_out)
{
  int i, ii = 0;
  int nz = ivec_nz (v, n);
  int * idx = ivec_new (nz);
  int * val = ivec_new (nz);

  for (i = 0 ; i < n ; i++) 
    if (v[i] != 0) {
      idx[ii] = i;
      val[ii] = v[i];
      ii++;
    }

  *idx_out = idx;
  *v_out = val;
  return nz;
}


float * spfvec_to_fvec (int * idx, float * v, int nz, int n)
{
  int i;
  float * ret = fvec_new_0 (n);
  for (i = 0 ; i < nz ; i++)
    ret[idx[i]] = v[i];

  return ret;
}


int * spivec_to_ivec (int * idx, int * v, int nz, int n)
{
  int i;
  int * ret = ivec_new_0 (n);
  for (i = 0 ; i < nz ; i++)
    ret[idx[i]] = v[i];

  return ret;
}


float spfvec_inner_product (int *idx1, float *val1, int nz1, 
			    int *idx2, float *val2, int nz2)
{
  double s = 0;
  long i1 = 0, i2 = 0;

  while (1) {

    if (i1 == nz1) 
      break;

    if (i2 == nz2) 
      break;

    if (idx1[i1] == idx2[i2]) {
      s += val1[i1] * val2[i2];
      i1++;
      i2++;
    }

    else {
      if (idx1[i1] < idx2[i2])
	i1++;

      else
	i2++;
    }
  }
  return s;
}


long ivec_index(const int * v, long n,int val) {
  long i;
  for(i=0;i<n;i++) if(v[i]==val) return i;
  return -1;
}

