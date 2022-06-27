#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>

#include "vector.h"
#include "matrix.h"


#define NEWA(type,n) (type*)malloc(sizeof(type)*(n))
#define NEWAC(type,n) (type*)calloc(sizeof(type),(n))
#define NEW(type) NEWA(type,1)


/* blas/lapack subroutines */

#define real float
#define integer int

int sgemm_ (char *transa, char *transb, integer * m, integer *
            n, integer * k, real * alpha, const real * a, integer * lda,
            const real * b, integer * ldb, real * beta, real * c__,
            integer * ldc);

int ssyev_ (char *jobz, char *uplo, integer * n, real * a,
            integer * lda, real * w, real * work, integer * lwork,
            integer * info);


int sgeqrf_ (integer * m, integer * n, real * a, integer * lda,
             real * tau, real * work, integer * lwork, integer * info);

int slarft_ (char *direct, char *storev, integer * n, integer *
             k, real * v, integer * ldv, real * tau, real * t, integer * ldt);

int slarfb_ (char *side, char *trans, char *direct, char *storev, integer * m,
             integer * n, integer * k, real * v, integer * ldv, real * t,
             integer * ldt, real * c__, integer * ldc, real * work,
             integer * ldwork);

#undef real
#undef integer


/*---------------------------------------------------------------------------*/
/* Standard operations                                                       */
/*---------------------------------------------------------------------------*/

/* Generate Gaussian random value, mean 0, variance 1  
   From Python source code. */

#define NV_MAGICCONST  1.71552776992141

static double gaussrand ()
{
  double z;
  while (1) {
    float u1, u2, zz;
    u1 = drand48 ();
    u2 = drand48 ();
    z = NV_MAGICCONST * (u1 - .5) / u2;
    zz = z * z / 4.0;
    if (zz < -log (u2))
      break;
  }
  return z;
}



/*---------------------------------------------------------------------------*/
/* Standard operations                                                       */
/*---------------------------------------------------------------------------*/

float *fmat_new (int nrow, int ncol)
{
  float *m = fvec_new (nrow * ncol);
  return m;
}


void fmat_mul (const float *left, const float *right,
               int n, int m, int k, float *mout)
{

  float alpha = 1;
  float beta = 0;
  char trans = 'n';

  sgemm_ (&trans, &trans, &k, &n, &m,
          &alpha, right, &k, left, &m, &beta, mout, &k);
}


float *fmat_new_mul (const float *left, const float *right, int n, int m,
                     int k)
{
  float *a = fmat_new (k, n);
  fmat_mul (left, right, n, m, k, a);
  return a;
}


void fmat_mul_tl (const float *left, const float *right,
                  int n, int m, int k, float *mout)
{

  float alpha = 1;
  float beta = 0;
  char transleft = 't';
  char transright = 'n';

  sgemm_ (&transright, &transleft, &k, &n, &m,
          &alpha, right, &k, left, &n, &beta, mout, &k);
}


float *fmat_new_mul_tl (const float *left, const float *right,
                        int n, int m, int k)
{
  float *a = fmat_new (k, n);
  fmat_mul_tl (left, right, n, m, k, a);
  return a;
}


void fmat_mul_tr (const float *left, const float *right,
                  int n, int m, int k, float *mout)
{

  float alpha = 1;
  float beta = 0;
  char transleft = 'n';
  char transright = 't';

  sgemm_ (&transright, &transleft, &k, &n, &m,
          &alpha, right, &m, left, &m, &beta, mout, &k);
}


float *fmat_new_mul_tr (const float *left, const float *right,
                        int n, int m, int k)
{
  float *a = fmat_new (k, n);
  fmat_mul_tr (left, right, n, m, k, a);
  return a;
}


void fmat_mul_tlr (const float *left, const float *right,
                   int n, int m, int k, float *mout)
{

  float alpha = 1;
  float beta = 0;
  char transleft = 't';
  char transright = 't';

  sgemm_ (&transright, &transleft, &k, &n, &m,
          &alpha, right, &m, left, &n, &beta, mout, &k);
}


float *fmat_new_mul_tlr (const float *left, const float *right,
                         int n, int m, int k)
{
  float *a = fmat_new (k, n);
  fmat_mul_tlr (left, right, n, m, k, a);
  return a;
}


/*! @brief Multiply a matrix by a vector */
float * fmat_mul_fvec (const float * a, const float * v, int nrow, int ncol)
{
  int i;
  float * res = malloc (nrow * sizeof (*res));
  for (i = 0 ; i < nrow ; i++)
    res[i] = fvec_inner_product (a + i * ncol, v, ncol);

  return res;
}


void fmat_display (const float *a, int nrow, int ncol)
{
  int i, j;

  printf ("[");
  for (i = 0; i < nrow; i++) {
    for (j = 0; j < ncol; j++)
      printf ("%.5g ", a[i * ncol + j]);
    if (i == nrow - 1)
      printf ("]\n");
    else
      printf (";\n");
  }
}


/*---------------------------------------------------------------------------*/
/* Matrix manipulation functions                                             */
/*---------------------------------------------------------------------------*/

float *fmat_get_submatrix (const float *a, int ncola, int r1, int c1, int r2,
                           int c2)
{
  int i, j;
  int nrow = r2 - r1;
  int ncol = c2 - c1;
  float *b = fmat_new (nrow, ncol);

  for (i = r1; i < r2; i++)
    for (j = c1; j < c2; j++)
      b[(i - r1) * ncol + (j - c1)] = a[i * ncola + j];

  return b;
}

float *fmat_get_rows (const float *a, int ncol, int nrowout, const int *rows)
{
  int i;
  float *b = fmat_new (nrowout, ncol);

  for (i = 0; i < nrowout; i++)
    memcpy (b + i * ncol, a + rows[i] * ncol, sizeof (*a) * ncol);

  return b;
}


/*---------------------------------------------------------------------------*/
/* Special matrices                                                          */
/*---------------------------------------------------------------------------*/
float *fmat_new_rand_gauss (int nrow, int ncol)
{
  int i;
  float *m = fmat_new (nrow, ncol);

  for (i = 0; i < nrow * ncol; i++)
    m[i] = gaussrand ();

  return m;
}


/* method: we compute the QR decomposition of a matrix with Gaussian
   values */
float *random_orthogonal_basis (int d)
{
  int i;


  /* generate a Gaussian matrix */
  float *x = fmat_new_rand_gauss (d, d);

  float *tau = NEWA (float, d);

  {                             /* compute QR decomposition */

    /* query work size */
    float lwork_query;
    int lwork = -1;
    int info;
    sgeqrf_ (&d, &d, x, &d, tau, &lwork_query, &lwork, &info);
    assert (info == 0);

    lwork = (int) lwork_query;
    float *work = NEWA (float, lwork);
    sgeqrf_ (&d, &d, x, &d, tau, work, &lwork, &info);
    assert (info == 0);

    free (work);
  }

  /* Decomposition now stored in x and tau. Apply to identity to get
     explicit matrix Q */

  float *q = NEWAC (float, d * d);
  {

    float *t = NEWA (float, d * d);

    slarft_ ("F", "C", &d, &d, x, &d, tau, t, &d);

    for (i = 0; i < d; i++)
      q[i + d * i] = 1;

    float *work = NEWA (float, d * d);

    slarfb_ ("Left", "N", "F", "C",
             &d, &d, &d, x, &d, t, &d, q, &d, work, &d);

    free (t);
    free (work);
  }

  free (tau);
  free (x);
  return q;
}


/* Construct a Hadamard matrix of dimension d using the Sylvester construction.
   d should be a power of 2 */
float *hadamard (int d)
{
  assert ((d & (d - 1)) == 0 || !"d must be power of 2");

  int i, j;
  float *had = fvec_new (d * d);

  if (d == 1) {
    had[0] = 1;
    return had;
  }

  /* Compute the Hadamard matrix of dimension d / 2 */
  int dd = d / 2;
  float *had_part = hadamard (dd);

  for (i = 0; i < dd; i++)
    for (j = 0; j < dd; j++) {
      had[i * d + j] = had_part[i * dd + j];
      had[i * d + j + dd] = had_part[i * dd + j];
      had[(i + dd) * d + j] = had_part[i * dd + j];
      had[(i + dd) * d + j + dd] = -had_part[i * dd + j];
    }

  free (had_part);
  return (had);
}




/*---------------------------------------------------------------------------*/
/* Statistical matrix operations                                             */
/*---------------------------------------------------------------------------*/


/* Input matrix: v(n,d) stored by rows.

   x is v data centered for each dimension 0<=j<n
   x = v - (1/n) * u * m 

   where :
   *   u(n,1) contains only 1's
   *   m = u' * v is the sum of values for each column of v 

   cov is the covariance matrix :

   cov = (1/n) x' * x
       = (1/n) v' * v - (1/n^2) * m' * m

   => no need to allocate an auxiliary array.
*/


float *compute_pca (int n, int d, float *v)
{

  double *sums = NEWAC (double, d);
  int i, j;

  for (i = 0; i < n; i++)
    for (j = 0; j < d; j++)
      sums[j] += v[i * d + j];

  float *cov = fvec_new (d * d);

  for (i = 0; i < d; i++)
    for (j = 0; j < d; j++)
      cov[i + j * d] = sums[i] * sums[j];

  free (sums);

  {
    float alpha = 1.0 / n, beta = -1.0 / (n * n);
    sgemm_ ("N", "T", &d, &d, &n, &alpha, v, &d, v, &d, &beta, cov, &d);
  }

  {
    int lwork = -1;
    float optimal_lwork;
    int info;

    /* query work size */
    ssyev_ ("V", "U", &d, NULL, &d, NULL, &optimal_lwork, &lwork, &info);
    assert (info == 0);
    lwork = (int) optimal_lwork;

    float *work = NEWA (float, lwork);
    float *eigenvals = NEWA (float, d);

    ssyev_ ("V", "U", &d, cov, &d, eigenvals, work, &lwork, &info);

    free (work);
    free (eigenvals);

    if (info != 0) {
      free (cov);
      return NULL;
    }

  }

  /* revert order of vectors to get vectors corresponding to the
     biggest eigenvalues first */

  for (i = 0; i < d / 2; i++) {
    int i2 = d - 1 - i;
    for (j = 0; j < d; j++) {
      float tmp = cov[i * d + j];
      cov[i * d + j] = cov[i2 * d + j];
      cov[i2 * d + j] = tmp;
    }
  }

  return cov;
}




void fmat_splat_separable(const float *a,int nrow,int ncol,
                          const int *row_assign,const int *col_assign,
                          int k,
                          float *accu) {
  int i,j;

  for(i=0;i<nrow;i++) for(j=0;j<ncol;j++) {
    accu[row_assign[i]*k+col_assign[j]]+=a[i*ncol+j];
  }

}

int *imat_joint_histogram(int n,int k,int *row_assign,int *col_assign) {
  int *hist=ivec_new_0(k*k);
  int i;

  for(i=0;i<n;i++) 
    hist[row_assign[i]*k+col_assign[i]]++;

  return hist;
}
