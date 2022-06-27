#ifndef __vector_h
#define __vector_h

#include <stdio.h>

/*---------------------------------------------------------------------------
 * Vectors are represented as C arrays of basic elements. Functions
 * operating on them are prefixed with:
 *
 * ivec_: basic type is int
 * fvec_: basic type is float
 *
 * Vector sizes are passed explicitly, as long int's to allow for
 * large arrays on 64 bit machines. Vectors can be free'd with free().
 *
 * Matrix lines are stored contiguously in memory. Element (i,j) of a
 * matrix vf of n lines and d columns can be accessed
 * with vf[ i * d + j ]
 *
 *---------------------------------------------------------------------------*/

/* some operations may be faster if input arrays are allocated/free'd
 * with this pair of functions (data is suitably aligned).
 */

/*! @brief Alloc a new aligned vector of floating point values. To be de-allocated with free. */
float * fvec_new (long n);

/*! @brief Alloc an int array. To be de-allocated with free. */
int *ivec_new(long n);

/*! @brief initialized with 0's */
float *fvec_new_0 (long n);
int *ivec_new_0(long n);

/*! @brief initialized with some value */
float *fvec_new_set(long n, float val);
int *ivec_new_set(long n, int val);

/*! initialize with uniformly drawn nbs in [0,1) */
float *fvec_new_rand(long n);

/* @brief new vector [a,a+1,...b-1] */
int * ivec_new_range (long a, long b);

/*! @brief new vector initialized with another vector */
int * ivec_new_copy (const int * v, long n);

/*! @brief new vector initialized with another vector */
float * fvec_new_copy (const float * v, long n);

/*! @brief random permutation of 0..n-1 */ 
int *ivec_new_random_perm(long n);

/*! count occurrences

   k is the range of the values that may be encountered (assuming start at 0)
   v is the vector of values to be histrogramized, of length n
*/
int * ivec_new_histogram (int k, int * v, int n);

/*! compute a hash value for the vector */
int ivec_hash(const int * v, long n);

/*! count occurences of a value in the vector */
int ivec_count_occurrences(const int * v, long n, int val);


/*---------------------------------------------------------------------------*/
/* Input/Output functions                                                    */
/* I/O of a single vector is supported only if it is smaller than 2^31       */ 
/*---------------------------------------------------------------------------*/

/*! @brief write a vector into an open file */
int ivec_fwrite(FILE *f, const int *v, int d);
int fvec_fwrite(FILE *f, const float *v, int d);

/* write without header */
int fvec_fwrite_raw(FILE *f, const float *v, long d);

/*! @brief write a set of vectors into an open file */
int ivecs_fwrite(FILE *f, int d, int n, const int *v);
int fvecs_fwrite (FILE *fo, int d, int n, const float *vf);

/*! @brief sevral integer vectors into an file */
int ivecs_write(const char *fname, int d, int n, const int *v);
int fvecs_write (const char *fname, int d, int n, const float *vf);



/*! @brief load float vectors from file.
 *
 * Returns nb of vectors read, or <0 on error
 */
int fvecs_new_read (const char *fname, int *d_out, float **vf);


int fvecs_new_fread_max (FILE *f, int *d_out, float **vf, long nmax);

/*! reads sparse vectors and return them as dense. d must be known */
int fvecs_new_read_sparse (const char *fname, int d, float **vf_out);

/*! @brief load float vector without allocating memory 
 *
 * Fills n*d array with as much vectors read from fname as possible.
 * Returns nb of clusters read, or <0 on error.
 */
int fvecs_read (const char *fname, int d, int n, float *a);


/*! @brief read a single vector from a file

 * Fill a with a single float vector from fname offset o_f into file a
 * Returns <0 on error
 */
int fvec_read (const char *fname, int d, float *a, int o_f);

/*! @brief load float vectors from an open file. Return the dimension */
int fvec_fread (FILE * f, float * v);


float *fvec_fread_raw(FILE * f, long n);


/*! @brief read and allocate a an integer vector file */
int * ivec_new_read(const char *fname, int *d_out);

/*! @brief read an integer vector file from an open file and return the dimension */
int ivec_fread (FILE *f, int * v);

/*! @brief read several integer vectors from an ivec file. Return number read */
int ivecs_new_read (const char *fname, int *d_out, int **vi);


/*! @brief display a float vector */
void fvec_print (const float * v, int n);
void fvec_fprintf (FILE * f, const float *v, int n, const char *fmt);

/*! @brief display an integer vector */
void ivec_print (const int * v, int n);
void ivec_fprintf (FILE * f, const int *v, int n, const char *fmt);

/* find first index of val (return -1 if not found) */
long ivec_index(const int * v, long n,int val);

/*---------------------------------------------------------------------------*/
/* Vector manipulation and elementary operations                             */
/*---------------------------------------------------------------------------*/

/*! @brief Set all the components of the vector v to 0 */
void fvec_0 (float * v, int n);
void ivec_0 (int * v, int n);

/*! @brief Set all the components of the vector v to the value val */
void fvec_set (float * v, int n, float val);
void ivec_set (int * v, int n, int val);

/*! @brief Increment or decrement a vector by a scalar value */
void fvec_incr (float * v, int n, double scal);
void fvec_decr (float * v, int n, double scal);

/*! @brief Multiply or divide a vector by a scalar */
void fvec_mul_by (float * v, int n, double scal);
void fvec_div_by (float * v, int n, double scal);

/*! @brief Add or subtract two vectors. The result is stored in v1. */
void fvec_add (float * v1, const float * v2, int n);
void fvec_sub (float * v1, const float * v2, int n);

/*! @brief Component-wise multiplication or division of two vectors (result in v1) */
void fvec_mul (float * v1, const float * v2, int n);
void fvec_div (float * v1, const float * v2, int n);

/*! @brief Normalize the vector for the given Minkowski norm. If the
  vector is all 0, it will be filled with NaNs. Bad luck! */
void fvec_normalize (float * v, int n, double norm);

/*! Standard functions: root square and square */
void fvec_sqrt (float * v, int n);
void fvec_sqr (float * v, int n);

/*! Replace the "Not a number values" by a given value */
int fvec_purge_nans(float * v, long n, float replace_value);

/*---------------------------------------------------------------------------*/
/* Vector measures and statistics                                            */
/*---------------------------------------------------------------------------*/

/*! @brief compute the sum of the value of the vector */
double fvec_sum (const float * v, int n);
long long ivec_sum (const int * v, int n);

/*! idem, squared */
long long ivec_sum_2 (const int * v, int n);


/*! @brief compute the norm of a given vector */
double fvec_norm (const float * v, int n, double norm);

/*! compute squared norm 2 */
double fvec_norm2sqr (const float * v, int n);

/*! @brief count the number of non-zeros elements */
int fvec_nz (const float * v, int n);
int ivec_nz (const int * v, int n);


/*! @brief compute the positions of the non-null positions.
  return the number of non-zeros positions. */
int fvec_find (const float *v, int n, int ** nzpos_out);
int ivec_find (const int *v, int n, int ** nzpos_out);

/*! @brief perform a random permutation on the elements of the vector */
void ivec_shuffle (int *v, long n);

/*! @brief entropy of the probability mass function represented by the vector */
double entropy (const float *pmf, int n);

/*! @brief entropy of a binary variable */
double binary_entropy (double p);

/*---------------------------------------------------------------------------*/
/* Distances and similarities                                                */
/*---------------------------------------------------------------------------*/

/*! @brief Return the Hamming distance (i.e., the number of different elements) */
int ivec_distance_hamming (const int * v1, const int * v2, int n);

/*! @brief Return the L2 distance between vectors */
double fvec_distance_L2 (const float * v1, const float * v2, int n);

/*! @brief Return the L1 distance between vectors */
double fvec_distance_L1 (const float * v1, const float * v2, int n);

/*! @brief Return the square L2 distance between vectors */
double fvec_distance_L2sqr (const float * v1, const float * v2, int n);

/*! @inner product between two vectors */
double fvec_inner_product (const float * v1, const float * v2, int n);

/*---------------------------------------------------------------------------
 * Sparse vector handling
 *
 * sparse vectors are represented with: int *idx, float *val, int nz.
 * 
 * Element idx[i] of vector is val[i], for i=0..nz-1
 *
 * for i=1..nz-1,  idx[i-1]<idx[i]
 *---------------------------------------------------------------------------*/


/*! @brief convert a vector to a sparse vector. 
  Return the number of non-zeros positions  */
int fvec_to_spfvec (float * v, int n, int ** idx_out, float ** v_out);
int ivec_to_spivec (int * v, int n, int ** idx_out, int ** v_out);


/*! @brief convert a sparse vector into a full vector */
float * spfvec_to_fvec (int * idx, float * v, int nz, int n);
int * spivec_to_ivec (int * idx, int * v, int nz, int n);

/*! @brief inner product between two sparse vectors */
float spfvec_inner_product (int *idx1, float *val1, int nz1, 
			    int *idx2, float *val2, int nz2);


/*---------------------------------------------------------------------------*/
/* Elaborate vector manipulations                                            */
/*---------------------------------------------------------------------------*/

/*! on output,

sl_out[0] =               v[      0] + ... + v[sl[0]-1]
sl_out[i] = sl_out[i-1] + v[sl[i-1]] + ... + v[sl[i]-1]

for 0<i<n

*/
void ivec_accumulate_slices(const int *v,int *sl,int n); 


/*! 
 * for i=0..n-1, do 
 *    accu[assign[i]] += a[i]
 */ 
void fvec_splat_add(const float *a,int n,
                    const int *assign,float *accu); 

/*! 
 * for i=0..n-1, do 
 *    accu[i] += a[assign[i]]
 */ 
void fvec_isplat_add(const float *a,int n,
                     const int *assign,float *accu); 



#endif
