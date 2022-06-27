#ifndef __matrix_h
#define __matrix_h

/*---------------------------------------------------------------------------*/
/*! @addtogroup matrix
 *  @{
 */
/*---------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*/
/* Standard operations                                                       */
/*---------------------------------------------------------------------------*/

/*! @brief Allocate a new nrow x ncol matrix */
float *fmat_new (int nrow, int ncol);


/*! @brief Matrix multiplication
 *
 * computes mout = left * right
 * where 
 *   mout     is n-by-k 
 *   left     is n-by-m
 *   right    is m-by-k
 * (all matrices stored by lines, like in C, and packed)
 */
void fmat_mul (const float *left, const float *right,
	       int n, int m, int k, float *mout);

float * fmat_new_mul (const float *left, const float *right,
		      int n, int m, int k);

/*! @brief Same as fmat_mul, but transpose left matrix (left of size m x n) */
void fmat_mul_tl (const float *left, const float *right,
		  int n, int m, int k, float *mout);

float *fmat_new_mul_tl (const float *left, const float *right, 
			int n, int m, int k);

/*! @brief Same as fmat_mul, but transpose right matrix (right of size k x m) */
void fmat_mul_tr (const float *left, const float *right,
		  int n, int m, int k, float *mout);

float *fmat_new_mul_tr (const float *left, const float *right, 
			int n, int m, int k);

/*! @brief Same as fmat_mul, but transpose both left and right matrices
  left is of size m * n and right of size k x m */
void fmat_mul_tlr (const float *left, const float *right,
		   int n, int m, int k, float *mout);

float *fmat_new_mul_tlr (const float *left, const float *right, 
			int n, int m, int k);

/*! @brief Multiply a matrix by a vector */
float * fmat_mul_fvec (const float * a, const float * v, int nrow, int ncol);

/*! @brief display the matrix in matlab-like format */
void fmat_display (const float *a, int nrow, int ncol);



/*---------------------------------------------------------------------------*/
/* Matrix manipulation functions                                             */
/*---------------------------------------------------------------------------*/

/*! @brief return the submatrix defined by left-upper corner (included) 
  and top-down corner (not included) */
float *fmat_get_submatrix (const float *a, int ncola, int r1, int c1, int r2, int c2);

/*! @brief produce a matrix composed of the rows indicated by the vector rows */
float *fmat_get_rows (const float *a, int ncol, int nrowout, const int *rows);


/*! 
 * a is ncol-by-nrow
 * accu is k-by-k
 *
 * for i=0..ncol-1,j=0..nrow-1, do 
 *    accu(row_assign[i],col_assign[j]) += a(i,j)
 *
 */ 
void fmat_splat_separable(const float *a,int nrow,int ncol,
                          const int *row_assign,const int *col_assign,
                          int k,
                          float *accu); 

int *imat_joint_histogram(int n,int k,int *row_assign,int *col_assign);

/*---------------------------------------------------------------------------*/
/* Special matrices                                                          */
/*---------------------------------------------------------------------------*/


/*! @brief produce a new matrix of size nrow x ncol, filled with gaussian values */
float * fmat_new_rand_gauss (int nrow, int ncol);

/*! @brief produce a random orthogonal basis matrix of size d*d */
float *random_orthogonal_basis (int d);

/*! @brief Construct a Hadamard matrix of dimension d using the Sylvester construction.
   d should be a power of 2 */
float * hadamard (int d);


/*---------------------------------------------------------------------------*/
/* Statistical matrix operations                                             */
/*---------------------------------------------------------------------------*/

/*! @brief Perform the Principal Component Analysis of a set of vectors,
 * v(n,d) stored per row
 *  return d*d matrix of eigenvectors */
float *compute_pca(int n,int d,float *v);


/*---------------------------------------------------------------------------*/
/*! @} */
/*---------------------------------------------------------------------------*/

#endif 

/*---------------------------------------------------------------------------*/

