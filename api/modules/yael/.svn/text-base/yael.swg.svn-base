// -*- c++ -*- 

%module yael;


/*****************************************************
 * a few useful functions & objects for everywhere:
 * - array objects 
 * - nearest neighbour and clustering functions
 */


%{

#include "nn.h"
#include "clustering.h"
#include "sorting.h"
#include "vector.h"
#include "matrix.h"
#include "machinedeps.h"

/* seed generators used in the rest */
void common_srandom(int seed) {
  srandom(seed);
  srand(seed);
  srand48(seed); 
}

typedef int * int_pointer_t;
typedef float * float_pointer_t;

%}

%include "common.swg"


%include <carrays.i>


/* 
   a=IntArray(4) constructs an array of 4 ints  
   a.cast() returns an int* usable as a function argument
   if x is a int* in C,  
   b=IntArray.frompointer(x) makes a[i] valid to access x[i]
   b=IntArray.acquirepointer(x) will also free(x) when b is deleted
   b.plus(2) returns a pointer to x+2
*/
%define array_class_plus(type,name)
%array_class(type,name)
%extend name {
  type *plus(long i) {
    return $self+i;
  }
  void clear(long n) {
    memset($self,0,sizeof(type)*n);
  }
  void copyfrom(type *src,long dest_ofs,long n) {
    memcpy($self+dest_ofs,src,sizeof(type)*n);
  } 
  PyObject *tostring(long n) {
    PyObject *str = PyString_FromStringAndSize(NULL, n*sizeof(type));
    memcpy(PyString_AsString(str),$self,sizeof(type)*n);
    return str;
  }
  void fromstring(PyObject*obj) {
    if(!PyString_Check(obj)) {
      fprintf(stderr,"need string!");
      return;
    }
    memcpy($self,PyString_AsString(obj),PyString_Size(obj));
  }
};
%pythoncode %{
name.acquirepointer=staticmethod(lambda ptr: _frompointer_and_acquire(name,ptr))
%}
%enddef

%pythoncode %{
def _frompointer_and_acquire(aclass,ptr):
  r=aclass.frompointer(ptr)
  if r: r.this.acquire()
  return r
%}


array_class_plus(double,DoubleArray)
array_class_plus(float,FloatArray)
array_class_plus(int,IntArray)

typedef int * int_pointer_t;
typedef float * float_pointer_t;

array_class_plus(int_pointer_t,IntPtrArray)
array_class_plus(float_pointer_t,FloatPtrArray)





/*************************************************************************
 * nn.h
 */



// fvecs_new_read and fvecs_new_fread_max returns a tuple 
// (floats,n,d)
// floats= float table
// d=dimension of the pts (128)
// n=nb of centroids or -1 on error

%typemap(in,numinputs=0) 
  (int *d_out,float **vf)
  (int d,float *vf) {
  $1=&d;
  $2=&vf;
}

%typemap(argout) (int *d_out,float **vf)  {
  $result=Py_BuildValue
        ("(NNi)",
         SWIG_NewPointerObj(*$2,SWIGTYPE_p_float,0),
         $result,*$1);
}


// ivecs_new_read returns a tuple 
// (ints,n,d)
// ints= float table
// d=dimension of the pts (128)
// n=nb of centroids or -1 on error

%typemap(in,numinputs=0) 
  (int *d_out,int **vi)
  (int d,int *vi) {
  $1=&d;
  $2=&vi;
}

%typemap(argout) (int *d_out,int **vi)  {
  $result=Py_BuildValue
        ("(NNi)",
         SWIG_NewPointerObj(*$2,SWIGTYPE_p_int,0),
         $result,*$1);
}

// for ivec_new_read

%typemap(in,numinputs=0) 
  (int *d_out)
  (int d) {
  $1=&d;
}

%typemap(argout) (int *d_out)  {
  $result=Py_BuildValue
        ("(Ni)",
         $result,*$1);
}




// slow functions with peek
PARALLEL_FUNCTION_PEEK(quantize_codebook)
PARALLEL_FUNCTION_PEEK(quantize_codebook_multiple)
PARALLEL_FUNCTION_PEEK(quantize_codebook_multiple_thread)
PARALLEL_FUNCTION_PEEK(quantize_codebook_thread)

PARALLEL_FUNCTION_PEEK(quantize_codebook_full)
PARALLEL_FUNCTION_PEEK(quantize_codebook_full_thread)

PARALLEL_FUNCTION(compute_cross_distances)
PARALLEL_FUNCTION(compute_cross_distances_alt)


%include "nn.h"


/*************************************************************************
 * clustering.h
 */


// wrap clustering_kmeans_assign to return (clust_assign,*cluster_size)
%typemap(in,numinputs=0)
        (int ** clust_assign_out) 
  (int *clust_assign) {
  $1=&clust_assign;
}

%typemap(argout) (int ** clust_assign_out)  {
  $result=Py_BuildValue
        ("(NN)",
         $result,
        SWIG_NewPointerObj(*$1,SWIGTYPE_p_int,0));
}




%include "clustering.h"

/*************************************************************************
 * sorting.h
 */

PARALLEL_FUNCTION(find_labels)
PARALLEL_FUNCTION(fvec_find_k_max)
PARALLEL_FUNCTION(fvec_find_k_min)
PARALLEL_FUNCTION(fvec_ranks_of)
PARALLEL_FUNCTION(fvec_ranks_inc_of)
PARALLEL_FUNCTION(fvec_median)
PARALLEL_FUNCTION(fvec_quantile)
PARALLEL_FUNCTION(ivec_sort)
PARALLEL_FUNCTION(ivec_sort_index)
PARALLEL_FUNCTION(ivec_sort_by_permutation)
PARALLEL_FUNCTION(merge_ordered_sets)


%typemap(in,numinputs=0)
        (int **labels_out,float **vals_out) 
  (int *labels,float *vals) {
  $1=&labels;
  $2=&vals;
}

%typemap(argout) (int **labels_out,float **vals_out) {
  $result=Py_BuildValue
        ("(NNN)",
         $result,
         SWIG_NewPointerObj(*$1,SWIGTYPE_p_int,0),
         SWIG_NewPointerObj(*$2,SWIGTYPE_p_float,0));
}




%include "sorting.h"

/*************************************************************************
 * vector.h
 */

PARALLEL_FUNCTION(fvecs_read)
PARALLEL_FUNCTION(fvecs_write)
PARALLEL_FUNCTION(fvec_fwrite)
PARALLEL_FUNCTION(ivec_fwrite)
PARALLEL_FUNCTION(ivec_count)
PARALLEL_FUNCTION(fvecs_new_read)
PARALLEL_FUNCTION(fvecs_new_fread_max)
PARALLEL_FUNCTION(fvec_read)
PARALLEL_FUNCTION(fvecs_new_read_sparse)
PARALLEL_FUNCTION(fvec_splat_add)
PARALLEL_FUNCTION(fvec_isplat_add)
PARALLEL_FUNCTION(fvec_distance_L2sqr)

%typemap(in,numinputs=0) 
  (int ** idx_out, float ** v_out)
  (int *idx,float *v) {
  $1=&idx;
  $2=&v;
}

%typemap(argout) (int ** idx_out, float ** v_out)  {
  $result=Py_BuildValue
        ("(NNN)",
         $result,
         SWIG_NewPointerObj(*$1,SWIGTYPE_p_int,0),
         SWIG_NewPointerObj(*$2,SWIGTYPE_p_float,0));
}


%typemap(in,numinputs=0) 
  (float ** vf_out)
  (float *vf) {
  $1=&vf;
}

%typemap(argout) (float ** vf_out)  {
  $result=Py_BuildValue
        ("(NN)",
         $result,
         SWIG_NewPointerObj(*$1,SWIGTYPE_p_float,0));
}




%include "vector.h"



/*************************************************************************
 * matrix.h
 */


%include "matrix.h"


/*************************************************************************
 * allows freeing some elements allocated by malloc/calloc in python
 */
void free (void *);


%include "machinedeps.h"

void common_srandom(int seed);

