// -*- c++ -*-

%module invfile


/* abuse of inlines in RV's code */
#define inline 



/*************************************************************************
 * Extra functions 
 */

%{

#include <siftgeo.h>
#include <invfile.h>
#include <clustering.h>

#include <it/io.h>
#include <it/mat.h>


PyObject *ivec_to_python(ivec v) {
  int i,n=ivec_length(v);
  PyObject *ret=PyList_New(n);
  for(i=0;i<n;i++) 
    PyList_SetItem(ret,i,PyInt_FromLong(v[i]));
  return ret;
}

PyObject *vec_to_python(vec v) {
  int i,n=vec_length(v);
  PyObject *ret=PyList_New(n);
  for(i=0;i<n;i++) 
    PyList_SetItem(ret,i,PyFloat_FromDouble(v[i]));
  return ret;
}

ivec python_to_ivec(PyObject *l) {
  ivec ret;
  int i;
  if(!PyList_Check(l))  {
    PyErr_SetString(PyExc_TypeError,"need list of integers");
    return NULL;
  }
  ret=ivec_new(PyList_Size(l));
  for(i=0;i<PyList_Size(l);i++) {
    PyObject *o=PyList_GetItem(l,i);
    if(! PyInt_Check(o)) {
      PyErr_SetString(PyExc_TypeError,"need list of integers");
      ivec_delete(ret);
      return NULL;
    }
    ret[i]=PyInt_AsLong(o);
  }
  return ret;
}

typedef struct {
  PyObject *callable;
  PyThreadState *_save;
} callback_t;

void call_peek_fun_thread(void *arg,double frac) {
  callback_t *callback=(callback_t*)arg;
  PyEval_RestoreThread(callback->_save);  
  {
    PyObject *obj=PyObject_CallFunction(callback->callable,"d",frac);
    if(!obj) {
      fprintf(stderr,"!!!exception in peek callback\n");
      PyErr_Print();
    } else
      Py_DECREF(obj);   
  }
  callback->_save=PyEval_SaveThread();
}

ihash_t *new_basic_ihash(int ncell) {
  return ihash_new (ncell, hash_basic_hash_function);
}

#define NS 23

void stress_mem(void) {
  int i;
  void *d[NS];
  for(i=0;i<NS;i++) 
    d[i]=calloc(1,1<<i);
  for(i=0;i<NS;i++) 
    free(d[i]);  
}

ivec imat_getitem(imat m,int i) {
  assert(i>=0 && i<imat_height(m));
  return m[i];
}

int ivec_getitem(ivec v,int i) {
  assert(i>=0 && i<ivec_length(v));
  return v[i];
}
void ivec_setitem(ivec v,int i,int a) {
  assert(i>=0 && i<ivec_length(v));
  v[i]=a;
}

void ivec_to_sparse (ivec v, ivec * svi_out, ivec * sv_out)
{
  int  n = ivec_length (v);
  int  nz = 0, i, j;
  ivec svi, sv;
  for (i = 0; i < n; i++)
    if (v[i] != 0)
      nz++;
  svi = ivec_new (nz);
  sv = ivec_new (nz);
  j = 0;
  for (i = 0; i < n; i++)
    if (v[i] != 0) {
      svi[j] = i;
      sv[j] = v[i];
      j++;
    }
  *svi_out = svi;
  *sv_out = sv;
}

ivec ivec_from_sparse (ivec svi, ivec sv)
{
  int n=ivec_max(svi);
  int nz=ivec_length(svi);
  ivec ret=ivec_new_zeros(n);
  int  i;
  
  for(i=0;i<nz;i++) 
    ret[svi[i]]=sv[i];
  return ret;
}


/* cancels renaming of callback_t *cb */
#define cb1 cb
#define cb2 cb
#define cb3 cb
#define cb4 cb
#define cb5 cb
#define cb6 cb
#define cb7 cb
#define cb8 cb
#define cb9 cb



%}

/*************************************************************************
 * Wrapper code
 */

%typemap(in) FILE* {
  if(!PyFile_Check($input)) {
    PyErr_SetString(PyExc_TypeError,"need regular file object");
    return NULL;
  }
  $1=((PyFileObject*)$input)->f_fp;
}

%include <it/vec.h>
%include <it/mat.h>


/* parallel functions with a peek callback need this */

%typemap(in, numinputs=1) 
    (void (*peek_fun)(void *arg,double frac),void *peek_arg) 
    (callback_t cb)
{
  if(Py_None==$input) {
    $1=NULL;
    $2=NULL; 
  } else {
    cb.callable=$input;
    $1=call_peek_fun_thread;
    $2=&cb;
  }
}


// function does not need python global interpreter lock

%define PARALLEL_FUNCTION(funcname)
%exception funcname {
  Py_BEGIN_ALLOW_THREADS
  $action
  Py_END_ALLOW_THREADS
}
%enddef

// same as before, but must be locked from time to time to 
// show progress
%define PARALLEL_FUNCTION_PEEK(funcname)
%exception funcname {
  cb._save = PyEval_SaveThread();
  $action
  PyEval_RestoreThread(cb._save);
}
%enddef


PARALLEL_FUNCTION_PEEK(invfile_query_L1_norm_peek)
PARALLEL_FUNCTION_PEEK(invfile_query_custom_norm_peek)

// for invfile_query_L1_norm 

%typemap(in,numinputs=0)
        (vec * dis) 
        (vec disp) {
  $1=&disp;
}

%typemap(argout) (vec *dis)  {
  PyObject *prev_res=$result;
  $result=Py_BuildValue
        ("(NN)",
        prev_res,
        SWIG_NewPointerObj(*$1,SWIGTYPE_p_double,0));
}



%ignore invfile_query_load_cdm_factors;
%ignore invfile2l_compute_medoid_norm;
%ignore cvec_new_reverse;
%ignore cmat_ones;

PARALLEL_FUNCTION(invfile_read)
PARALLEL_FUNCTION(invfile_write)
PARALLEL_FUNCTION(invfile_compute_vec_norm_L1)
PARALLEL_FUNCTION(invfile_compute_vec_norm_L1_tfidf)



%include <invfile.h>


PyObject *ivec_to_python(ivec v);
PyObject *vec_to_python(vec v);
ivec python_to_ivec(PyObject *l);

// wrap a few i/o functions 
//%include <it/io.h>
/* Write and read matrix in a pack format (dimension+data)              */

  void vec_fwrite (FILE * stream, vec v);
  void fvec_fwrite (FILE * stream, vec v);
  void bvec_fwrite (FILE * stream, bvec v);
  void ivec_fwrite (FILE * stream, ivec v);
  void cvec_fwrite (FILE * stream, cvec v);


/* Return 1 if the vector is successfully read, 0 otherwise */
  int  vec_fread (FILE * stream, vec v);
  int  fvec_fread (FILE * stream, vec v);
  int  bvec_fread (FILE * stream, bvec v);
  int  ivec_fread (FILE * stream, ivec v);
  int  cvec_fread (FILE * stream, cvec v);

  /* return NULL if the vector can't be read */
  vec  vec_new_fread (FILE * stream);
  vec  fvec_new_fread (FILE * stream);
  bvec bvec_new_fread (FILE * stream);
  ivec ivec_new_fread (FILE * stream);
  cvec cvec_new_fread (FILE * stream);

  void mat_fwrite (FILE * stream, mat m);
  void bmat_fwrite (FILE * stream, bmat m);
  void imat_fwrite (FILE * stream, imat m);
  void cmat_fwrite (FILE * stream, cmat m);

  mat  mat_new_fread (FILE * stream);
  bmat bmat_new_fread (FILE * stream);
  imat imat_new_fread (FILE * stream);
  cmat cmat_new_fread (FILE * stream);


%include <cpointer.i>

%pointer_functions(ivec,ivecptr)

%typemap(in,numinputs=0)
        (ivec *svi_out,ivec *sv_out) 
        (ivec svi,ivec sv) {
  $1=&sv;
  $2=&svi;
}

%typemap(argout) (ivec *svi_out,ivec *sv_out)  {
  Py_DECREF($result);
  $result=Py_BuildValue
        ("(NN)",
        SWIG_NewPointerObj(*$1,SWIGTYPE_p_int,0),
        SWIG_NewPointerObj(*$2,SWIGTYPE_p_int,0));
}
%ignore ann_weights_inner_product;

PARALLEL_FUNCTION_PEEK(assign_vw_ann_peek)
PARALLEL_FUNCTION_PEEK(assign_vw_peek)


%typemap(in,numinputs=1)
        (const char *string,int string_length) {
  if(!PyString_Check($input)) {
    PyErr_SetString(PyExc_TypeError,"need string");
    return NULL;
  }  
  $2=PyString_Size($input);
  $1=PyString_AsString($input);
}


%include <siftgeo.h>

%ignore clustering_random;
%ignore clustering_random_pos;

%include <clustering.h>
%include <ann.h>


%include <hash.h>


ihash_t *new_basic_ihash(int ncell);

void stress_mem(void);

ivec imat_getitem(imat m,int i);
int ivec_getitem(ivec,int i);
void ivec_setitem(ivec,int i,int v);


void ivec_to_sparse (ivec v, ivec * svi_out, ivec * sv_out);
ivec ivec_from_sparse (ivec svi, ivec sv);
