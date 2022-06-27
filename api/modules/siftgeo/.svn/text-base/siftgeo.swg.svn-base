// -*- c++ -*-

%module siftgeo

/*************************************************************************
 * wrapping of siftgeo.h, nn.h and libit
 * 
 * All the C functions are wrapped (even if they are not
 * used). Consequence: when a function is defined in a .h but not
 * implemented in a .c, it generates linking errors when the Python
 * module is imported.
 *
 * Most of the SWIG code is typemaps for functions that return results
 * by pointer arguments: they are converted to Python tuples.
 * 
 * The only object for which memory is managed automatically is
 * pointset_t. The other ones should be allocated and deallocated
 * manually in Python.
 *
 */


%{

#include "siftgeo.h"
#include "siftgeo_binarize.h"
#include "siftgeo_and_nn.h"
#include "siftgeo_agg.h"

typedef pointset_t* pointset_t_p;

int compute_point_size(point_t * corners, int nc, int fmt)
{
  int size;

  if(nc==0) return 0;

  switch (fmt) {
  case 0:
    /* This assumes that all descriptors are filled and all have the same size */
    size = sizeof(corners->dim) + sizeof(corners->geom) + corners->dim;
    break;
  case 1:
    size = sizeof(corners->vw) + sizeof(corners->geom);
    break;
  case 2:
    size = sizeof(corners->vw) + sizeof(corners->geom) + sizeof (corners->binsign);
    break;
  default:
    fprintf (stderr, "compute_point_size: invalid format %d\n", fmt);
    size = -1;
  }

  return size;
}

pointset_t *pointset_showup(pointset_t *a) {
  return a;
}

vwgeoset_t *vwgeoset_cast(void *ptr) {
  return (vwgeoset_t *) ptr;
}

%}


%include <yael/common.swg>


/*************************************************************************
 * siftgeo.h
 */


// read_points_add_with_mask takes a point_t * as input, returns 
// (point_t *,int) as output

%typemap(in) (point_t ** corners_io) (point_t *corners) {
  void *vcorners;
  int trycv=SWIG_ConvertPtr($input,&vcorners,SWIGTYPE_p_point_t, 0);
  if(!SWIG_IsOK(trycv)) {
    PyErr_SetString(PyExc_TypeError,"need point_t array");
    return NULL;    
  }
  corners=(point_t*)vcorners;
  $1=&corners;
}

%typemap(argout) (point_t ** corners_io)  {
  $result=Py_BuildValue
        ("(NN)",
         SWIG_NewPointerObj(*$1,SWIGTYPE_p_point_t,0),
         $result);
}


// read_points_string takes a python string
%typemap(in,numinputs=1)
        (const char *string,int string_length) {
  if(!PyString_Check($input)) {
    PyErr_SetString(PyExc_TypeError,"need string");
    return NULL;
  }  
  $2=PyString_Size($input);
  $1=PyString_AsString($input);
}


// read_points_* returns a tuple
// (points,npoint)

%typemap(in,numinputs=0) 
  (point_t ** corners_out)
  (point_t *corners) {
  $1=&corners;
}

%typemap(argout) (point_t ** corners_out)  {
  $result=Py_BuildValue
        ("(NN)",         
         SWIG_NewPointerObj(*$1,SWIGTYPE_p_point_t,0),
         $result);
}


/* hack to access both dim and vw fields. Generates a warning */ 

typedef struct {
  geom_t geom;

  int dim;                    /* dimension of the descriptor (for siftgeo) */
  int vw;                     /* visual word (for vwgeo) */

  unsigned char *desc;
  
  binsign_t binsign;
} point_t;


// pointset

%newobject pointset_new;
%newobject pointset_read;
%newobject pointset_read_file_max;
%newobject pointset_read_cornerness;
%newobject pointset_dup;
%newobject pointset_from_vw;
%newobject pointset_alloc;


/* parallelize some i/o functions (may be useful for slow NFS links) */

PARALLEL_FUNCTION(pointset_read);
PARALLEL_FUNCTION(pointset_read_cornerness);
PARALLEL_FUNCTION(siftgeo_read_dispatch);
PARALLEL_FUNCTION(read_points);
PARALLEL_FUNCTION(count_points);
PARALLEL_FUNCTION(read_points_file);
PARALLEL_FUNCTION(pointset_read_file_max);
PARALLEL_FUNCTION(pointset_write);
PARALLEL_FUNCTION(write_points_file);
PARALLEL_FUNCTION(write_points);
PARALLEL_FUNCTION(pointset_crop);
PARALLEL_FUNCTION(vwgeoset_sort);
PARALLEL_FUNCTION(pointset_dup);
PARALLEL_FUNCTION(pointset_append);
PARALLEL_FUNCTION(pointset_affine_transform);
PARALLEL_FUNCTION(pointset_crop_polygon);
PARALLEL_FUNCTION(pointset_crop);
PARALLEL_FUNCTION(pointset_filter_cornerness)
PARALLEL_FUNCTION(vwgeoset_filter_n_cornerness_max)


%include "siftgeo.h"

vwgeoset_t *vwgeoset_cast(void *ptr);

%extend pointset_t {
  pointset_t() {return pointset_new(); }

  /*  hack: we want to be able to manage pointers comming from other
     modules, so shortcircuit the type system... */
  pointset_t(void *ps) {return ps; }

  ~pointset_t() {pointset_delete($self); }
  point_t *__getitem__(int i) {
    assert(i>=0 && i<$self->n);
    return &$self->pts[i];
  }
  PyObject *tostring(int fmt) {
    int       size = $self->n * compute_point_size($self->pts, $self->n, fmt);
    PyObject *str = PyString_FromStringAndSize(NULL, size);
    char     *w   = (char *) PyString_AsString(str);

    if (write_points_string(w, size, $self->pts, $self->n, fmt) < 0) {
      Py_INCREF(Py_None);
      return Py_None;
    }
    return str;
  }
  void fromstring(PyObject *str, int fmt) {
    if (!PyString_Check(str)) {
      return;
    }

    char    *string = (char *) PyString_AsString(str);
    int      size   = PyString_Size(str);

    $self->n = read_points_string(string, size, &$self->pts, fmt);
  }
  int hasdescs() {
    /* returns: 0 for empty set
     *          1 for set with only vwgeos
     *          2 for set with only siftgeos
     *          3 for mixed vwgeo/siftgeo set
     */
    int hd=0,i;
    for(i=0;i<$self->n;i++) 
      hd|=$self->pts[i].desc ? 2 : 1;
    return hd;
  }
};

%include <carrays.i>

// see the PointsetArray class @ end
%array_functions(pointset_t*,PointsetArray);



// siftgeo_to_ffq returns a tuple
// (floats,d,n)
// floats= float table
// d=dimension of the pts (128)
// n=nb of pts or -1 on error
%typemap(in,numinputs=0) 
  (int * d_out)
  (int d) {
  $1=&d;
}

%typemap(argout) (int* d_out)  {
  $result=Py_BuildValue
        ("(Ni)",
         $result,*$1);
}

%include "siftgeo_and_nn.h"


/*************************************************************************
 * siftgeo_binarize.h
 */


PARALLEL_FUNCTION(siftgeo_binarize_fill_medians);
PARALLEL_FUNCTION(siftgeo_binarize_ffq);
PARALLEL_FUNCTION(siftgeo_binarize);

%include "siftgeo_binarize.h"


/*************************************************************************
 * siftgeo_agg.h
 */


// wrap clustering_kmeans_assign to return (clust_assign,*cluster_size)
%typemap(in,numinputs=0)
  (int ** idx_out,  float **v_out) 
  (int * idx, float * v) {
  $1 = &idx;
  $2 = &v;
}

%typemap(argout) (int ** idx_out,  float **v_out)  {
  $result=Py_BuildValue
        ("(NNN)",
         $result,
	 SWIG_NewPointerObj(*$1,SWIGTYPE_p_int,0),
	 SWIG_NewPointerObj(*$2,SWIGTYPE_p_float,0));
}



PARALLEL_FUNCTION(ldc_saliency);
PARALLEL_FUNCTION(vwsgeo_to_bof);
PARALLEL_FUNCTION(atoms_signature);
PARALLEL_FUNCTION(atoms_signature_fvec);
PARALLEL_FUNCTION(hesp_signature);

%include "siftgeo_agg.h"



/*************************************************************************
 * Locally defined functions
 */



%pythoncode %{

class PointsetArray:
  """ array owns its pointsets """

  def __init__(self,n):
    self.this=new_PointsetArray(n)
    self.n=n
    for i in range(n):
      x=pointset_t(); x.thisown=False
      PointsetArray_setitem(self.this,i,x)
      
  def __getitem__(self,i):
    return PointsetArray_getitem(self.this,i)

  def __setitem__(self,i,ps):
    x=PointsetArray_getitem(self.this,i)
    x.thisown=True # dealloc old one
    PointsetArray_setitem(self.this,i,ps)        
    ps.thisown=False

  def __del__(self):
    for i in range(self.n):
      x=PointsetArray_getitem(self.this,i)
      x.thisown=True
    delete_PointsetArray(self.this)


def unpickle_pointset_t(fmt,s,d):
  ps=pointset_t()
  ps.fromstring(s,fmt)
  return ps

def pointset_t_reduce(self):
  hd=self.hasdescs()
  if hd in (0,1): fmt=2
  elif hd==2: fmt=0
  else: raise RuntimeError("cannot pickle mixed pointset")  
  return (unpickle_pointset_t,(fmt,self.tostring(fmt),{}))

pointset_t.__reduce__=pointset_t_reduce


%}
