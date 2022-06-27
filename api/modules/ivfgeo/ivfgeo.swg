// -*- c++ -*- 

%module ivfgeo;




%{

#include "ivfgeo.h"
#include "ivfgeo_pq.h"

%}

%include <yael/common.swg>



PARALLEL_FUNCTION_PEEK(ivfgeo_query_peek)
PARALLEL_FUNCTION(ivfgeo_query)
PARALLEL_FUNCTION(ivfgeo_add)

// so that invfile_query_* functions return
// (labels,distances) 
// labels=C table of image labels
// distances=C table of distances
%typemap(in,numinputs=0)
        (distype_t ** retdis) 
        (distype_t *disp) {
  $1=&disp;
}

%typemap(argout) (distype_t ** retdis)  {
  PyObject *prev_res=$result;
  $result=Py_BuildValue
        ("(NN)",
        prev_res,
        SWIG_NewPointerObj(*$1,SWIGTYPE_p_float,0));
}



typedef struct pointset_s pointset_t;

typedef pointset_t vwgeoset_t;

%include "ivfgeo.h"



typedef struct ivfpq_t ivfpq_t;
typedef struct pq_s pq_t;

PARALLEL_FUNCTION(ivfgeo_pq_add) 
PARALLEL_FUNCTION(ivfgeo_pq_query)

PARALLEL_FUNCTION(ivfgeo_pq_add_with_vw)

%include "ivfgeo_pq.h"
