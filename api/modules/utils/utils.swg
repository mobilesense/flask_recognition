// -*- c++ -*- 

%module utils;


/*****************************************************
 * a few useful functions & objects for everywhere:
 * - array objects 
 * - nearest neighbour and clustering functions
 */


%{

#include "hkm.h"
#include "pq.h"
#include "pqr.h"
#include "ivfpq.h"
#include "trees.h"
#include "annvw.h"

%}


%include <yael/common.swg>


// slow functions with peek


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


%include "hkm.h"



/*************************************************************************
 * trees.h
 */

%typemap(in,numinputs=0) 
  (long *k_out,float *v_out)
  (long k,float v) {
  $1=&k;
  $2=&v;
}

%typemap(argout) (long *k_out,float *v_out)  {
  Py_DECREF($result);
  $result=Py_BuildValue("(lf)",*$1,(double)(*$2));
}



%include "trees.h"


/*************************************************************************
 * pq.h
 */

PARALLEL_FUNCTION(pq_encode)
PARALLEL_FUNCTION(pq_decode)
PARALLEL_FUNCTION(pq_asym_nns)
PARALLEL_FUNCTION(pq_asym_fc_nns)
PARALLEL_FUNCTION(pq_nns)

%include "pq.h"

/*************************************************************************
 * pqr.h
 */

PARALLEL_FUNCTION(pqr_encode)
PARALLEL_FUNCTION(pqr_decode)
PARALLEL_FUNCTION(pqr_nns)

%include "pqr.h"


/*************************************************************************
 * ivfpq.h
 */


typedef struct pq_s  pq_t;

PARALLEL_FUNCTION(ivfpq_query)
PARALLEL_FUNCTION(ivfpq_add)


%include "ivfpq.h"

/*************************************************************************
 * annvw.h
 */


typedef struct annvw_s  annvw_t;

PARALLEL_FUNCTION_PEEK(quantize_codebook_annvw_full_thread)
PARALLEL_FUNCTION(quantize_codebook_annvw_multiple_thread)
PARALLEL_FUNCTION(quantize_codebook_annvw_multiple)
PARALLEL_FUNCTION(quantize_codebook_annvw_thread)

PARALLEL_FUNCTION_PEEK(assign_vw_ann_peek)
PARALLEL_FUNCTION_PEEK(assign_vw_peek)
PARALLEL_FUNCTION_PEEK(assign_vw_raw)

%include "annvw.h"
