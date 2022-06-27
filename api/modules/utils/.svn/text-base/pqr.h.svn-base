/* product quantizer with residual coding */

#ifndef __pqr_h
#define __pqr_h

#include "pq.h"


/* the residual coding structure */
typedef struct pqr_s {
  pq_t * coa;        /* Coarse product quantizer */
  pq_t * resi;       /* Product quantizer for the residual vector */
} pqr_t;


/* create (learn) /free the residual quantizer */
pqr_t * pqr_learn (int d, int nsq_coa, int nsq_bits_coa, int nsq_resi, int nsq_bits_resi, 
		   const real * v, int n, int nb_iter_max);

void pqr_free (pqr_t * pqr);



/*--- encode/decode product quantizer representations ---*/

/* encode the input vectors */
void pqr_encode (const pqr_t * pqr, int n, const float * v, 
		 pqcode_t * codes_coa_out, pqcode_t * codes_resi_out);


/* decode the input codes. The residual vector after coarse quantization 
   must be allocated externaly (size=n * d * sizeof (float))*/
void pqr_decode (const pqr_t * pqr, int n, const pqcode_t * codes_coa, 
		 const pqcode_t * codes_resi, float * v_out);


/*--- nearest neighbor search (asymetric version, unbiased) */

/* find the k-nearest neighbor using refinement for kver of the nearest neighbors */
void pqr_nns (const pqr_t * pqr, const real * vquery, 
	      const pqcode_t * codes_coa, const pqcode_t * codes_resi, 
	      int n1, int n2, int k, int kverif, 
	      int * nns_out, real * nns_dis_out);



#endif
