#include <stdlib.h> 
#include <assert.h>

#include <yael/sorting.h>
#include <yael/nn.h>
#include <yael/vector.h>

#include "pq.h"
#include "pqr.h"

/* create (learn) /free the residual quantizer */
pqr_t * pqr_learn (int d, int nsq_coa, int nsq_bits_coa, int nsq_resi, int nsq_bits_resi, 
		   const real * v, int n, int nb_iter_max)
{
  pqr_t * pqr = (pqr_t *) malloc (sizeof (*pqr));

  printf ("pqr_learn: learn coarse product quantizer\n");
  pqr->coa = pq_learn (d, nsq_coa, nsq_bits_coa, v, n, nb_iter_max);

  /* compute the residual and learn the refinement quantizer on it */
  real * resi = malloc (d * n * sizeof (*resi));
  pq_residual (pqr->coa, v, n, resi);


  printf ("pqr_learn: learn the fine quantizer\n");
  pqr->resi = pq_learn (d, nsq_resi, nsq_bits_resi, resi, n, nb_iter_max);
  free (resi);
  
  return pqr;
}

void pqr_free (pqr_t * pqr)
{
  pq_free (pqr->coa);
  pq_free (pqr->resi);
  free (pqr);
}


/* encode the input vectors */
void pqr_encode (const pqr_t * pqr, int n, const float * v, 
		 pqcode_t * codes_coa_out, pqcode_t * codes_resi_out)
{
  int d = pqr->coa->d;
  pq_encode (pqr->coa, v, n, codes_coa_out);

  /* compute and encode the residual */
  real * resi = malloc (d * n * sizeof (*resi));
  pq_residual_from_codes (pqr->coa, v, codes_coa_out, n, resi);

  pq_encode (pqr->resi, resi, n, codes_resi_out);
  free (resi);
}


/* decode the input codes. 
   The residual vector (of dimension d) must be allocated externaly */
void pqr_decode (const pqr_t * pqr, int n, const pqcode_t * codes_coa, 
		 const pqcode_t * codes_resi, float * v_out)
{
  int d = pqr->coa->d;
  int i;
  real * resi = (real *) malloc (n * d * sizeof (*resi));
  pq_decode (pqr->coa, codes_coa, n, v_out);
  pq_decode (pqr->resi, codes_resi, n, resi);

  for (i = 0 ; i < n * d ; i++)
    v_out[i] += resi[i];
 
  free (resi);
}


/*--- distance computation and nearest neighbor search (asymetric version, unbiased) */

/* compute the estimate distance. kver is the number 
   of neighbors for which the fine quantizer is used to refine the distance. */

void pqr_nns (const pqr_t * pqr, const real * vquery, 
	      const pqcode_t * codes_coa, const pqcode_t * codes_resi,
	      int n1, int n2, int k, int kverif,
	      int * nns, real * nns_dis)
{
  assert (kverif >= k && kverif <= n2);

  int i, j, d =  pqr->coa->d;
  int off_coa = pqcode_sizeof (pqr->coa);
  int off_resi = pqcode_sizeof (pqr->resi);

  real * dis = (real *) malloc (n2 * sizeof (*dis));

  /* vectors to be refined */
  int * idxtorefine = (int *) malloc (kverif * sizeof (*idxtorefine));
  real * dec = (real *) malloc (d * sizeof (*dec));
  real * disfine = (real *) malloc (kverif * sizeof (*disfine));

  for (i = 0 ; i < n1 ; i++) {
    /* find kverif closest codes with the coarse PQ */
    pq_asym_L2sq (pqr->coa, vquery + i * d, codes_coa, n2, dis);
    fvec_find_k_min (dis, n2, idxtorefine, kverif);

    /* refine only those kverif distances */
    for (j = 0 ; j < kverif ; j++) {
      int off1 = idxtorefine[j] * off_coa;
      int off2 = idxtorefine[j] * off_resi;
      pqr_decode (pqr, 1, codes_coa + off1, codes_resi + off2, dec);
      disfine[j] = fvec_distance_L2sqr (vquery + i * d, dec, d);
    }

    /* find k closest points out of kverif */ 
    fvec_find_k_min (disfine, kverif, nns + i * k, k);
    
    for (j = 0 ; j < k ; j++)
      nns[i * k + j] = idxtorefine[nns[i * k + j]];

    if (nns_dis != NULL)
      for (j = 0; j < k; j++)
        nns_dis[i * k + j] = dis[nns[i * k + j]];

  }

  free (dis);
  free (idxtorefine);
  free (dec);
  free (disfine);    
}
