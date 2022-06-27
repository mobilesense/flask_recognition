#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include <yael/machinedeps.h>
#include <yael/sorting.h>
#include <yael/vector.h>
#include <yael/clustering.h>
#include <yael/nn.h>

#include "pq.h"


/*--- Declaration of internal types ---*/

typedef unsigned short int pqcode_std_t;     /* Standard internal representation of codes */
typedef unsigned long long int pqcode8x8_t;  /* 8x8 code representation */


inline long pqcode_sizeof (const pq_t * pq)
{
  if (pq->mode == PQMODE_STD)
    return sizeof (pqcode_std_t) * pq->nsq;
  else if (pq->mode == PQMODE_8x8)
    return sizeof (pqcode8x8_t);
  else
    assert (0);
  return 0;
}


int pqcode_sizeof_packed (const pq_t * pq)
{
  return (pq->nsq * pq->nsq_bits + 7) / 8;
}


pqcode_t *pqcode_new (const pq_t * pq, int n)
{
  pqcode_t * codes = (pqcode_t *) memalign (16, n * pqcode_sizeof (pq));
  return codes; 
}


pqcode_t *pqcode_realloc (const pq_t * pq, pqcode_t * code, int n)
{
  return realloc (code, pqcode_sizeof (pq) * n);
}


pqcode_t *pqcode_offset (const pq_t * pq, const pqcode_t * code, int i)
{
  return (pqcode_t *) ((char *) code + i * pqcode_sizeof(pq));
}


void pqcode_free (pqcode_t * codes)
{
  free (codes);
}


void pqcode_display (const pq_t * pq, const pqcode_t * codes, int n)
{
  int i, q;
  
  if (pq->mode == PQMODE_STD)
    for (i = 0 ; i < n ; i++) {
      pqcode_std_t * c = (pqcode_std_t *) codes;
      for (q = 0 ; q < pq->nsq ; q++)
	printf ("%d ", c[i * pq->nsq + q]);
      printf (" / ");
    }

  else if (pq->mode == PQMODE_8x8) {
    for (i = 0 ; i < n ; i++) {
      pqcode8x8_t c = ((pqcode8x8_t *) codes)[i];
      for (q = 0 ; q < pq->nsq ; q++) {
	printf ("%d ", (int) (c & 0xff));
	c >>= 8;
      }
      printf (" / ");
    }
  }

  else assert (0);
}


/* The two next functions serve for the I/O of pqcodes. 
   It might be intesting to export it, but I see no reason at the ime  */
static unsigned char * pqcode_pack (const pq_t * pq, const pqcode_t * codes, 
				    int n, int * nbytes_out)
{
  /* compute the number of bits/bytes required to store the codes */
  int nsq = pq->nsq;
  int nsq_bits = pq->nsq_bits;
  int packed_size = pqcode_sizeof_packed (pq);

  /* in optimized modes, the memory representation is let unchanged */
  if (pq->mode != PQMODE_STD) {
    *nbytes_out = n * pqcode_sizeof (pq);
    unsigned char * packed = (unsigned char *) malloc (*nbytes_out);
    memcpy (packed, codes, *nbytes_out);
    return packed;
  }

  unsigned char * packed = (unsigned char *) calloc (packed_size * n, sizeof (*packed));

  int i, q, j, bitpos;

  for (i = 0 ; i < n ; i++) {
    bitpos = packed_size * i * 8;

    for (q = 0 ; q < nsq ; q++) {
      pqcode_std_t code = ((pqcode_std_t *) codes)[i * nsq + q];

      for (j = 0 ; j < nsq_bits ; j++) {
	int bit = (code >> j) & 1;              /* the bit to write */
	int bytepos = bitpos / 8;
	int bitmod = bitpos % 8;
	
	packed[bytepos] |= (unsigned char) (bit << bitmod);
	bitpos++;
      }
    }
  }
  *nbytes_out = n * packed_size;
  return packed;
}


static pqcode_t * pqcode_unpack (const pq_t * pq, const unsigned char * packed, 
				 int nbytes, int * n_out)
{
  /* compute the number of bits/bytes required to store the codes */
  int nsq = pq->nsq;
  int nsq_bits = pq->nsq_bits;
  int packed_size = pqcode_sizeof_packed (pq);

  /* in optimized modes, the memory representation is let unchanged */
  if (pq->mode != PQMODE_STD) {
    pqcode_t * codes = (pqcode_t *) malloc (nbytes);
    memcpy (codes, packed, nbytes);
    *n_out = nbytes / pqcode_sizeof (pq);
    return codes;
  }

  /* this is the number for symbols to decode */
  int n = nbytes / packed_size;          
  assert (n * packed_size == nbytes);

  pqcode_std_t * codes = (pqcode_std_t *) pqcode_new (pq, n);

  int i, q, j, bitpos;

  for (i = 0 ; i < n ; i++) {
    bitpos = packed_size * i * 8;

    for (q = 0 ; q < nsq ; q++) {
      pqcode_std_t code = 0;

      for (j = 0 ; j < nsq_bits ; j++) {
	int bytepos = bitpos / 8;
	int bitmod = bitpos % 8;
	int bit = (packed[bytepos] >> bitmod) & 1;

	code |= bit << j;
	bitpos++;
      }
      codes[i * nsq + q] = code;
    }
  }

  *n_out = n;
  return (pqcode_t *) codes;
}



void pqcode_fwrite (FILE * f, const pq_t * pq, const pqcode_t * codes, int n)
{
  int nbytes;
  unsigned char * packed = pqcode_pack (pq, codes, n, &nbytes);
  
  fwrite (packed, nbytes, 1, f);
  free (packed);
}


void pqcode_write (const char * filename, const pq_t * pq, const pqcode_t * codes, int n)
{
  FILE * f = fopen (filename, "w");

  pqcode_fwrite (f, pq, codes, n);
  fclose (f);
}


pqcode_t * pqcode_new_fread (FILE * f, const pq_t * pq, int n)
{
  int nbytes = n * pqcode_sizeof_packed (pq);
  int n_check;
  unsigned char * packed = (unsigned char *) malloc (nbytes);

  int ret = fread (packed, 1, nbytes, f);
  if (ret != nbytes) {
    fprintf (stderr, "# pcode_new_fread error : unable to read from input file\n");
    free (packed);
    return NULL;
  }

  pqcode_t * codes = pqcode_unpack (pq, packed, nbytes, &n_check);
  assert (nbytes == n_check * pqcode_sizeof_packed (pq));

  free (packed);
  return codes;
}


pqcode_t * pqcode_new_read (const char * filename, const pq_t * pq)
{
  int ret;

  int packed_size = pqcode_sizeof_packed (pq);
  int ncodes_alloc = 100; 
  int ncodes = 0;
  
  unsigned char * buf = malloc (packed_size);
  unsigned char * packed = (unsigned char *) malloc (ncodes_alloc * packed_size);
  
  
  FILE * f = fopen (filename, "r");
  if (!f) {
    fprintf (stderr, "# pcode_new_read error 1: unable to read file %s\n", filename);
    free (buf);
    free (packed);
    return NULL;
  }
  

  while (1) {
    ret = fread (buf, 1, packed_size, f);
    if (ret != packed_size) {
      if (feof (f))
	break;
      else {
	fprintf (stderr, "# pcode_new_read error 2: file %s opened, but incorrectly read\n", filename);
	free (buf);
	free (packed);
	return NULL;
      }
    }
    
    if (ncodes_alloc == ncodes) {                     /* need to re-alloc ? */
      ncodes_alloc = (10 + ncodes_alloc * 14) / 10;   /* geometric re-allocation with factor 1.4 */
      packed = realloc (packed, ncodes_alloc * packed_size);
      assert (packed);
    }

    memcpy (packed + packed_size * ncodes, buf, packed_size);
    ncodes++;
  }

  int n_check;
  pqcode_t * codes = pqcode_unpack (pq, packed, ncodes * packed_size, &n_check);
  assert (ncodes == n_check);
  free (packed);
  free (buf);
  fclose (f);

  return codes;
}


/*--- General functions ----*/

pq_t *pq_new (int d, int nsq, int nsq_bits)
{
  int i, q;
  pq_t *pq = (pq_t *) malloc (sizeof (pq_t));
  pq->d = d;
  pq->nsq = nsq;
  assert (pq->d % pq->nsq == 0);        /* dimension must be a power of nb nsq */

  pq->nsq_bits = nsq_bits;
  pq->ksq = 1 << pq->nsq_bits;

  pq->centroids = malloc (nsq * sizeof (*pq->centroids));
  pq->centroids_dis = malloc (nsq * sizeof (*pq->centroids_dis));
  pq->sberr = malloc (nsq * sizeof (*pq->sberr));
  for (q = 0; q < pq->nsq; q++) {
    pq->centroids[q] = NULL;
    pq->centroids_dis[q] = (real *) calloc (pq->ksq * pq->ksq,
                                            sizeof (**pq->centroids_dis));
    pq->sberr[q] = (real *) calloc (pq->ksq, sizeof (**pq->sberr));
  }

  pq->mode = PQMODE_STD;

  /* by default, use the normal order for the components */
  pq->extract_map = (int *) malloc (pq->d * sizeof (*pq->extract_map));
  for (i = 0 ; i < pq->d ; i++)
    pq->extract_map[i] = i;

  return pq;
}


void pq_free (pq_t * pq)
{
  int q;
  free (pq->centroids[0]);
  for (q = 0; q < pq->nsq; q++) {
    free (pq->centroids_dis[q]);
    free (pq->sberr[q]);
  }
  free (pq->centroids);
  free (pq->centroids_dis);
  free (pq->sberr);
  free (pq);
}


/* This static function is used to select which components are aggregated */
static inline void pq_extract_subvectors (const pq_t * pq, const real * v, int n,
					  int subset, real * x_out)
{
  int i, j;
  int dsq = pq->d / pq->nsq;
  const real *buf = v;
  real * x = x_out;

  for (i = 0 ; i < n ; i++) {
    for (j = 0 ; j < dsq ; j++) 
      x[j] = buf[pq->extract_map[j + dsq * subset]];
    
    buf += pq->d;
    x += dsq;
  }
}




pq_t *pq_learn (int d, int nsq, int nsq_bits,
                const real * v, int n, int nb_iter_max) {
  
  return pq_learn_with_extract_map (d, nsq, nsq_bits, v, NULL, n, nb_iter_max); 
}


static real  pq_learn_1_quantizer(pq_t *pq, const float *v, int n, int nb_iter_max, int n_thread, 
                                 int q) {

  real toterr = 0;
  int i, c;                  /* i: for v, c: for centroids */

  int dsq = pq->d / pq->nsq;
  int ksq = pq->ksq;
  
  real *x = (real *) malloc (n * dsq * sizeof (*x));

  /* first create the set of lower dimensional vectors */
  pq_extract_subvectors (pq, v, n, q, x);

  /* k-means clustering */
  int *clust_assign;
  fprintf (stderr, "                                            "
           "                  q=%2d/%d\r", q + 1, pq->nsq);
  
  real * cents=clustering_kmeans_assign_with_score (n, dsq, x, ksq,
                                                    nb_iter_max, 0, n_thread,
                                                    NULL, &clust_assign);
  
  pq->centroids[q] = pq->centroids[0] + q * dsq * ksq;

  memcpy(pq->centroids[q], cents, sizeof(real) * dsq * ksq);

  free(cents);

  /* count the number of vectors assigned to each cluster */
  int *nb_assign = ivec_new_histogram (ksq, clust_assign, n);

  /* compute the quantization error of a vector with its centroids */
  real **qerr = malloc (ksq * sizeof (*qerr));
  int *posinqerr = calloc (ksq, sizeof (*posinqerr));
  for (c = 0; c < ksq; c++)
    qerr[c] = calloc (nb_assign[c], sizeof (**qerr));

  for (i = 0; i < n; i++) {
    c = clust_assign[i];
    real *centroid = pq->centroids[q] + c * dsq;
    real err = fvec_distance_L2sqr (x + i * dsq, centroid, dsq);
    qerr[c][posinqerr[c]++] += err;
    toterr += err;
  }

  for (c = 0; c < ksq; c++) {
    for (i = 0; i < nb_assign[c]; i++)
      pq->sberr[q][c] += qerr[c][i];
    pq->sberr[q][c] = pq->sberr[q][c] / nb_assign[c]; /* biased estimation */
    free (qerr[c]);
  }

  free (posinqerr);
  free (qerr);
  free (clust_assign);
  free (nb_assign);
  
  /* compute the distance between the centroids */
  compute_cross_distances (dsq, ksq, ksq,
                           pq->centroids[q], pq->centroids[q],
                           pq->centroids_dis[q]);
  
  free(x);

  return toterr;
} 

typedef struct {
  pq_t *pq; 
  const float *v;
  int n;
  int nb_iter_max;
  int n_thread;

  float *toterrs;  
} pq_learn_t;  


static void task_pq_learn_1_quantizer (void *arg, int tid, int q) {
  pq_learn_t *l=arg;
  l->toterrs[q]=pq_learn_1_quantizer(l->pq, l->v, l->n, l->nb_iter_max, l->n_thread, q);
}

pq_t *pq_learn_with_extract_map (int d, int nsq, int nsq_bits,
                                 const real * v, const int *extract_map,
                                 int n, int nb_iter_max)
{
  int q;                  /* i: for v, c: for centroids, q: for quantizer */

  pq_t *pq = pq_new (d, nsq, nsq_bits);

  if(extract_map) memcpy(pq->extract_map,extract_map,sizeof(*extract_map)*d);

  int dsq = pq->d / pq->nsq;
  int ksq = pq->ksq;

  printf ("ksq = %d, d = %d, dsq = %d\n", ksq, d, dsq);

  pq->centroids[0]=fvec_new(pq->nsq * dsq * ksq);

  real toterr;
  int ncpu=count_cpu();

  if(2*ncpu>pq->nsq) { /* fine-grain threading */

    toterr=0;
    for (q = 0; q < pq->nsq; q++) 
      toterr+=pq_learn_1_quantizer(pq,v,n,nb_iter_max,ncpu,q);

  } else { /* coarse */

    pq_learn_t pql={pq,v,n,nb_iter_max,1 | CLUSTERING_KMEANS_QUIET,fvec_new(pq->nsq)};

    compute_tasks (pq->nsq,count_cpu(),&task_pq_learn_1_quantizer,&pql);

    toterr=fvec_sum(pql.toterrs,pq->nsq);
    free(pql.toterrs);

  }
  
  fprintf (stdout, "\ntoterr = %8g\n", toterr / n);


  if (nsq == 8 && nsq_bits == 8)
    pq->mode = PQMODE_8x8;

  return pq;
}


void pq_display (const pq_t * pq, int verbose)
{
  int q, i;
  assert (pq);
  printf ("d        (dimension of original vector) = %d\n", pq->d);
  printf ("nsq      (number of subquantizers)      = %d\n", pq->nsq);
  printf ("nsq_bits (number of bits per subquant)  = %d\n", pq->nsq_bits);
  printf ("n_bits   (total bit vector length)      = %d\n",
          pq->nsq * pq->nsq_bits);
  printf ("ksq      (nb centroid - subquantizer)   = %d\n", pq->ksq);

  if (verbose >= 1) {
    for (q = 0; q < pq->nsq; q++) {
      printf ("[---- q = %d ----]\n", q);
      for (i = 0; i < pq->ksq; i++)
        printf ("%-6.1f ", pq->sberr[q][i]);
      printf ("\n");
    }
  }
}


/*--- encode/decode product quantizer representations ---*/
void pq_encode_std (const pq_t * pq, const real * v, int n,
		    pqcode_std_t * codes)
{
  int i, q;
  int nsq = pq->nsq;
  int dsq = pq->d / nsq;


  /* first allocate the number of bits needed to store the indexes */
  int *codes_tmp = malloc (n * sizeof (*codes_tmp));

  real *x = (real *) malloc (n * dsq * sizeof (*x));

  /* re-organize the data */
  for (q = 0; q < nsq; q++) {

    /* first create the set of lower dimensional vectors */
    pq_extract_subvectors (pq, v, n, q, x);

    quantize_codebook (n, pq->ksq, dsq, pq->centroids[q],
                       x, codes_tmp, NULL, NULL);

    for (i = 0 ; i < n ; i++)
      codes[i * nsq + q] = codes_tmp[i];
  }

  free (codes_tmp);
  free (x);
}


void pq_encode_8x8 (const pq_t * pq, const float *v,
		    int n, pqcode8x8_t * codes)
{
  int i, q;
  pqcode_std_t *idx = malloc (8 * sizeof (*idx) * n);
  pq_encode_std (pq, v, n, (pqcode_std_t *) idx);

  assert (pq->nsq == 8 && pq->nsq_bits == 8);
  memset (codes, 0, n * sizeof (*codes));

  for (i = 0; i < n; i++) {
    for (q = pq->nsq - 1; q >= 0; q--) {
      codes[i] <<= 8;
      codes[i] += idx[i * pq->nsq + q];
    }
  }
  free (idx);
}


void pq_encode (const pq_t * pq, const real * v, int n, pqcode_t * codes)
{
  if (pq->mode == PQMODE_STD)
    pq_encode_std (pq, v, n, (pqcode_std_t *) codes);

  else if (pq->mode == PQMODE_8x8)
    pq_encode_8x8 (pq, v, n, (pqcode8x8_t *) codes);

  else
    assert (0);
}


void pq_residual_from_codes (const pq_t * pq, const real * v, const pqcode_t * codes, 
			     int n, real * resi)
{
  int i;
  pq_decode (pq, codes, n, resi);
  for (i = 0 ; i < n * pq->d ; i++)
    resi[i] = v[i] - resi[i];
}


void pq_residual (const pq_t * pq, const real * v, int n, real * resi)
{
  pqcode_t * codes = pqcode_new (pq, n);
  pq_encode (pq, v, n, codes);
  
  pq_residual_from_codes (pq, v, codes, n, resi);
  pqcode_free (codes);
}


static inline void pq_decode_std (const pq_t * pq, const pqcode_std_t * codes,
                                  int n, real * v)
{
  int q, i, j;
  int d = pq->d;
  int nsq = pq->nsq;
  int dsq = d / nsq;

  for (i = 0; i < n; i++)
    for (q = 0; q < nsq; q++)
      for (j = 0; j < dsq; j++)
        v[i * d + q * dsq + j] =
            pq->centroids[q][j + codes[i * nsq + q] * dsq];
}


static inline void pq_decode_8x8 (const pq_t * pq, const pqcode8x8_t * codesb,
                                  int n, real * v)
{
  int q, i, j;
  int d = pq->d;
  int nsq = pq->nsq;
  int dsq = d / nsq;
  pqcode8x8_t codeb, codetmp;

  for (i = 0; i < n; i++) {
    codeb = codesb[i];
    for (q = 0; q < nsq; q++) {
      codetmp = codeb & 0xffLL;
      codeb >>= 8;
      for (j = 0; j < dsq; j++)
        v[i * d + q * dsq + j] = pq->centroids[q][j + codetmp * dsq];
    }
  }
}


void pq_decode (const pq_t * pq, const pqcode_t * codes, int n, real * v)
{
  if (pq->mode == PQMODE_STD)
    pq_decode_std (pq, (pqcode_std_t *) codes, n, v);
  else if (pq->mode == PQMODE_8x8)
    pq_decode_8x8 (pq, (pqcode8x8_t *) codes, n, v);
  else
    assert (0);
}


/* add a code to a float vector */
static inline void pq_fvec_add_std (const pq_t * pq, float *v,
                                    const pqcode_std_t * code)
{
  int q, j, i = 0;
  int d = pq->d;
  int nsq = pq->nsq;
  int dsq = d / nsq;

  for (q = 0; q < nsq; q++)
    for (j = 0; j < dsq; j++)
      v[i++] += pq->centroids[q][j + code[q] * dsq];
}


static inline void pq_fvec_add_8x8 (const pq_t * pq, float *v,
                                    const pqcode8x8_t * code)
{
  int q, j;
  int d = pq->d;
  int nsq = pq->nsq;
  int dsq = d / nsq;

  pqcode8x8_t c = *code;

  for (q = 0; q < nsq; q++) {
    for (j = 0; j < dsq; j++)
      v[q * dsq + j] += pq->centroids[q][j + (c & 0xffLL) * dsq];
    c >>= 8;
  }
}


void pq_fvec_add (const pq_t * pq, float *v, const pqcode_t * code)
{
  if (pq->mode == PQMODE_STD)
    pq_fvec_add_std (pq, v, (pqcode_std_t *) code);

  else if (pq->mode == PQMODE_8x8)
    pq_fvec_add_8x8 (pq, v, (pqcode8x8_t *) code);

  else
    assert (0);
}


/*--- Expected distance and search: the symmetric case ----*/


static inline void pq_L2sq_cor_std (const pq_t * pq, const pqcode_std_t * codes1,
                                    const pqcode_std_t * codes2, int n,
                                    real * dis)
{
  int i, q;
  memset (dis, 0, n * sizeof (*dis));

  for (i = 0; i < n; i++)
    for (q = 0; q < pq->nsq; q++)
      dis[i] += pq->sberr[q][codes1[q]]
          + pq->sberr[q][codes2[q + i * pq->nsq]]
          + pq->centroids_dis[q][codes1[q] * pq->ksq +
                                 codes2[q + i * pq->nsq]];
}


void pq_L2sq_cor (const pq_t * pq, const pqcode_t * codes1,
                  const pqcode_t * codes2, int n, real * dis)
{
  if (pq->mode == PQMODE_STD)
    pq_L2sq_cor_std (pq, (pqcode_std_t *) codes1, (pqcode_std_t *) codes2, n, dis);

  else if (pq->mode == PQMODE_8x8) {
    fprintf (stderr,
             "# function pq_L2sq_cor not implemented in optimized mode\n");
    exit (0);
  } else
    assert (0);
}


static inline void pq_L2sq_std (const pq_t * pq, const pqcode_std_t * codes1,
                                const pqcode_std_t * codes2, int n, real * dis)
{
  int i, q, tmp;
  memset (dis, 0, n * sizeof (*dis));

  for (i = 0; i < n; i++) {
    tmp = i * pq->nsq;

    for (q = 0; q < pq->nsq; q++)
      dis[i] += pq->centroids_dis[q][codes1[q] * pq->ksq + codes2[q + tmp]];
  }
}


static inline void pq_L2sq_8x8 (const pq_t * pq, const pqcode8x8_t * codeq,
                                const pqcode8x8_t * codesb, int n, real * dis)
{
  int i;
  real distmp;
  pqcode8x8_t cq = *codeq;

  real *cdis0 = pq->centroids_dis[0] + ((cq >> 0) & 0xffLL) * 256;
  real *cdis1 = pq->centroids_dis[1] + ((cq >> 8) & 0xffLL) * 256;
  real *cdis2 = pq->centroids_dis[2] + ((cq >> 16) & 0xffLL) * 256;
  real *cdis3 = pq->centroids_dis[3] + ((cq >> 24) & 0xffLL) * 256;
  real *cdis4 = pq->centroids_dis[4] + ((cq >> 32) & 0xffLL) * 256;
  real *cdis5 = pq->centroids_dis[5] + ((cq >> 40) & 0xffLL) * 256;
  real *cdis6 = pq->centroids_dis[6] + ((cq >> 48) & 0xffLL) * 256;
  real *cdis7 = pq->centroids_dis[7] + ((cq >> 56) & 0xffLL) * 256;

  pqcode8x8_t codeb;
  assert (pq->nsq == 8 && pq->nsq_bits == 8);
  memset (dis, 0, n * sizeof (*dis));

  for (i = 0; i < n; i++) {
    codeb = codesb[i];
    distmp = cdis0[codeb & 0xffLL];
    codeb >>= 8;
    distmp += cdis1[codeb & 0xffLL];
    codeb >>= 8;
    distmp += cdis2[codeb & 0xffLL];
    codeb >>= 8;
    distmp += cdis3[codeb & 0xffLL];
    codeb >>= 8;
    distmp += cdis4[codeb & 0xffLL];
    codeb >>= 8;
    distmp += cdis5[codeb & 0xffLL];
    codeb >>= 8;
    distmp += cdis6[codeb & 0xffLL];
    codeb >>= 8;
    dis[i] = distmp + cdis7[codeb & 0xffLL];
  }
}


void pq_L2sq (const pq_t * pq, const pqcode_t * code1,
              const pqcode_t * codes2, int n, real * dis)
{
  if (pq->mode == PQMODE_STD)
    pq_L2sq_std (pq, (pqcode_std_t *) code1, (pqcode_std_t *) codes2, n, dis);
  else if (pq->mode == PQMODE_8x8)
    pq_L2sq_8x8 (pq, (pqcode8x8_t *) code1, (pqcode8x8_t *) codes2, n, dis);
  else
    assert (0);
}


void pq_nns_cor (const pq_t * pq, const pqcode_t * codes1,
                 const pqcode_t * codes2, int n1, int n2, int k,
                 int *nns, real * nns_dis)
{
  int i, j;
  int offset = pqcode_sizeof (pq);
  real *dis = (real *) malloc (n2 * sizeof (*dis));

  for (i = 0; i < n1; i++) {
    pq_L2sq_cor (pq, codes1 + i * offset, codes2, n2, dis);
    fvec_find_k_min (dis, n2, nns + i * k, k);

    if (nns_dis != NULL)
      for (j = 0; j < k; j++)
        nns_dis[i * k + j] = dis[nns[i * k + j]];
  }

  free (dis);
}


void pq_nns (const pq_t * pq, const pqcode_t * codes1,
             const pqcode_t * codes2, int n1, int n2, int k,
             int *nns, real * nns_dis)
{
  int i, j;
  int offset = pqcode_sizeof (pq);
  real *dis = (real *) malloc (n2 * sizeof (*dis));

  for (i = 0; i < n1; i++) {
    pq_L2sq (pq, codes1 + i * offset, codes2, n2, dis);

    fvec_find_k_min (dis, n2, nns + i * k, k);

    if (nns_dis != NULL)
      for (j = 0; j < k; j++)
        nns_dis[i * k + j] = dis[nns[i * k + j]];
  }

  free (dis);
}



/*---- Expected distance and search: asymetric case (query not quantized) */


/* compute the expected distance between two points */
void pq_asym_L2sq_cor_std (const pq_t * pq, const real * v,
                           const pqcode_std_t * codes, int n, real * dis)
{
  int i, c, q;
  int dsq = pq->d / pq->nsq;
  int ksq = pq->ksq;

  real *discentroids = malloc (ksq * sizeof (*discentroids));
  real *x = (real *) malloc (dsq * sizeof (*x));

  memset (dis, 0, n * sizeof (*dis));

  for (q = 0; q < pq->nsq; q++) {

    /* compute the distance from v to all the possible centroids */
    pq_extract_subvectors (pq, v, 1, q, x);

    for (c = 0; c < ksq; c++) {
      real *centroid = pq->centroids[q] + c * dsq;
      discentroids[c] = fvec_distance_L2sqr (x, centroid, dsq);
    }

    for (i = 0; i < n; i++) {
      int c = codes[q + i * pq->nsq];
      dis[i] += discentroids[c] + pq->sberr[q][c];
    }
  }

  free (discentroids);
  free (x);
}


void pq_asym_L2sq_cor (const pq_t * pq, const real * v,
                       const pqcode_t * codes, int n, real * dis)
{
  if (pq->mode == PQMODE_STD)
    pq_asym_L2sq_cor_std (pq, v, (pqcode_std_t *) codes, n, dis);
  else if (pq->mode == PQMODE_8x8) {
    fprintf (stderr,
             "# function pq_asym_L2sq_cor_8x8 not implemented in optimized mode\n");
    exit (0);
  } else
    assert (0);
}



static inline void pq_asym_L2sq_std (const pq_t * pq, const real * vquery,
				     const pqcode_std_t * codesb, int n, real * dis)
{
  int i, c, q;
  int dsq = pq->d / pq->nsq;
  int ksq = pq->ksq;

  real *discentroids = malloc (pq->nsq * ksq * sizeof (*discentroids));
  real *x = (real *) malloc (dsq * sizeof (*x));

  memset (dis, 0, n * sizeof (*dis));

  for (q = 0; q < pq->nsq; q++) {

    /* compute the distance from v to all the possible centroids */
    pq_extract_subvectors (pq, vquery, 1, q, x);

    for (c = 0; c < ksq; c++) {
      real *centroid = pq->centroids[q] + c * dsq;
      discentroids[c + ksq * q] = fvec_distance_L2sqr (x, centroid, dsq);
    }
  }

  for (i = 0; i < n; i++) {
    for (q = 0; q < pq->nsq; q++) {
      int c = codesb[q + i * pq->nsq];
      dis[i] += discentroids[c + ksq * q];
    }
  }

  free (discentroids);
  free (x);
}


static inline void pq_asym_L2sq_8x8 (const pq_t * pq, const real * v,
                                     const pqcode8x8_t * codesb, int n,
                                     real * dis)
{
  int i, c, q;
  int dsq = pq->d / pq->nsq;
  int ksq = pq->ksq;
  real distmp;

  real *discentroids = malloc (pq->nsq * ksq * sizeof (*discentroids));
  real *x = (real *) malloc (dsq * sizeof (*x));

  real *cdis0 = discentroids;
  real *cdis1 = discentroids + ksq * 1;
  real *cdis2 = discentroids + ksq * 2;
  real *cdis3 = discentroids + ksq * 3;
  real *cdis4 = discentroids + ksq * 4;
  real *cdis5 = discentroids + ksq * 5;
  real *cdis6 = discentroids + ksq * 6;
  real *cdis7 = discentroids + ksq * 7;

  pqcode8x8_t codeb;
  assert (pq->nsq == 8 && pq->nsq_bits == 8);


  


  for (q = 0; q < pq->nsq; q++) {

    /* compute the distance from v to all the possible centroids */
    pq_extract_subvectors (pq, v, 1, q, x);

    for (c = 0; c < ksq; c++) {
      real *centroid = pq->centroids[q] + c * dsq;
      discentroids[c + ksq * q] = fvec_distance_L2sqr (x, centroid, dsq);
    }

    /*    compute_cross_distances(dsq, 1, ksq, x, pq->centroids[q], discentroids+q*ksq); */

  }



  for (i = 0; i < n; i++) {
    codeb = codesb[i];
    distmp = cdis0[codeb & 0xffLL];
    codeb >>= 8;
    distmp += cdis1[codeb & 0xffLL];
    codeb >>= 8;
    distmp += cdis2[codeb & 0xffLL];
    codeb >>= 8;
    distmp += cdis3[codeb & 0xffLL];
    codeb >>= 8;
    distmp += cdis4[codeb & 0xffLL];
    codeb >>= 8;
    distmp += cdis5[codeb & 0xffLL];
    codeb >>= 8;
    distmp += cdis6[codeb & 0xffLL];
    codeb >>= 8;
    dis[i] = distmp + cdis7[codeb & 0xffLL];
  }

  free (discentroids);
  free (x);
}


void pq_asym_L2sq (const pq_t * pq, const real * vquery,
                   const pqcode_t * codesb, int n, real * dis)
{
  if (pq->mode == PQMODE_STD)
    pq_asym_L2sq_std (pq, vquery, (pqcode_std_t *) codesb, n, dis);
  else if (pq->mode == PQMODE_8x8)
    pq_asym_L2sq_8x8 (pq, vquery, (pqcode8x8_t *) codesb, n, dis);
  else
    assert (0);
}




void pq_asym_nns_cor (const pq_t * pq,
                      const real * v, const pqcode_t * codes2,
                      int n1, int n2, int k, int *nns, real * nns_dis)
{
  int i, j;
  real *dis = (real *) malloc (n2 * sizeof (*dis));

  for (i = 0; i < n1; i++) {
    pq_asym_L2sq_cor (pq, v + i * pq->d, codes2, n2, dis);
    fvec_find_k_min (dis, n2, nns + i * k, k);

    if (nns_dis != NULL)
      for (j = 0; j < k; j++)
        nns_dis[i * k + j] = dis[nns[i * k + j]];
  }

  free (dis);
}


void pq_asym_nns (const pq_t * pq,
                  const real * vquery, const pqcode_t * codesb,
                  int n1, int n2, int k, int *nns, real * nns_dis)
{
  int i, j;
  real *dis = (real *) malloc (n2 * sizeof (*dis));

  for (i = 0; i < n1; i++) {
    pq_asym_L2sq (pq, vquery + i * pq->d, codesb, n2, dis);
    fvec_find_k_min (dis, n2, nns + i * k, k);

    if (nns_dis != NULL)
      for (j = 0; j < k; j++)
        nns_dis[i * k + j] = dis[nns[i * k + j]];
  }

  free (dis);
}


static inline void pq_asym_fc_nns_std (const pq_t * pq,
                                       const pqcode_std_t * codesq,
                                       const real * vbase, int n1, int n2,
                                       int k, int *nns, real * nns_dis)
{
  int i = 0, j, c, q;
  int nsq = pq->nsq;
  int dsq = pq->d / nsq;
  int ksq = pq->ksq;

  real *dissubcent = malloc (n2 * nsq * ksq * sizeof (*dissubcent));
  real *xb = (real *) malloc (n2 * dsq * sizeof (*xb));
  real *dis = (real *) malloc (n2 * sizeof (*dis));

  for (q = 0; q < nsq; q++) {

    pq_extract_subvectors (pq, vbase, n2, q, xb);

    for (c = 0; c < ksq; c++) {
      real *centroid = pq->centroids[q] + c * dsq;

      for (j = 0; j < n2; j++)
        dissubcent[i++] = fvec_distance_L2sqr (xb + j * dsq, centroid, dsq);
    }
  }

  for (i = 0; i < n1; i++) {
    memset (dis, 0, n2 * sizeof (*dis));

    for (q = 0; q < nsq; q++) {
      int c = codesq[q + i * nsq];
      real *disptr = dissubcent + n2 * (c + ksq * q);

      for (j = 0; j < n2; j++)
        dis[j] += disptr[j];
    }
    fvec_find_k_min (dis, n2, nns + i * k, k);

    if (nns_dis != NULL)
      for (j = 0; j < k; j++)
        nns_dis[i * k + j] = dis[nns[i * k + j]];
  }

  free (dissubcent);
  free (xb);
  free (dis);
}


void pq_asym_fc_nns_8x8 (const pq_t * pq, const pqcode8x8_t * codesq,
                         const real * vbase,
                         int n1, int n2, int k, int *nns, real * nns_dis)
{
  int i = 0, j, c, q;
  int nsq = pq->nsq;
  int dsq = pq->d / nsq;
  int ksq = pq->ksq;

  real *dissubcent = malloc (n2 * nsq * ksq * sizeof (*dissubcent));
  real *xb = (real *) malloc (n2 * dsq * sizeof (*xb));
  real *dis = (real *) malloc (n2 * sizeof (*dis));
  real distmp;

  for (q = 0; q < nsq; q++) {

    pq_extract_subvectors (pq, vbase, n2, q, xb);

    for (c = 0; c < ksq; c++) {
      real *centroid = pq->centroids[q] + c * dsq;

      for (j = 0; j < n2; j++)
        dissubcent[i++] = fvec_distance_L2sqr (xb + j * dsq, centroid, dsq);
    }
  }

  for (i = 0; i < n1; i++) {
    pqcode8x8_t codeq = codesq[i];
    //    memset (dis, 0, n2 * sizeof (*dis));

    real *disptr0 = dissubcent + n2 * ((codeq & 0xffLL) + ksq * 0);
    codeq >>= 8;
    real *disptr1 = dissubcent + n2 * ((codeq & 0xffLL) + ksq * 1);
    codeq >>= 8;
    real *disptr2 = dissubcent + n2 * ((codeq & 0xffLL) + ksq * 2);
    codeq >>= 8;
    real *disptr3 = dissubcent + n2 * ((codeq & 0xffLL) + ksq * 3);
    codeq >>= 8;
    real *disptr4 = dissubcent + n2 * ((codeq & 0xffLL) + ksq * 4);
    codeq >>= 8;
    real *disptr5 = dissubcent + n2 * ((codeq & 0xffLL) + ksq * 5);
    codeq >>= 8;
    real *disptr6 = dissubcent + n2 * ((codeq & 0xffLL) + ksq * 6);
    codeq >>= 8;
    real *disptr7 = dissubcent + n2 * ((codeq & 0xffLL) + ksq * 7);
    codeq >>= 8;

    for (j = 0; j < n2; j++) {
      distmp = disptr0[j];
      distmp += disptr1[j];
      distmp += disptr2[j];
      distmp += disptr3[j];
      distmp += disptr4[j];
      distmp += disptr5[j];
      distmp += disptr6[j];
      dis[j] = distmp + disptr7[j];
    }

    fvec_find_k_min (dis, n2, nns + i * k, k);

    if (nns_dis != NULL)
      for (j = 0; j < k; j++)
        nns_dis[i * k + j] = dis[nns[i * k + j]];
  }

  free (dissubcent);
  free (xb);
  free (dis);
}



void pq_asym_fc_nns (const pq_t * pq, const pqcode_t * codesq,
                     const real * vbase,
                     int n1, int n2, int k, int *nns, real * nns_dis)
{
  if (pq->mode == PQMODE_STD)
    pq_asym_fc_nns_std (pq, (pqcode_std_t *) codesq, vbase, n1, n2, k, nns, nns_dis);
  else if (pq->mode == PQMODE_8x8)
    pq_asym_fc_nns_8x8 (pq, (pqcode8x8_t *) codesq, vbase, n1, n2, k, nns,
                        nns_dis);
  else
    assert (0);
}


/*--- Input/Output ---*/

void pq_fwrite (FILE * f, const pq_t * pq)
{
  int q, ret;
  int dsq = pq->d / pq->nsq;
  int ksq = pq->ksq;
  int pq_headersize = 5 * sizeof (int);


  ret = fwrite (pq, pq_headersize, 1, f);
  assert (ret == 1);

  for (q = 0; q < pq->nsq; q++) {
    ret = fwrite (pq->centroids[q], sizeof (**pq->centroids), dsq * ksq, f);
    assert (ret == dsq * ksq);
  }

  for (q = 0; q < pq->nsq; q++) {
    ret = fwrite (pq->centroids_dis[q], sizeof (**pq->centroids_dis), 
		  ksq * ksq, f);
    assert (ret == pq->ksq * pq->ksq);
  }

  for (q = 0; q < pq->nsq; q++) {
    ret = fwrite (pq->sberr[q], sizeof (**pq->sberr), ksq, f);
    assert (ret == ksq);
  }

  ret = fwrite (pq->extract_map, sizeof (*pq->extract_map), pq->d, f);
  assert (ret == pq->d);
}


void pq_write (const char *fname, const pq_t * pq)
{
  FILE *f = fopen (fname, "w");
  if (!f) {
    fprintf (stderr, "# Unable to open pq file %s\n", fname);
    exit (1);
  }

  pq_fwrite (f, pq);
  fclose (f);
}


pq_t  *pq_fread (FILE * f)
{
  int q, ret;
  int pq_headersize = 5 * sizeof (int);

  pq_t *pq = (pq_t *) malloc (sizeof (pq_t));

  ret = fread (pq, pq_headersize, 1, f);
  assert (ret == 1);
  int dsq = pq->d / pq->nsq;

  pq->centroids = malloc (pq->nsq * sizeof (*pq->centroids));
  pq->centroids_dis = malloc (pq->nsq * sizeof (*pq->centroids_dis));
  pq->sberr = malloc (pq->nsq * sizeof (*pq->sberr));

  float *allcentroids=fvec_new(pq->nsq * dsq * pq->ksq);

  for (q = 0; q < pq->nsq; q++) {
    pq->centroids[q] = allcentroids + q * dsq * pq->ksq;
    ret = fread (pq->centroids[q], sizeof (**pq->centroids), dsq * pq->ksq, f);
    assert (ret == dsq * pq->ksq);
  }

  
  for (q = 0; q < pq->nsq; q++) {
    pq->centroids_dis[q] = (real *) calloc (pq->ksq * pq->ksq,
                                            sizeof (**pq->centroids_dis));
    ret = fread (pq->centroids_dis[q], sizeof (**pq->centroids_dis),
		 pq->ksq * pq->ksq, f);
    assert (ret == pq->ksq * pq->ksq);
  }

  for (q = 0; q < pq->nsq; q++) {
    pq->sberr[q] = (real *) calloc (pq->ksq, sizeof (**pq->sberr));
    ret = fread (pq->sberr[q], sizeof (**pq->sberr), pq->ksq, f);
    assert (ret == pq->ksq);
  }

  pq->extract_map = malloc (sizeof (*pq->extract_map) * pq->d);
  ret = fread (pq->extract_map, sizeof (*pq->extract_map), pq->d, f);
  assert (ret == pq->d);


  return pq;
}


pq_t *pq_read (const char *fname)
{
  pq_t *pq;
  FILE *f = fopen (fname, "r");
  assert (f);
  pq = pq_fread (f);
  fclose (f);
  return pq;
}




/*--- k-means clustering using codes ---*/

float *pq_kmeans_clustering (const pq_t * pq, const pqcode_t * codes,
                             int n, int k, int nb_iter_max,
                             int *clust_assign_out)
{
  int i, j, c, iter = 0;
  int d = pq->d;
  int offset = pqcode_sizeof (pq);

  float *centroids = fvec_new (k * d);
  float *centroids_new = fvec_new (k * d);
  float *tmp;

  /* First select a subset of vectors */
  int *randperm = ivec_new_random_perm (n);
  for (i = 0; i < k; i++)
    pq_decode (pq, codes + randperm[i] * offset, 1, centroids + i * d);
  free (randperm);

  /* vector counting the number of vectors in each centroids */
  int *cluster_size = (int *) malloc (sizeof (int) * n);

  /* to which cluster a vector is assigned ? */
  int *clust_assign = (int *) malloc (sizeof (int) * n);

  /* hash key to detect the convergence */
  double hashconvergence, hashconvergence_new = 1e99;


  fprintf (stderr, "n = %d, k = %d, nbitermax = %d\n", n, k, nb_iter_max);

  /* launch the construction iterating process */
  do {
    iter++;
    printf ("Iteration [%d / %d]\n", iter, nb_iter_max);

    hashconvergence = hashconvergence_new;
    memset (centroids_new, 0, sizeof (float) * k * d);
    memset (cluster_size, 0, sizeof (int) * k);

    pq_asym_fc_nns (pq, codes, centroids, n, k, 1, clust_assign, NULL);

    for (i = 0; i < n; i++) {
      pq_fvec_add (pq, centroids_new + clust_assign[i] * d,
                   codes + i * offset);
      cluster_size[clust_assign[i]]++;
    }

    /* check here because we want clust_assign to be in sync with centroids */
    if (iter == nb_iter_max)
      break;

    tmp = centroids_new;
    for (i = 0; i < k; i++) {
      if (cluster_size[i] > 0)
        for (j = 0; j < d; j++) {
          *tmp /= cluster_size[i];
          tmp++;
        }

      /* re-allocation strategy */
      else {
        /* Find a non-void cluster */
        do
          c = floor (drand48 () * k);
        while (cluster_size[c] < 10);   /* minimum size for splitting a cluster */

        /* Add a very small gaussian noise to the previous cluster */
        float noise_amount = 0.00001;
        {
          float *ci = centroids_new + i * d;
          float *cc = centroids_new + c * d;
          int j;
          for (j = 0; j < d; j++) {
            double r = drand48 ();
            ci[j] = cc[j] + (r - 0.5) * noise_amount;
            cc[j] = cc[j] - (r - 0.5) * noise_amount;
          }
        }
        fprintf (stderr, "r");
      }
    }
    //    ivec_print (cluster_size, k);     

    hashconvergence_new = fvec_sum (centroids_new, k * d);

    tmp = centroids_new;
    centroids_new = centroids;
    centroids = tmp;

  } while (hashconvergence != hashconvergence_new);
  fprintf (stderr, "\n");

  if (clust_assign_out)
    memcpy (clust_assign_out, clust_assign, n * sizeof (int));

  free (clust_assign);
  free (cluster_size);

  return centroids;
}


/****************************************************************************
 * Blocked and multithreading NN versions
 */

typedef struct {
  const pq_t *pq;

  const pqcode_t *codes1;
  const real *v;
  const pqcode_t *codes2;

  int n1, n2, k;
  int nslice;

  int *nns_out;
  float *nns_dis_out;

  void *distfunc;

} task_params_t;

static void compute_nn_block (void *arg, int slice, int ignored);

void pq_nns_thread (const pq_t * pq, const pqcode_t * codes1,
                    const pqcode_t * codes2, int n1, int n2, int k,
                    int nthread, int *nns_out, real * nns_dis_out)
{

  task_params_t params = {
    pq, codes1, NULL, codes2,
    n1, n2, k, nthread, nns_out, nns_dis_out,
    pq_L2sq
  };

  compute_tasks (nthread, nthread, &compute_nn_block, &params);

}

void pq_asym_nns_thread (const pq_t * pq,
                         const real * v, const pqcode_t * codes2,
                         int n1, int n2, int k, int nthread,
                         int *nns_out, real * nns_dis_out)
{

  task_params_t params = {
    pq, NULL, v, codes2,
    n1, n2, k, nthread, nns_out, nns_dis_out,
    pq_asym_L2sq
  };

  compute_tasks (nthread, nthread, &compute_nn_block, &params);

}


#define MIN(a,b) ((a)<(b) ? (a) : (b))

static void compute_nn_block (void *arg, int ignored, int slice)
{
  task_params_t *t = arg;

  int n1 = t->n1, n2 = t->n2, k = t->k;
  const pqcode_t *codes1 = (pqcode_t *) t->codes1;
  const real *v = t->v;
  int *nns_out = t->nns_out;
  real *nns_dis_out = t->nns_dis_out;

  if (t->nslice > 1) {          /* restrict to a query slice */
    int i0 = n1 * slice / t->nslice;
    int i1 = n1 * (slice + 1) / t->nslice;

    //    printf("restrict to %d-%d/%d\n",i0,i1,n1);

    n1 = i1 - i0;
    codes1 = pqcode_offset (t->pq, codes1, i0);
    v += i0 * t->pq->d;
    nns_out += k * i0;
    nns_dis_out += k * i0;
  }

  int i, j;

  int step1 = MIN (n1, 64), step2 = MIN (n2, 1<<20);

  real *dists = fvec_new (step1 * step2);

  /* allocate all heaps at once */
  long oneh = sizeof (maxheap_t) + k * sizeof (heap_entry_t);
  char *minbuf = malloc ((oneh) * step1);

#define MINS(i) ((maxheap_t*)(minbuf + oneh * i))

  int i1, i2;
  

  for (i1 = 0; i1 < n1; i1 += step1) {  /* loop over query blocks */

    int m1 = MIN (step1, n1 - i1);

    /* clear heaps */

    for (i = 0; i < m1; i++) {
      MINS (i)->n = k;
      MINS (i)->i = 0;
    }

    for (i2 = 0; i2 < n2 ; i2 += step2) {        /* loop over base blocks */

      int m2 = MIN (step2, n2 - i2);

      /* compute distances */

      if (t->distfunc == pq_L2sq) {

        for (i = 0; i < m1; i++)
          pq_L2sq (t->pq, pqcode_offset (t->pq, codes1, i + i1),
                   pqcode_offset (t->pq, (pqcode_t *) t->codes2, i2), m2,
                   dists + i * m2);

      } else if (t->distfunc == pq_asym_L2sq) {
        int d = t->pq->d;
        const real *q = v + i1 * d;

        for (i = 0; i < m1; i++) {
          pq_asym_L2sq (t->pq, q,
                        pqcode_offset (t->pq, (pqcode_t *) t->codes2, i2), m2,
                        dists + i * m2);
          q += d;
        }
      } else
        assert(0);

      /* update minima */
      for (i = 0 ; i < m1 ; i++)
        maxheap_add_multiple (MINS (i), i2, m2, dists + i * m2);

    }

    /* collect minima for this query block */
    for (i = 0; i < m1; i++) {
      maxheap_t *mh = MINS (i);
      int *nns = nns_out + (i + i1) * k;
      real *nns_dis = nns_dis_out + (i + i1) * k;

      maxheap_sort (mh);
      
      for (j = 0 ; j < mh->i ; j++) {
        nns[j] = mh->elts[j].label;
        nns_dis[j] = mh->elts[j].val;
      }

      for (; j < k ; j++)
        nns[j] = -1;
    }

  }
  free (minbuf);
  free (dists);

#undef MINS
}
