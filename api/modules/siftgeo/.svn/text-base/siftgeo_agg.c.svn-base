#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include <yael/vector.h>
#include <yael/sorting.h>

#include "siftgeo_agg.h"



/*---------------------------------------------------------------------------*/
/* Elementary operations on atoms                                            */
/*---------------------------------------------------------------------------*/

atoms_t *atoms_new (int n, int d)
{
  atoms_t *atoms = (atoms_t *) malloc (sizeof (*atoms));
  atoms->n = n;
  atoms->d = d;
  atoms->nz = ivec_new (n);
  atoms->idx = (int **) malloc (n * sizeof (*atoms->idx));
  atoms->val = (float **) malloc (n * sizeof (*atoms->val));

  memset (atoms->nz, n * sizeof (*atoms->nz), 0);
  memset (atoms->idx, n * sizeof (*atoms->idx), 0);
  memset (atoms->val, n * sizeof (*atoms->val), 0);

  return atoms;
}


void atoms_free (atoms_t * atoms)
{
  int i;
  for (i = 0; i < atoms->n; i++) {
    free (atoms->idx[i]);
    free (atoms->val[i]);
  }

  free (atoms->idx);
  free (atoms->val);
  free (atoms->nz);
  free (atoms);
}


void atoms_add_spvec (atoms_t * atoms, int no, int *idx, float *val, int n)
{
  /*  printf ("%d %d %d\n", no, atoms->n, n); */
  assert (no >= 0 && no < atoms->n);
  atoms->nz[no] = n;
  atoms->idx[no] = malloc (n * sizeof (**atoms->idx));
  atoms->val[no] = malloc (n * sizeof (**atoms->val));

  memcpy (atoms->idx[no], idx, n * sizeof (**atoms->idx));
  if (val)
    memcpy (atoms->val[no], val, n * sizeof (**atoms->val));
  else {
    int i;
    for (i = 0; i < n; i++)
      atoms->val[no][i] = 1.0;
  }
}


#define ATOMS_READ_ERROR(ret, expected_ret)                             \
  {									\
    if (ret != (expected_ret)) {					\
      fprintf (stderr, "# Unable to read the atoms from file %s\n",     \
	       filename);						\
      return NULL;							\
    }									\
  }

#define ATOMS_WRITE_ERROR(ret, expected_ret)                            \
  {									\
    if (ret != (expected_ret)) {					\
      fprintf (stderr, "# Unable to write the atoms into file %s\n",    \
	       filename);						\
      return;	                                                        \
    }									\
  }


atoms_t *atoms_read (const char *filename)
{
  int ret, i, n, d;
  atoms_t *atoms;

  FILE *f = fopen (filename, "r");
  if (!f) {
    fprintf (stderr, "# Unable to open atoms file %s for reading\n",
             filename);
    return NULL;
  }

  /* read the number of atoms and their dimension */
  ret = fread (&n, sizeof (n), 1, f);
  ATOMS_READ_ERROR (ret, 1);
  ret = fread (&d, sizeof (d), 1, f);
  ATOMS_READ_ERROR (ret, 1);
  assert (n >= 1 && d >= 1);

  /* Create the structure and read the number of non-zeros values per atom */
  atoms = atoms_new (n, d);
  ret = fread (atoms->nz, sizeof (*atoms->nz), n, f);
  ATOMS_READ_ERROR (ret, n);

  /* Allocate and read the atoms */
  for (i = 0; i < atoms->n; i++) {
    atoms->idx[i] = malloc (atoms->nz[i] * sizeof (**atoms->idx));
    assert (atoms->idx[i]);

    atoms->val[i] = malloc (atoms->nz[i] * sizeof (**atoms->val));
    assert (atoms->val[i]);

    ret = fread (atoms->idx[i], sizeof (**atoms->idx), atoms->nz[i], f);
    ATOMS_READ_ERROR (ret, atoms->nz[i]);

    ret = fread (atoms->val[i], sizeof (**atoms->val), atoms->nz[i], f);
    ATOMS_READ_ERROR (ret, atoms->nz[i]);
  }

  fclose (f);
  return atoms;
}


void atoms_write (const char *filename, const atoms_t * atoms)
{
  int i, ret;

  FILE *f = fopen (filename, "w");
  if (!f) {
    fprintf (stderr, "# Unable to create the atom file %s\n", filename);
    return;
  }

  ret = fwrite (&atoms->n, sizeof (atoms->n), 1, f);
  ATOMS_WRITE_ERROR (ret, 1);
  ret = fwrite (&atoms->d, sizeof (atoms->d), 1, f);
  ATOMS_WRITE_ERROR (ret, 1);
  assert (atoms->n >= 1 && atoms->d >= 1);

  ret = fwrite (atoms->nz, sizeof (*atoms->nz), atoms->n, f);
  ATOMS_WRITE_ERROR (ret, atoms->n);

  /* Allocate and read the atoms */
  for (i = 0; i < atoms->n; i++) {
    ret = fwrite (atoms->idx[i], sizeof (**atoms->idx), atoms->nz[i], f);
    ATOMS_WRITE_ERROR (ret, atoms->nz[i]);

    ret = fwrite (atoms->val[i], sizeof (**atoms->val), atoms->nz[i], f);
    ATOMS_WRITE_ERROR (ret, atoms->nz[i]);
  }
  fclose (f);
}


void atoms_display (const atoms_t * atoms, int verboselevel)
{
  int i;
  printf ("n = %d, d = %d\n", atoms->n, atoms->d);

  if (verboselevel >= 1) {
    printf ("nz = ");
    ivec_print (atoms->nz, atoms->n);
  }

  if (verboselevel >= 2) {
    for (i = 0; i < atoms->n; i++) {
      printf ("atoms[%d] (%d nz) = ", i, atoms->nz[i]);
      ivec_print (atoms->idx[i], atoms->nz[i]);
      fvec_print (atoms->val[i], atoms->nz[i]);
    }
  }
}


float *atoms_mul_fvec (const atoms_t * atoms, float * v)
{
  int i, j;
  int n = atoms->n;

  float * res = fvec_new_0 (n);

  for (i = 0 ; i < n ; i++)
    for (j = 0 ; j < atoms->nz[i] ; j++) 
      res[i] += atoms->val[i][j] * v[atoms->idx[i][j]];

  return res;
}


void atoms_sort_indices(atoms_t * atoms) 
{
  int i,j;
  for(i=0;i<atoms->n;i++) {
    int *idx=atoms->idx[i];
    float *val=atoms->val[i];
    int nz=atoms->nz[i];

    int *perm=ivec_new(nz);
    ivec_sort_index(idx,nz,perm);

    atoms->idx[i]=ivec_new(nz);
    atoms->val[i]=fvec_new(nz);    
    for(j=0;j<nz;j++) {
      atoms->idx[i][j]=idx[perm[j]];
      atoms->val[i][j]=val[perm[j]];
    }

    free(idx);
    free(val);
    free(perm);
  }

}

void atoms_crop(atoms_t *atoms,int n0,int n1) {
  assert(0<=n0 && n0<=n1 && n1<=atoms->n);
  int i;
  for(i=0;i<atoms->n;i++) if(!(n0<=i || i<n1)) {
    free(atoms->idx[i]);
    free(atoms->val[i]);
  }
      
  memmove(atoms->idx,atoms->idx+n0,sizeof(*atoms->idx)*(n1-n0));
  memmove(atoms->nz, atoms->nz+n0, sizeof(*atoms->nz )*(n1-n0));
  memmove(atoms->val,atoms->val+n0,sizeof(*atoms->val)*(n1-n0));

  atoms->n=n1-n0;
  
}


/*---------------------------------------------------------------------------*/
/* Miscellaneous                                                             */
/*---------------------------------------------------------------------------*/

static double sqrdiff (int x, int y)
{
  double diff = x - y;
  return diff * diff;
}


/* Compute the distances between a descriptor and a set of descriptors */
void point_to_points_dis (const point_t * set, const point_t * des, int n,
                          float *disout)
{
  int i, j;
  int d = set[0].dim;

  for (i = 0; i < n; i++) {
    double dissqr = 0;
    for (j = 0; j < d; j++)
      dissqr += sqrdiff (des->desc[j], set[i].desc[j]);
    disout[i] = dissqr;
  }
}


#include <math.h>
/* return the weighting of the scale */
static float scalew (const point_t * pt)
{
  return 1;
  /*  return log (10 + pt->geom.scale); */
}


/* Convert a Bag-of-features into a packed representation */
int vwsgeo_to_bof (vwgeoset_t * vgs, int **idx_out, float **v_out)
{
  if (vgs->n == 0) {
    *idx_out = NULL;
    *v_out = NULL;
    return 0;
  }

  int i, ii = 0, nz = 1;

  point_t *desc = vgs->pts;
  vwgeoset_sort (vgs);

  /* first count the number of non-zeros positions */
  for (i = 1; i < vgs->n; i++)
    if (desc[i].vw != desc[i - 1].vw)
      nz++;

  /* Allocate the bof */
  int *idx = ivec_new (nz);
  float *v = fvec_new (nz);

  idx[0] = desc[0].vw;
  v[0] = scalew(desc + 0);

  for (i = 1; i < vgs->n; i++) {
    assert (ii < nz);

    if (desc[i].vw > idx[ii]) {
      ii++;
      idx[ii] = desc[i].vw;
      v[ii] = scalew(desc + i);
    } else
      v[ii] += scalew(desc + i);
  }

  *idx_out = idx;
  *v_out = v;
  return nz;
}

/*---------------------------------------------------------------------------*/
/* Atoms based on Hadamard basis vectors intersection                        */
/*---------------------------------------------------------------------------*/


/* the Hadamard base intersection does not change if h[i] is replaced
 * with h[i]^h[j] for i!=j.  The function initializes masks and
 * reduces h to be in the following format (x stands for any bit, y
 * for uninitialized):
 *
 *           h              masks
 *       lsb     msb     lsb     msb
 * 0     xxxxxxxx1       111111111
 * .     xxxxxx100       111111100
 * .     xxxxx1000       111111000
 * l2-1  x10000000       110000000
 * .     000000000       yyyyyyyyy
 * l-1   000000000       yyyyyyyyy
 * 
 * The function returns l2.
 */
static int canonicalize_h (int l, int *h, int *masks)
{
  int i, j;

  /* normalize h */
  for (i = 0; i < l; i++) {

    /* find maximum */
    int m = 0, m_index = -1;
    for (j = i; j < l; j++)
      if (h[j] > m) {
        m_index = j;
        m = h[j];
      }

    {                           /* swap to index i */
      h[m_index] = h[i];
      h[i] = m;
    }

    if (m == 0)
      break;

    /* which bit does this maximum cancel? */
    int bit_to_cancel = 1;
    while (m >= bit_to_cancel)
      bit_to_cancel <<= 1;
    bit_to_cancel = bit_to_cancel >> 1;

    masks[i] = bit_to_cancel | (bit_to_cancel - 1);

    /* cancel that bit on other h's */
    for (j = i + 1; j < l; j++)
      if (h[j] >= bit_to_cancel)
        h[j] ^= m;

  }

  return i;
}

static int count_bits (int x)
{
  int n = 0;
  while (x) {
    x &= x - 1;
    n++;
  }
  return n;
}


/* Produce combinations of bits 
 * 
 *    mask is 1..10..0b with n0 0's and n1 1's 
 * 
 * return all values v that verify 
 * 
 *    nbit(v & ci) is even 
 *    (v & ~mask) == val
 * 
 * tab is filled with v's, return nb of filled-in values
 */
static int bit_combinations (int mask, int val, int ci, int *tab)
{

  int n0 = 0, n1;
  while (!((mask >> n0) & 1))
    n0++;
  n1 = n0;
  while ((mask >> n1) & 1)
    n1++;
  n1 -= n0;

  assert (mask == ((1 << n1) - 1) << n0);

  int nt = 0, i;
  for (i = 0; i < (1 << n1); i++) {
    int v = val | i << n0;
    if ((count_bits (v & ci) & 1) == 0)
      tab[nt++] = v;
  }

  return nt;
}


/* Produce the intersection of hadamard vectors h[0]..h[l-1]. 
 * 
 * On output the h[i]'s are in canonical format and the number l2 of
 * non-0 l2 h entries is returned. The size of sup_indices is d>>l2.
 */
static int hadamard_intersection (int d, int l, int *h, int *sup_indices)
{

  /* printf("h="); ivec_print(h,l); */

  int *masks = ivec_new (l);
  int l2 = canonicalize_h (l, h, masks);

  /* printf("h="); ivec_print(h,l); */

  int tot_nsup = d >> l2;
  int *buf = ivec_new (2 * tot_nsup), *sup = buf, *sup1 = buf + tot_nsup;
  int nsup = 1;
  sup[0] = 0;

  int prev_m = 0;
  int i, j;

  /* loop over bits from low to high and generate combinations that
     have an even nb of 1's when masked with h[i] */

  for (i = l2 - 1; i >= 0; i--) {
    /* generate combinations of bits in masks[i]^prev_m */
    int nsup1 = 0;
    for (j = 0; j < nsup; j++) {
      nsup1 +=
          bit_combinations (masks[i] ^ prev_m, sup[j], h[i], sup1 + nsup1);
      assert (nsup1 <= tot_nsup);
    }
    prev_m = masks[i];
    nsup = nsup1;
    {                           /* swap */
      int *tmp = sup;
      sup = sup1;
      sup1 = tmp;
    }
  }

  if (prev_m != d - 1) {
    /* generate last upper bits */
    int nsup1 = 0;
    for (j = 0; j < nsup; j++) {
      nsup1 += bit_combinations ((d - 1) ^ prev_m, sup[j], 0, sup1 + nsup1);
      assert (nsup1 <= tot_nsup);
    }
    nsup = nsup1;
    {                           /* swap */
      int *tmp = sup;
      sup = sup1;
      sup1 = tmp;
    }
  }
  assert (nsup == tot_nsup);

  /* copy support */
  memcpy (sup_indices, sup, nsup * sizeof (int));
  free (buf);
  free (masks);
  return l2;
}


atoms_t *atoms_new_hadamard_inter (int d, int l, int *h_indices_out)
{
  assert ((d & (d - 1)) == 0 || !"d must be power of 2");

  int *h_indices = h_indices_out ? h_indices_out : ivec_new (d * l);

  int *baseperm = ivec_new_random_perm (d - 1);

  atoms_t *atoms = atoms_new (d - 1, d);

  int i, j;
  int *sup = ivec_new (d);
  for (i = 0; i < d - 1; i++) {
    int *h = h_indices + i * l;

    /* make initial h intersection */
    for (j = 0; j < l; j++)
      h[j] = 1 + baseperm[(i + ((1 << j) >> 1)) % (d - 1)];

    int l2 = hadamard_intersection (d, l, h, sup);

    if (l2 < l) {
      printf
          ("atoms_new_hadamard_inter vector %d: actual intersection %d < l=%d\n",
           i, l2, l);
    }

    atoms_add_spvec (atoms, i, sup, NULL, d >> l2);

  }
  free (sup);
  free (baseperm);
  if (!h_indices_out)
    free (h_indices);

  return atoms;
}


static int ivec_sorted_contains(int n,const int *tab,int v) {
  if(v<tab[0]) return 0;
  int i0=0,i1=n;
  /* tab[i0]<=v<tab[i1] */
  while(i0+1<i1) {
    int med=(i0+i1)/2;
    if(tab[med]<=v) i0=med;
    else            i1=med;
  }
  return tab[i0]==v;
}

/* perms is a (n+1)-by-d table 
 *
 *    |          i0      imax        d-1
 * ---+---------------------------------
 * 0  |  0   1   2   3   4   5   6   7
 * .  |  4   7   0   1   5   3   2   6
 * n  |  6   5   2   4   7   0   1   3
 * 
 * each line is a permutation. The function tries to swap elements in
 * (n,i0:imax-1) so that two elements do not appear together in
 * distinct columns. It returns a boolean that indicates whether this
 * succeded.
 *
 * iperms is a n-by-d table where line i is the inverse permutation of
 * line i of the perms table.
 */
static int sample_assignement(int d,int n,int *perms,const int *iperms,int i0,int imax,
                              int *iter_remain) {

  if(*iter_remain) (*iter_remain)--;
  else return 0;

  if(i0==imax) return 1;

  int *perm=perms+n*d;
  int j,k,ip;
  
  /* conflicting elements. Should never be paired with perm[i0] */
  int *conflicts=ivec_new(n);
  for(j=0;j<n;j++) 
    conflicts[j]=perms[j*d+i0];
  ivec_sort(conflicts,n);

  for(ip=i0;ip<d;ip++) {
    int i=perm[ip];
    
    if(ivec_sorted_contains(n,conflicts,i)) goto next_ip;
      
    for(j=0;j<n;j++) {
      int i1=iperms[i+j*d];
      if(i1>=imax) /* don't care */
        continue;

      assert(i1!=i0);
      
      int n1=i1<i0 ? n+1 : n;

      for(k=0;k<n1;k++) {
        if(k!=j && 
           ivec_sorted_contains(n,conflicts,perms[k*d+i1]))
          goto next_ip;
      }      
    }
    
    /* we found one: switch it to i0 */
    perm[ip]=perm[i0]; 
    perm[i0]=i;
    
    int found=sample_assignement(d,n,perms,iperms,i0+1,imax,iter_remain);
    
    if(found || ! *iter_remain) {
      free(conflicts);
      return found;
    }

    /* failed, switch back */     
    perm[i0]=perm[ip];
    perm[ip]=i;

  next_ip:;

  }

  free(conflicts);

  return 0;
}


atoms_t *atoms_new_random_by_support (int n, int d, int nsup)
{

  /* assert(nsup*nsup <= d || !"will not work with such big support"); */
  assert(n <= d);

  int *perms=ivec_new(nsup*d);
  int *iperms=ivec_new(nsup*d);
  
  int i,j;

  for(i=0;i<d;i++) 
    perms[i]=iperms[i]=i;

  /*
  for(j=0;j<d;j++) 
    printf("%2d ",perms[j]);
  printf("\n"); 
  */

  int max_iter=20*1000;
  int backtrack_remain=8;

  for(i=1;i<nsup;i++) {

    int *perm=perms+i*d;
    {
      int *init_perm=ivec_new_random_perm(d);
      memcpy(perm,init_perm,sizeof(int)*d);
      free(init_perm);
    }

    int iter_remain=max_iter;

    int sa=sample_assignement(d,i,perms,iperms,0,n,&iter_remain);
    
    if(!sa) {
      backtrack_remain--;
      fprintf(stderr,"atoms_new_random_by_support: support %d/%d, "
             "iterations %d/%d backtrack_remain %d.\n",
             i,nsup,max_iter-iter_remain,max_iter,backtrack_remain);
      if(backtrack_remain==0) {
        fprintf(stderr,"atoms_new_random_by_support: setting nsup to %d\n",i);
        nsup=i;
        break;
      }
      i=i-1;
      continue;
    }

    printf("support %d ok in %d iterations\n",i,max_iter-iter_remain);
    
    int *iperm=iperms+i*d;

    for(j=0;j<d;j++) iperm[perm[j]]=j;
/*    
    for(j=0;j<n;j++) 
      printf("%2d ",perm[j]);
    printf("\n");
*/
  }

  atoms_t *atoms = atoms_new (n, d);

  int *sup=ivec_new(nsup);

  for(i=0;i<n;i++) {
    for(j=0;j<nsup;j++) 
      sup[j]=perms[j*d+i];

    atoms_add_spvec(atoms,i,sup,NULL,nsup);    
  }

  free(sup);

  free(perms);
  free(iperms);

  
  return atoms;
}


atoms_t *atoms_new_random_by_sparsity_rate (int n, int d, float rate)
{
  int i, j, nz;
  atoms_t *atoms = atoms_new (n, d);
  int *support = ivec_new (d);

  for (i = 0; i < n; i++) {
    nz = 0;

    for (j = 0 ; j < d ; j++) 
      if (drand48 () < rate) {
	support[nz] = j;
	nz++;
      }

    atoms->nz[i] = nz;
    atoms->idx[i] = ivec_new_copy (support, nz);
    atoms->val[i] = fvec_new_set (nz, 1);
  }

  free (support);
  return atoms;
}


atoms_t * atoms_new_aggregate_components (int n, int d, int nz)
{
  int i, j; 
  int m = (nz * n) / d;
  int v = n / m;        /* v = number of elements per aggregated vocabulary */

  atoms_t *atoms = atoms_new (n, d);
  ivec_set (atoms->nz, n, nz);

  int * order = ivec_new_range (0, d);

  for (j = 0 ; j < m ; j++) {
    for (i = 0 ; i < v ; i++) {
      atoms->idx[j * v + i] = ivec_new_copy (order + i * nz, nz);
      atoms->val[j * v + i] = fvec_new_set (nz, 1);
    }
    ivec_shuffle (order, d);
  }

  return atoms;
}


/*---------------------------------------------------------------------------*/
/* Sparse Hamming Embedding                                                  */
/*---------------------------------------------------------------------------*/

float * atoms_signature_fvec (const atoms_t * atoms, const float*bof, const float * tfidf) {

  /* Optionally, apply the tf-idf scheme */
  if (tfidf)
    fvec_mul (bof, tfidf, atoms->d);

  /* Normalize the bof */
  fvec_normalize (bof, atoms->d, 2);

  /* projection onto the atoms */
  float * v = atoms_mul_fvec (atoms, bof);

  return v;
}



float * atoms_signature (const atoms_t * atoms, vwgeoset_t * vgs, const float * tfidf)
{
  int * idx;
  float * val;
  int nz = vwsgeo_to_bof (vgs, &idx, &val);

  /* Construct the bof */
  float * bof = spfvec_to_fvec (idx, val, nz, atoms->d);

  free (idx);
  free (val);

  float *ret=atoms_signature_fvec(atoms,bof,tfidf);

  free (bof);

  return ret;
}


int * hesp_signature (const atoms_t * atoms, vwgeoset_t * vgs, const float * tfidf)
{
  int i;

  /* projection onto the atoms */
  float * v = atoms_signature (atoms, vgs, tfidf);

  /* find the atoms that have the higher value */
  int * ranks = ivec_new (atoms->n);
  fvec_find_k_max (v, atoms->n, ranks, atoms->n / 2);

  /* convert into a binary sequence */
  int * ret = ivec_new_0 (atoms->n);

  for (i = 0 ; i < atoms->n / 2 ; i++)
    ret[ranks[i]] = 1;

  free (v);
  free (ranks);
  return ret;
}


