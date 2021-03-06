
#ifndef SORTING_H_INCLUDED
#define SORTING_H_INCLUDED

/* Various sorting functions + a few simple array functions that can
   be called from python efficiently */

/*! 
 tab is a n-element table
 fills maxes such that

 tab[maxes[0]] >= tab[maxes[1]] >= ... >= tab[maxes[k-1]] >= tab[i] 

 for all i not in maxes.
 

*/
void fvec_find_k_max(const float *tab,int n,
		     int *maxes, int k);


void fvec_find_k_min(const float *tab, int n,
		     int *maxes, int k);


/*! finds the ranks of vals[i] for i=0..nval-1 in tab if it was sorted
 * by *decreasing* order
 * minranks[i]-1 is the highest index of values > vals[i]
 * maxranks[i] is the lowest index of values < vals[i]
 * both may be NULL if you're not interested
 */ 
void fvec_ranks_of(const float *tab,int n,
                     const float *vals,int nval,
                     int *minranks,int *maxranks);

/* idem but ranks in increasing array */
void fvec_ranks_inc_of(const float *tab, int n,
                         const float *vals, int nval,
                         int *minranks, int *maxranks);


/*---------------------------------------------------------------------------*/
/* Simple index functions (useful to call from C)                            */
/*---------------------------------------------------------------------------*/

/*! @brief Replace ilabels[i] with the location of ilabels[i] in the table labels.
 *
 *  on input: labels[nres],ilabels[nilabels]\n
 *  on output: labels[ilabels_out[i]]=ilabels[i] for 0<=i<nilabels or -1 if there is none
*/
void find_labels (int *labels, int nres, int *ilabels, int nilabels);

/*! count nb of 0s in array */
int fvec_count_0(const float *val,int n); 

float fvec_min(const float *f, long n);
float fvec_max(const float *f, long n);


/*! computes the median of a float array. Array modified on output! */
float fvec_median (float *f,int n);

/*! computes the arg min of a float array */
int fvec_arg_min (const float *f, int n);



/*! find quantile so that q elements are <= this quantile. On ouput
  the 0..q-1 elements of f are below this quantile */
float fvec_quantile(float *f,int n,int q);


/*! in-place sort */
void ivec_sort(int *tab, int n);

/*! return permutation to sort an array. Is stable. */
void ivec_sort_index (const int *tab, int n, int *perm);

/* fill-in iperm so that iperm[perm[i]]=i for i=0..n-1 */
void ivec_invert_perm(const int *perm, int n, int *iperm); 

/*! return permutation to sort an array. Is stable. */
void fvec_sort_index (const float *tab, int n, int *perm);

/*! sort according to the input permutation. The permutation is 
   typically generated using the ivec_sort_index function. In that 
   case the function perform the sort accordingly. 
*/
void ivec_sort_by_permutation (int * v, const int * order, int n);

/*! count occurrences of val in sorted vector */
int ivec_sorted_count_occurrences(const int *v,int n,int val);

/*! find index of highest value <= val (=-1 if all values are > val) */
int ivec_sorted_find(const int *v,int n,int val);

/*! count unique occurrences  */
int ivec_sorted_count_unique(const int *v,int n);


/*---------------------------------------------------------------------------
 Maxheap functions. (see
 http://en.wikipedia.org/wiki/Binary_heap). When adding elements to
 the maxheap, the elements with smallest val seen so far are in
 
 mh->elts[0] ... mh->elts[mh->i-1] 
  
---------------------------------------------------------------------------*/

typedef struct { /* label & corresponding value */
  float val;
  int label; 
} heap_entry_t;

typedef struct {
  int i; /* number of filled elements */
  int n; /* capacity of the heap */
  heap_entry_t elts[0]; /* to allocate in 1 malloc */
} maxheap_t;

maxheap_t *maxheap_new(int n);
void maxheap_delete(maxheap_t *mh);
int maxheap_add(maxheap_t *mh,int label,float val); 
void maxheap_pop (maxheap_t * mh);

/*! maxheap_add (i+label_base,vals[i]) for i=0..n-1 */
void maxheap_add_multiple(maxheap_t *mh,int label_base,int n,const float *vals); 

void maxheap_add_multiple_labels(maxheap_t *mh,const int *labels,const float *vals,int n); 

/* sort entries in elts by vals in increasing order. WARNING: don't do
   maxheap_add's afterwards */
void maxheap_sort(maxheap_t *mh);

/* sort entries by labels in elts in increasing order. WARNING: don't do
   maxheap_add's afterwards */
void maxheap_sort_labels(maxheap_t *mh);





/* merge k ordered sets defined by 
 * 
 * [(lists[i][j],vals[i][j]),j=0..sizes[i]-1]
 * 
 * for i=0..k-1 
 * 
 * the individual lists are supposes to be ordered already.
 * 
 * returns total number of elements (=sum(sizes[i],i=0..k-1))
 */
int merge_ordered_sets (const int **labels, const float **vals,
                        const int *sizes, int k,
                        int **labels_out, float **vals_out); 


/* finds the smallest value m of vals, compresses array labels by
removing labels[i] for which vals[i] < m * ratio returns new size
of labels array */
int compress_labels_by_disratio(int *labels, const float *vals, int n, float ratio); 


#endif
