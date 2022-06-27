/*---------------------------------------------------------------------------*/

#ifndef fast_alloc_h_INCLUDED
#define fast_alloc_h_INCLUDED

/*---------------------------------------------------------------------------*/
/* Fonctions for fast allocation/deallocation of "units" of fixed            */
/* size. Accesses are not thread-safe                                        */
/*---------------------------------------------------------------------------*/
/*! @addtogroup utils
 *  @{
 */
/*---------------------------------------------------------------------------*/

/*! @brief Pool of allocated units, grouped in blocks */
typedef struct fast_alloc_t fast_alloc_t;

/*---------------------------------------------------------------------------*/

/*! @brief Make a new pool. 
 *
 * us: unit size (must be >= sizeof(void*))\n
 * upb: units per block
 */
fast_alloc_t *fa_alloc(int us,int upb);

void *fa_allocunit(fast_alloc_t*);

void fa_freeunit(fast_alloc_t*,void*unit);

/*! @brief Free units without returning memory to the system */
void fa_clear(fast_alloc_t*);

void fa_free(fast_alloc_t*);

/*! @brief debugging function */
void fa_describe(fast_alloc_t*);

/*! @brief debugging function */
void fa_stats(fast_alloc_t*,int *n_block);

/*---------------------------------------------------------------------------*/
/*! @} */
/*---------------------------------------------------------------------------*/

#endif

/*---------------------------------------------------------------------------*/

