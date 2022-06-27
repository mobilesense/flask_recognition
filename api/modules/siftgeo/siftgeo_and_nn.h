/*---------------------------------------------------------------------------*/

#ifndef siftgeo_and_nn_INCLUDED
#define siftgeo_and_nn_INCLUDED

/*---------------------------------------------------------------------------*/
/* includes functions that depend both on nn.h and siftgeo.h. They are       */
/* not included in siftgeo.c as this would require to link with Lapack       */
/* & friends every time siftgeo.o is used                                    */
/*---------------------------------------------------------------------------*/

#include "siftgeo.h"

/*---------------------------------------------------------------------------*/
/*! @addtogroup siftgeo
 *  @{
 */

/*! @brief Converts to ffq array (see utils/nn.h) */
float *siftgeo_to_fvecs (const point_t * corners, int nc, int *d_out);

/*! @brief Converts existing ffq array, returns dim */
int pointset_into_fvecs (const pointset_t *ps,float *f);



/* @} */

/*---------------------------------------------------------------------------*/

#endif

/*---------------------------------------------------------------------------*/

