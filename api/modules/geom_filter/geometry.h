#ifndef GEOMETRY_INCLUDED
#define GEOMETRY_INCLUDED


#include <siftgeo/siftgeo.h>




/*-----------------------------------------------------------------
 * Geometry related functions
 *-----------------------------------------------------------------*/

struct pointmatch_t;

/* return <0 for bad estimation */
int estimate_affine_transform(point_t *qpts,point_t *dbpts,
                              struct pointmatch_t **pms,int npm,
                              double aff[6]);


int find_agreeing_matches(point_t *qpts,point_t *dbpts,
                          double aff[6],
                          struct pointmatch_t *pm,
                          double thresh,
                          struct pointmatch_t **agreeing_matches, int* ids);

/* Count distinct query or dbpoints. When count reaches stopat, 
 * stop counting
 * algorithm is crappy as result depends on order of matches.
 */

int count_distinct_points(point_t *pts,
                          struct pointmatch_t **pms,int npm,
                          int use_query,int stopat,
                          double radius);


#endif
