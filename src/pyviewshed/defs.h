//
// Created by jhewers on 11/07/22.
//

#ifndef PYVIEWSHED_DEFS_H
#define PYVIEWSHED_DEFS_H

#define _(str) (str)

#define _FREE free
#define _MALLOC malloc
//
// typedef struct _viewoptions {
//    float obsElev;
//    float tgtElev;
//    float maxDist;
//    float horizontal_angle_min;
//    float horizontal_angle_max;
//} _ViewOptions;

struct Cell_head {
    Cell_head(){};
    Cell_head(int rows, int cols, double ew_res, double ns_res, double north,
              double east, double south, double west)
        : rows(rows),
          cols(cols),
          ew_res(ew_res),
          ns_res(ns_res),
          north(north),
          east(east),
          south(south),
          west(west){};

    int format;

    int rows;

    int cols;

    int proj;

    double ew_res;

    double ns_res;

    double north;

    double south;

    double east;

    double west;
};

/* From Grass docs:
 *
 * - PROJECTION_XY      0 - x,y (Raw imagery)
 * - PROJECTION_UTM     1 - UTM   Universal Transverse Mercator
 * - PROJECTION_SP      2 - State Plane (in feet) - not used, removed
 * - PROJECTION_LL      3 - Latitude-Longitude
 * - PROJECTION_OTHER  99 - others
 */
#define PROJECTION_XY 0
#define PROJECTION_LL 3

#define PROJECTION PROJECTION_XY

#define NOVALUE -1

#endif  // PYVIEWSHED_DEFS_H
