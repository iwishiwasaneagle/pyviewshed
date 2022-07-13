
/****************************************************************************
 *
 * MODULE:       r.viewshed
 *
 * AUTHOR(S):    Laura Toma, Bowdoin College - ltoma@bowdoin.edu
 *               Yi Zhuang - yzhuang@bowdoin.edu

 *               Ported to GRASS by William Richard -
 *               wkrichar@bowdoin.edu or willster3021@gmail.com
 *               Markus Metz: surface interpolation
 *
 * Date:         april 2011
 *
 * PURPOSE: To calculate the viewshed (the visible cells in the
 * raster) for the given viewpoint (observer) location.  The
 * visibility model is the following: Two points in the raster are
 * considered visible to each other if the cells where they belong are
 * visible to each other.  Two cells are visible to each other if the
 * line-of-sight that connects their centers does not intersect the
 * terrain. The terrain is NOT viewed as a tesselation of flat cells,
 * i.e. if the line-of-sight does not pass through the cell center,
 * elevation is determined using bilinear interpolation.
 * The viewshed algorithm is efficient both in
 * terms of CPU operations and I/O operations. It has worst-case
 * complexity O(n lg n) in the RAM model and O(sort(n)) in the
 * I/O-model.  For the algorithm and all the other details see the
 * paper: "Computing Visibility on * Terrains in External Memory" by
 * Herman Haverkort, Laura Toma and Yi Zhuang.
 *
 * COPYRIGHT: (C) 2008 by the GRASS Development Team
 *
 * This program is free software under the GNU General Public License
 * (>=v2). Read the file COPYING that comes with GRASS for details.
 *
 *****************************************************************************/

#include "eventlist.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "defs.h"
#include "logging.h"

/* forced to use this because DistanceCompare::compare has troubles if
   i put it inside the class */
Viewpoint globalVP;

/*///////////////////////////////////////////////////////////////////// */
double calculate_angle(double eventX, double eventY, double viewpointX,
                       double viewpointY) {
    double angle = atan(fabs(eventY - viewpointY) / fabs(eventX - viewpointX));

    /*M_PI is defined in math.h to represent 3.14159... */
    if (viewpointY == eventY && eventX > viewpointX) {
        return 0; /*between 1st and 4th quadrant */
    } else if (eventX > viewpointX && eventY < viewpointY) {
        /*first quadrant */
        return angle;

    } else if (viewpointX == eventX && viewpointY > eventY) {
        /*between 1st and 2nd quadrant */
        return (M_PI) / 2;

    } else if (eventX < viewpointX && eventY < viewpointY) {
        /*second quadrant */
        return (M_PI - angle);

    } else if (viewpointY == eventY && eventX < viewpointX) {
        /*between 1st and 3rd quadrant */
        return M_PI;

    } else if (eventY > viewpointY && eventX < viewpointX) {
        /*3rd quadrant */
        return (M_PI + angle);

    } else if (viewpointX == eventX && viewpointY < eventY) {
        /*between 3rd and 4th quadrant */
        return (M_PI * 3.0 / 2.0);
    } else if (eventX > viewpointX && eventY > viewpointY) {
        /*4th quadrant */
        return (M_PI * 2.0 - angle);
    }
    assert(eventX == viewpointX && eventY == viewpointY);
    return 0;
}

/* ------------------------------------------------------------ */
/* calculate the exact position of the given event, and store them in x
   and y.
   quadrants:  1 2
   3 4
   ----->x
   |
   |
   |
   V y
 */
void calculate_event_position(AEvent e, dimensionType viewpointRow,
                              dimensionType viewpointCol, double *y,
                              double *x) {
    assert(x && y);
    *x = 0;
    *y = 0;

    if (e.eventType == CENTER_EVENT) {
        /*FOR CENTER_EVENTS */
        *y = e.row;
        *x = e.col;
        return;
    }

    if (e.row < viewpointRow && e.col < viewpointCol) {
        /*first quadrant */
        if (e.eventType == ENTERING_EVENT) {
            /*if it is ENTERING_EVENT */
            *y = e.row - 0.5;
            *x = e.col + 0.5;
        } else {
            /*otherwise it is EXITING_EVENT */
            *y = e.row + 0.5;
            *x = e.col - 0.5;
        }

    } else if (e.col == viewpointCol && e.row < viewpointRow) {
        /*between the first and second quadrant */
        if (e.eventType == ENTERING_EVENT) {
            /*if it is ENTERING_EVENT */
            *y = e.row + 0.5;
            *x = e.col + 0.5;
        } else {
            /*otherwise it is EXITING_EVENT */
            *y = e.row + 0.5;
            *x = e.col - 0.5;
        }

    } else if (e.col > viewpointCol && e.row < viewpointRow) {
        /*second quadrant */
        if (e.eventType == ENTERING_EVENT) {
            /*if it is ENTERING_EVENT */
            *y = e.row + 0.5;
            *x = e.col + 0.5;
        } else { /*otherwise it is EXITING_EVENT */
            *y = e.row - 0.5;
            *x = e.col - 0.5;
        }

    } else if (e.row == viewpointRow && e.col > viewpointCol) {
        /*between the second and the fourth quadrant */
        if (e.eventType == ENTERING_EVENT) {
            /*if it is ENTERING_EVENT */
            *y = e.row + 0.5;
            *x = e.col - 0.5;
        } else {
            /*otherwise it is EXITING_EVENT */
            *y = e.row - 0.5;
            *x = e.col - 0.5;
        }

    } else if (e.col > viewpointCol && e.row > viewpointRow) {
        /*fourth quadrant */
        if (e.eventType == ENTERING_EVENT) {
            /*if it is ENTERING_EVENT */
            *y = e.row + 0.5;
            *x = e.col - 0.5;
        } else {
            /*otherwise it is EXITING_EVENT */
            *y = e.row - 0.5;
            *x = e.col + 0.5;
        }

    } else if (e.col == viewpointCol && e.row > viewpointRow) {
        /*between the third and fourth quadrant */
        if (e.eventType == ENTERING_EVENT) {
            /*if it is ENTERING_EVENT */
            *y = e.row - 0.5;
            *x = e.col - 0.5;
        } else {
            /*otherwise it is EXITING_EVENT */
            *y = e.row - 0.5;
            *x = e.col + 0.5;
        }

    } else if (e.col < viewpointCol && e.row > viewpointRow) {
        /*third quadrant */
        if (e.eventType == ENTERING_EVENT) {
            /*if it is ENTERING_EVENT */
            *y = e.row - 0.5;
            *x = e.col - 0.5;
        } else {
            /*otherwise it is EXITING_EVENT */
            *y = e.row + 0.5;
            *x = e.col + 0.5;
        }

    } else if (e.row == viewpointRow && e.col < viewpointCol) {
        /*between first and third quadrant */
        if (e.eventType == ENTERING_EVENT) { /*if it is ENTERING_EVENT */
            *y = e.row - 0.5;
            *x = e.col + 0.5;
        } else {
            /*otherwise it is EXITING_EVENT */
            *y = e.row + 0.5;
            *x = e.col + 0.5;
        }
    } else {
        /*must be the viewpoint cell itself */
        assert(e.row == viewpointRow && e.col == viewpointCol);
        *x = e.col;
        *y = e.row;
    }

    assert(fabs(*x - e.col) < 1 && fabs(*y - e.row) < 1);
    /*
    if ((fabs(*x -e.col) >=1) || (fabs(*y -e.row) >=1)) {
       _LOG_WARNING("x-e.col=%f, y-e.row=%f ", fabs(*x -e.col), fabs(*y
    -e.row)); print_event(e, 0); _LOG_WARNING("vp=(%d, %d), x=%.3f, y=%.3f",
    viewpointRow, viewpointCol, *x, *y); exit(1);
       }
    */
    return;
}

void calculate_event_row_col(AEvent e, dimensionType viewpointRow,
                             dimensionType viewpointCol, int *y, int *x) {
    assert(x && y);
    *x = 0;
    *y = 0;

    if (e.eventType == CENTER_EVENT) {
        _ABORT_F(
            "calculate_event_row_col() must not be called for CENTER events");
    }

    if (e.row < viewpointRow && e.col < viewpointCol) {
        /*first quadrant */
        if (e.eventType == ENTERING_EVENT) {
            /*if it is ENTERING_EVENT */
            *y = e.row - 1;
            *x = e.col + 1;
        } else {
            /*otherwise it is EXITING_EVENT */
            *y = e.row + 1;
            *x = e.col - 1;
        }

    } else if (e.col == viewpointCol && e.row < viewpointRow) {
        /*between the first and second quadrant */
        if (e.eventType == ENTERING_EVENT) {
            /*if it is ENTERING_EVENT */
            *y = e.row + 1;
            *x = e.col + 1;
        } else {
            /*otherwise it is EXITING_EVENT */
            *y = e.row + 1;
            *x = e.col - 1;
        }

    } else if (e.col > viewpointCol && e.row < viewpointRow) {
        /*second quadrant */
        if (e.eventType == ENTERING_EVENT) {
            /*if it is ENTERING_EVENT */
            *y = e.row + 1;
            *x = e.col + 1;
        } else { /*otherwise it is EXITING_EVENT */
            *y = e.row - 1;
            *x = e.col - 1;
        }

    } else if (e.row == viewpointRow && e.col > viewpointCol) {
        /*between the second and the fourth quadrant */
        if (e.eventType == ENTERING_EVENT) {
            /*if it is ENTERING_EVENT */
            *y = e.row + 1;
            *x = e.col - 1;
        } else {
            /*otherwise it is EXITING_EVENT */
            *y = e.row - 1;
            *x = e.col - 1;
        }

    } else if (e.col > viewpointCol && e.row > viewpointRow) {
        /*fourth quadrant */
        if (e.eventType == ENTERING_EVENT) {
            /*if it is ENTERING_EVENT */
            *y = e.row + 1;
            *x = e.col - 1;
        } else {
            /*otherwise it is EXITING_EVENT */
            *y = e.row - 1;
            *x = e.col + 1;
        }

    } else if (e.col == viewpointCol && e.row > viewpointRow) {
        /*between the third and fourth quadrant */
        if (e.eventType == ENTERING_EVENT) {
            /*if it is ENTERING_EVENT */
            *y = e.row - 1;
            *x = e.col - 1;
        } else {
            /*otherwise it is EXITING_EVENT */
            *y = e.row - 1;
            *x = e.col + 1;
        }

    } else if (e.col < viewpointCol && e.row > viewpointRow) {
        /*third quadrant */
        if (e.eventType == ENTERING_EVENT) {
            /*if it is ENTERING_EVENT */
            *y = e.row - 1;
            *x = e.col - 1;
        } else {
            /*otherwise it is EXITING_EVENT */
            *y = e.row + 1;
            *x = e.col + 1;
        }

    } else if (e.row == viewpointRow && e.col < viewpointCol) {
        /*between first and third quadrant */
        if (e.eventType == ENTERING_EVENT) { /*if it is ENTERING_EVENT */
            *y = e.row - 1;
            *x = e.col + 1;
        } else {
            /*otherwise it is EXITING_EVENT */
            *y = e.row + 1;
            *x = e.col + 1;
        }
    } else {
        /*must be the viewpoint cell itself */
        _LOG_F(1, "calculate_event_row_col() called for viewpoint cell itself");
        assert(e.row == viewpointRow && e.col == viewpointCol);
        *x = e.col;
        *y = e.row;
    }

    /* assert(fabs(*x - e.col) <= 1 && fabs(*y - e.row) <= 1); */

    if ((abs(*x - e.col) > 1) || (abs(*y - e.row) > 1)) {
        _LOG_WARNING("calculate_event_row_col() :");
        _LOG_WARNING("x-e.col=%d, y-e.row=%d", abs(*x - e.col),
                     abs(*y - e.row));
        print_event(e);
        _LOG_WARNING("vp=(%d, %d), x=%d, y=%d", viewpointRow, viewpointCol, *x,
                     *y);
        exit(1);
    }

    return;
}

/* ------------------------------------------------------------ */
void print_event(AEvent a) {
    char c = '0';

    if (a.eventType == ENTERING_EVENT) c = 'E';
    if (a.eventType == EXITING_EVENT) c = 'X';
    if (a.eventType == CENTER_EVENT) c = 'Q';

    _LOG_VERBOSE("ev=[(%3d, %3d), e=%8.1f a=%4.2f t=%c] ", a.row, a.col,
                 a.elev[1], a.angle, c);
    return;
}

/* ------------------------------------------------------------ */
/*computes the distance from the event to the viewpoint. Note: all 3
   //events associate to a cell are considered at the same distance, from
   //the center of the cell to the viewpoint */
double get_square_distance_from_viewpoint(const AEvent &a,
                                          const Viewpoint &vp) {
    double eventy, eventx, dist;

    calculate_event_position(a, vp.row, vp.col, &eventy, &eventx);

#if PROJECTION == PROJECTION_LL
    struct Cell_head window;

    Rast_get_window(&window);

    dist = G_distance(Rast_col_to_easting(vp.col + 0.5, &window),
                      Rast_row_to_northing(vp.row + 0.5, &window),
                      Rast_col_to_easting(eventx + 0.5, &window),
                      Rast_row_to_northing(eventy + 0.5, &window));

    dist = dist * dist;
#else
    /*don't take sqrt, it is expensive; suffices for comparison */
    dist = (eventx - vp.col) * (eventx - vp.col) +
           (eventy - vp.row) * (eventy - vp.row);
#endif
    return dist;
}

/* ------------------------------------------------------------ */
/*determines if the point at row,col is outside the maximum distance
   limit.  Return 1 if the point is outside limit, 0 if point is inside
   limit. */
int is_point_outside_max_dist(Viewpoint vp, GridHeader hd, dimensionType row,
                              dimensionType col, float maxDist) {
    /* it is not too smart to compare floats */
    if ((int)maxDist == INFINITY_DISTANCE) return 0;

    float dist;
#if PROJECTION == PROJECTION_LL
    struct Cell_head window;

    Rast_get_window(&window);

    dist = G_distance(Rast_col_to_easting(vp.col + 0.5, &window),
                      Rast_row_to_northing(vp.row + 0.5, &window),
                      Rast_col_to_easting(col + 0.5, &window),
                      Rast_row_to_northing(row + 0.5, &window));

    dist = dist * dist;
#else
    /*don't take sqrt, it is expensive; suffices for comparison */
    dist = (col - vp.col) * (col - vp.col) + (row - vp.row) * (row - vp.row);
#endif

    /* multiplication is faster than sqrt */
    if ((maxDist * maxDist) < dist) {
        return 1;
    }

    return 0;
}

/* ------------------------------------------------------------
   //note: this is expensive because distance is not stored in the event
   //and must be computed on the fly */
int DistanceCompare::compare(const AEvent &a, const AEvent &b) {
    /*calculate distance from viewpoint
       //don't take sqrt, it is expensive; suffices for comparison */
    double da, db;

    da = get_square_distance_from_viewpoint(a, globalVP);
    db = get_square_distance_from_viewpoint(b, globalVP);

    if (da > db) {
        return 1;
    } else if (da < db) {
        return -1;
    } else {
        return 0;
    }
    return 0;
}

/* ------------------------------------------------------------ */
int RadialCompare::compare(const AEvent &a, const AEvent &b) {
    if (a.row == b.row && a.col == b.col && a.eventType == b.eventType)
        return 0;

    assert(a.angle >= 0 && b.angle >= 0);

    if (a.angle > b.angle) {
        return 1;
    } else if (a.angle < b.angle) {
        return -1;
    } else {
        /*a.angle == b.angle */
        if (a.eventType == EXITING_EVENT) return -1;
        if (b.eventType == EXITING_EVENT) return 1;
        if (a.eventType == ENTERING_EVENT) return 1;
        if (b.eventType == ENTERING_EVENT) return -1;
        return 0;
    }
}

/* ------------------------------------------------------------ */
/* a copy of the function above is needed by qsort, when the
   computation runs in memory */

int radial_compare_events(const void *x, const void *y) {
    AEvent *a, *b;

    a = (AEvent *)x;
    b = (AEvent *)y;
    if (a->row == b->row && a->col == b->col && a->eventType == b->eventType)
        return 0;

    assert(a->angle >= 0 && b->angle >= 0);

    if (a->angle > b->angle) {
        return 1;
    } else if (a->angle < b->angle) {
        return -1;
    } else {
        /*a->angle == b->angle */
        if (a->eventType == EXITING_EVENT)
            return -1;
        else if (a->eventType == ENTERING_EVENT)
            return 1;
        return 0;
    }
}
