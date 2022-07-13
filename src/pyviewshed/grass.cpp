
/****************************************************************************
 *
 * MODULE:       r.viewshed
 *
 * AUTHOR(S):    Laura Toma, Bowdoin College - ltoma@bowdoin.edu
 *               Yi Zhuang - yzhuang@bowdoin.edu
 *
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

#include "grass.h"

#include <stdio.h>
#include <stdlib.h>

#include <cassert>

#include "defs.h"
#include "logging.h"
#include "visibility.h"

/* calculate ENTER and EXIT event elevation (bilinear interpolation) */
surface_type calculate_event_elevation(AEvent e, int nrows, int ncols,
                                       dimensionType vprow, dimensionType vpcol,
                                       Eigen::Matrix<float, 3, -1> raster) {
    int row1, col1;
    surface_type event_elev;
    G_SURFACE_T elev1, elev2, elev3, elev4;

    calculate_event_row_col(e, vprow, vpcol, &row1, &col1);
    if (row1 >= 0 && row1 < nrows && col1 >= 0 && col1 < ncols) {
        elev1 = raster.row(row1 - e.row + 1)[col1];
        elev2 = raster.row(row1 - e.row + 1)[e.col];
        elev3 = raster.row(1)[col1];
        elev4 = raster.row(1)[e.col];
        if (elev1 == -1 || elev2 == -1 || elev3 == -1 || elev4 == -1)
            event_elev = raster.row(1)[e.col];
        else {
            event_elev = (elev1 + elev2 + elev3 + elev4) / 4.;
        }
    } else
        event_elev = raster.row(1)[e.col];

    return event_elev;
}

/*  ************************************************************ */
/* input: an array capable to hold the max number of events, a raster
   name, a viewpoint and the viewOptions; action: figure out all events
   in the input file, and write them to the event list. data is
   allocated and initialized with all the cells on the same row as the
   viewpoint. it returns the number of events. initialize and fill
   AEvent* with all the events for the map.  Used when solving in
   memory, so the AEvent* should fit in memory.  */
size_t init_event_list_in_memory(AEvent *eventList, Eigen::MatrixXf raster,
                                 Viewpoint *vp, GridHeader *hd,
                                 ViewOptions viewOptions, surface_type ***data,
                                 MemoryVisibilityGrid *visgrid) {
    _LOG_INFO(_("Computing events..."));

    /*buffer to hold 3 rows */
    Eigen::Matrix<surface_type, 3, Eigen::Dynamic> inrast;
    inrast.resize(3, raster.cols());
    inrast.fill(-1);
    _LOG_F(2, "inrast has shape (%d, %d)", inrast.rows(), inrast.cols());

    int nrows = raster.rows();
    int ncols = raster.cols();
    _LOG_VERBOSE("raster has shape (%d, %d)", nrows, ncols);

    /*alloc data ; data is used to store all the cells on the same row
       as the viewpoint. */
    *data = (surface_type **)_MALLOC(3 * sizeof(surface_type *));
    assert(*data);
    (*data)[0] = (surface_type *)_MALLOC(3 * ncols * sizeof(surface_type));
    assert((*data)[0]);
    (*data)[1] = (*data)[0] + ncols;
    (*data)[2] = (*data)[1] + ncols;

    int isnull = 0;

    /*keep track of the number of events added, to be returned later */
    size_t nevents = 0;

    /*scan through the raster data */
    dimensionType i, j;
    double ax, ay;
    AEvent e;

    /* read first row */
    inrast.row(2) = raster.row(0).array();

    e.angle = -1;
    for (i = 0; i < nrows; i++) {
        /*read in the raster row */

        _LOG_SCOPE_F(1, "Reading row %d", i);

        Eigen::Matrix<surface_type, 1, -1> tmprast = inrast.row(0);
        inrast.row(0) = inrast.row(1);
        inrast.row(1) = inrast.row(2);
        inrast.row(2) = tmprast;

        if (i < nrows - 1)
            inrast.row(2) = raster.row(i + 1);
        else
            inrast.row(2).array() = -1;

        /*fill event list with events from this row */
        for (j = 0; j < ncols; j++) {
            e.row = i;
            e.col = j;

            /*read the elevation value into the event */
            isnull = inrast.row(1)[j] == -1;
            e.elev[1] = inrast.row(1)[j];

            /*write it into the row of data going through the viewpoint */
            if (i == vp->row) {
                (*data)[0][j] = e.elev[1];
                (*data)[1][j] = e.elev[1];
                (*data)[2][j] = e.elev[1];
            }

            /* set the viewpoint, and don't insert it into eventlist */
            if (i == vp->row && j == vp->col) {
                set_viewpoint_elev(vp, e.elev[1] + viewOptions.obsElev);
                if (viewOptions.tgtElev > 0)
                    vp->target_offset = viewOptions.tgtElev;
                else
                    vp->target_offset = 0.;
                if (isnull) {
                    /*what to do when viewpoint is NODATA ? */
                    _LOG_WARNING(_("Viewpoint is NODATA."));
                    _LOG_INFO(_("Will assume its elevation is = %f"), vp->elev);
                }

                add_result_to_inmem_visibilitygrid(visgrid, i, j, 180);
                continue;
            }

            /*don't insert in eventlist nodata cell events */
            if (isnull) {
                /* record this cell as being NODATA; this is necessary so
                   that we can distingush invisible events, from nodata
                   events in the output */
                add_result_to_inmem_visibilitygrid(visgrid, i, j,
                                                   hd->nodata_value);
                continue;
            }
            /* if point is outside maxDist, do NOT include it as an
               event */
            if (is_point_outside_max_dist(*vp, *hd, i, j, viewOptions.maxDist))
                continue;

            /* if it got here it is not the viewpoint, not NODATA, and
               within max distance from viewpoint; generate its 3 events
               and insert them */

            /* get ENTER elevation */
            e.eventType = ENTERING_EVENT;
            e.elev[0] = calculate_event_elevation(e, nrows, ncols, vp->row,
                                                  vp->col, inrast);

            /* get EXIT elevation */
            e.eventType = EXITING_EVENT;
            e.elev[2] = calculate_event_elevation(e, nrows, ncols, vp->row,
                                                  vp->col, inrast);

            /*write adjusted elevation into the row of data going through the
             * viewpoint */
            if (i == vp->row) {
                (*data)[0][j] = e.elev[0];
                (*data)[1][j] = e.elev[1];
                (*data)[2][j] = e.elev[2];
            }

            /*put event into event list */
            e.eventType = ENTERING_EVENT;
            calculate_event_position(e, vp->row, vp->col, &ay, &ax);
            e.angle = calculate_angle(ax, ay, vp->col, vp->row);
            eventList[nevents] = e;
            nevents++;

            e.eventType = CENTER_EVENT;
            calculate_event_position(e, vp->row, vp->col, &ay, &ax);
            e.angle = calculate_angle(ax, ay, vp->col, vp->row);
            eventList[nevents] = e;
            nevents++;

            e.eventType = EXITING_EVENT;
            calculate_event_position(e, vp->row, vp->col, &ay, &ax);
            e.angle = calculate_angle(ax, ay, vp->col, vp->row);
            eventList[nevents] = e;
            nevents++;
        }
    }

    return nevents;
}