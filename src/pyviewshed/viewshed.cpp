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
 ****************************************************************************/

#include "viewshed.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include <cstdlib>

#include "defs.h"
#include "eventlist.h"
#include "grass.h"
#include "logging.h"
#include "statusstructure.h"
#include "visibility.h"

#define VIEWSHEDDEBUG if (0)
#define INMEMORY_DEBUG if (0)

/* ------------------------------------------------------------ */
static void print_statusnode(StatusNode sn) {
    _LOG_F(3, "processing (row=%d, col=%d, dist=%f, grad=%f)", sn.row, sn.col,
           sn.dist2vp, sn.gradient[1]);
    return;
}

/* ------------------------------------------------------------ */
/* allocates the eventlist array used by kreveled_in_memory; it is
   possible that the amount of memory required is more than that
   supported by the platform; for e.g. on a 32-bt platform cannot
   allocate more than 4GB. Try to detect this situation.  */
AEvent *allocate_eventlist(GridHeader *hd) {
    AEvent *eventList;

    long long totalsize = hd->ncols * hd->nrows * 3;

    totalsize *= sizeof(AEvent);
    _LOG_F(1, "total size of eventlist is %lld B (%d MB);  ", totalsize,
           (int)(totalsize >> 20));

    /* what's the size of size_t on this machine? */
    int sizet_size;

    sizet_size = (int)sizeof(size_t);
    _LOG_F(1, "size_t is %d B", sizet_size);

    if (sizet_size >= 8) {
        _LOG_F(1, "64-bit platform, great.");
    } else {
        /* this is the max value of size_t */
        long long m = ((long long)1 << (sizeof(size_t) * 8 - 1)),
                  maxsizet = m - 1;

        maxsizet += m;

        _LOG_F(1, "max size_t is %lld", maxsizet);

        /* checking whether allocating totalsize causes an overflow */
        if (totalsize > maxsizet) {
            _ABORT_F(
                _("Running the program in-memory mode requires "
                  "memory beyond the capability of the platform. "
                  "Use external mode, or a 64-bit platform."));
        }
    }

    _LOG_F(1, "allocating eventList...");
    eventList = (AEvent *)_MALLOC(totalsize);

    assert(eventList);
    _LOG_F(1, "...ok");

    return eventList;
}

/*///////////////////////////////////////////////////////////
   ------------------------------------------------------------ run
   Viewshed's sweep algorithm on the grid stored in the given file, and
   with the given viewpoint.  Create a visibility grid and return
   it. The computation runs in memory, which means the input grid, the
   status structure and the output grid are stored in arrays in
   memory.


   The output: A cell x in the visibility grid is recorded as follows:

   if it is NODATA, then x  is set to NODATA
   if it is invisible, then x is set to INVISIBLE
   if it is visible,  then x is set to the vertical angle wrt to viewpoint

 */
float **viewshed_in_memory(Eigen::MatrixXf raster, GridHeader *hd,
                           Viewpoint *vp, ViewOptions viewOptions) {
    _LOG_VERBOSE(_("Start sweeping."));

    /* ------------------------------ */
    /* create the visibility grid  */
    MemoryVisibilityGrid *visgrid;

    visgrid = create_inmem_visibilitygrid(*hd, *vp);
    /* set everything initially invisible */
    set_inmem_visibilitygrid(visgrid, INVISIBLE);
    assert(visgrid);
    _LOG_F(1, "visibility grid size:  %d x %d x %d B (%d MB)", hd->nrows,
           hd->ncols, (int)sizeof(float),
           (int)(((long long)(hd->nrows * hd->ncols * sizeof(float))) >> 20));

    /* ------------------------------ */
    /* construct the event list corresponding to the given input file
       and viewpoint; this creates an array of all the cells on the
       same row as the viewpoint */
    surface_type **data;
    size_t nevents;

    AEvent *eventList = allocate_eventlist(hd);

    nevents = init_event_list_in_memory(eventList, raster, vp, hd, viewOptions,
                                        &data, visgrid);

    _LOG_F(1, "actual nb events is %lu", (long unsigned int)nevents);

    /* ------------------------------ */
    /*sort the events radially by angle */
    _LOG_VERBOSE(_("Sorting events..."));
    fflush(stdout);

    /*this is recursive and seg faults for large arrays
       //qsort(eventList, nevents, sizeof(AEvent), radial_compare_events);

       //this is too slow...
       //heapsort(eventList, nevents, sizeof(AEvent), radial_compare_events);

       //iostream quicksort */
    std::qsort(eventList, nevents, sizeof(AEvent), radial_compare_events);

    _LOG_VERBOSE("Done.");
    fflush(stdout);

    /* ------------------------------ */
    /*create the status structure */
    StatusList *status_struct = create_status_struct();

    /*Put cells that are initially on the sweepline into status structure */
    StatusNode sn;

    for (dimensionType i = vp->col + 1; i < hd->ncols; i++) {
        AEvent e;
        double ax, ay;

        sn.col = i;
        sn.row = vp->row;
        e.col = i;
        e.row = vp->row;
        e.elev[0] = data[0][i];
        e.elev[1] = data[1][i];
        e.elev[2] = data[2][i];

        if (!is_nodata(data[1][i]) &&
            !is_point_outside_max_dist(*vp, *hd, sn.row, sn.col,
                                       viewOptions.maxDist)) {
            /*calculate Distance to VP and Gradient, store them into sn */
            /* need either 3 elevation values or
             * 3 gradients calculated from 3 elevation values */
            /* need also 3 angles */
            e.eventType = ENTERING_EVENT;
            calculate_event_position(e, vp->row, vp->col, &ay, &ax);
            sn.angle[0] = calculate_angle(ax, ay, vp->col, vp->row);
            calculate_event_gradient(&sn, 0, ay, ax, e.elev[0], vp, *hd);

            e.eventType = CENTER_EVENT;
            calculate_event_position(e, vp->row, vp->col, &ay, &ax);
            sn.angle[1] = calculate_angle(ax, ay, vp->col, vp->row);
            calculate_dist_n_gradient(&sn, e.elev[1], vp, *hd);

            e.eventType = EXITING_EVENT;
            calculate_event_position(e, vp->row, vp->col, &ay, &ax);
            sn.angle[2] = calculate_angle(ax, ay, vp->col, vp->row);
            calculate_event_gradient(&sn, 2, ay, ax, e.elev[2], vp, *hd);

            assert(sn.angle[1] == 0);

            if (sn.angle[0] > sn.angle[1]) sn.angle[0] -= 2 * M_PI;

            _LOG_F(2, "inserting: ");
            print_statusnode(sn);
            /*insert sn into the status structure */
            insert_into_status_struct(sn, status_struct);
        }
    }
    _FREE(data[0]);
    _FREE(data);

    /* ------------------------------ */
    /*sweep the event list */
    long nvis = 0; /*number of visible cells */
    AEvent *e;

    _LOG_INFO(_("Computing visibility..."));

    for (size_t i = 0; i < nevents; i++) {
        /*get out one event at a time and process it according to its type */
        e = &(eventList[i]);

        sn.col = e->col;
        sn.row = e->row;
        // sn.elev = e->elev;

        /*calculate Distance to VP and Gradient */
        calculate_dist_n_gradient(&sn, e->elev[1] + vp->target_offset, vp, *hd);
        _LOG_F(3, "event: ");
        print_event(*e);
        _LOG_F(3, "sn.dist=%f, sn.gradient=%f", sn.dist2vp, sn.gradient[1]);

        switch (e->eventType) {
            case ENTERING_EVENT:
                double ax, ay;
                /*insert node into structure */
                _LOG_F(3, "..ENTER-EVENT: insert");

                /* need either 3 elevation values or
                 * 3 gradients calculated from 3 elevation values */
                /* need also 3 angles */
                calculate_event_position(*e, vp->row, vp->col, &ay, &ax);
                // sn.angle[0] = calculate_angle(ax, ay, vp->col, vp->row);
                sn.angle[0] = e->angle;
                calculate_event_gradient(&sn, 0, ay, ax, e->elev[0], vp, *hd);

                e->eventType = CENTER_EVENT;
                calculate_event_position(*e, vp->row, vp->col, &ay, &ax);
                sn.angle[1] = calculate_angle(ax, ay, vp->col, vp->row);
                calculate_dist_n_gradient(&sn, e->elev[1], vp, *hd);

                e->eventType = EXITING_EVENT;
                calculate_event_position(*e, vp->row, vp->col, &ay, &ax);
                sn.angle[2] = calculate_angle(ax, ay, vp->col, vp->row);
                calculate_event_gradient(&sn, 2, ay, ax, e->elev[2], vp, *hd);

                e->eventType = ENTERING_EVENT;

                if (e->angle < M_PI) {
                    if (sn.angle[0] > sn.angle[1]) sn.angle[0] -= 2 * M_PI;
                } else {
                    if (sn.angle[0] > sn.angle[1]) {
                        sn.angle[1] += 2 * M_PI;
                        sn.angle[2] += 2 * M_PI;
                    }
                }

                insert_into_status_struct(sn, status_struct);
                break;

            case EXITING_EVENT:
                /*delete node out of status structure */
                _LOG_F(3, "..EXIT-EVENT: delete");
                /* need only distance */
                delete_from_status_struct(status_struct, sn.dist2vp);
                break;

            case CENTER_EVENT:
                _LOG_F(3, "..QUERY-EVENT: query");
                /*calculate visibility */
                double max;

                /* consider current angle and gradient */
                max = find_max_gradient_in_status_struct(
                    status_struct, sn.dist2vp, e->angle, sn.gradient[1]);

                /*the point is visible: store its vertical angle  */
                if (max <= sn.gradient[1]) {
                    float vert_angle = get_vertical_angle(
                        *vp, sn, e->elev[1] + vp->target_offset, false);
                    _LOG_F(3, "vert_angle: %f", vert_angle);
                    add_result_to_inmem_visibilitygrid(visgrid, sn.row, sn.col,
                                                       vert_angle);
                    assert(vert_angle >= 0);
                    /* when you write the visibility grid you assume that
                       visible values are positive */
                    nvis++;
                }
                // else {
                /* cell is invisible */
                /*  the visibility grid is initialized all invisible */
                // visgrid->grid->grid_data[sn.row][sn.col] = INVISIBLE;
                // }
                break;
        }
    }

    _LOG_VERBOSE(_("Sweeping done."));
    _LOG_VERBOSE(
        _("Total cells %ld, visible cells %ld (%.1f percent)."),
        (long)visgrid->grid->hd->nrows * visgrid->grid->hd->ncols, nvis,
        (float)((float)nvis * 100 /
                (float)(visgrid->grid->hd->nrows * visgrid->grid->hd->ncols)));

    /*cleanup */
    _FREE(eventList);

    return visgrid->grid->grid_data;
}
