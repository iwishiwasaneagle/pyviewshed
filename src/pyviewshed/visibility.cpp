
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
#include "visibility.h"

#include <stdio.h>

#include <cassert>

#include "defs.h"
#include "grass.h"
#include "grid.h"
#include "logging.h"

/* ------------------------------------------------------------ */
/* viewpoint functions */
void print_viewpoint(Viewpoint vp) {
    _LOG_VERBOSE("vp=(%d, %d, %.1f) ", vp.row, vp.col, vp.elev);
    return;
}

/* ------------------------------------------------------------ */
void set_viewpoint_coord(Viewpoint *vp, dimensionType row, dimensionType col) {
    assert(vp);
    vp->row = row;
    vp->col = col;
    return;
}

/* ------------------------------------------------------------ */
void set_viewpoint_elev(Viewpoint *vp, float elev) {
    assert(vp);
    vp->elev = elev;
    return;
}

/* ------------------------------------------------------------ */
/*copy from b to a */
void copy_viewpoint(Viewpoint *a, Viewpoint b) {
    assert(a);
    a->row = b.row;
    a->col = b.col;
    a->elev = b.elev;
    return;
}

/* ------------------------------------------------------------ */
/* MemoryVisibilityGrid functions */

/* create and return a grid of the sizes specified in the header */
MemoryVisibilityGrid *create_inmem_visibilitygrid(const GridHeader &hd,
                                                  Viewpoint vp) {
    MemoryVisibilityGrid *visgrid;

    visgrid = (MemoryVisibilityGrid *)_MALLOC(sizeof(MemoryVisibilityGrid));

    assert(visgrid);

    /* create the grid  */
    visgrid->grid = create_empty_grid();
    assert(visgrid->grid);

    /* create the header */
    visgrid->grid->hd = (GridHeader *)_MALLOC(sizeof(GridHeader));

    assert(visgrid->grid->hd);

    /* set the header */
    copy_header(visgrid->grid->hd, hd);

    /* allocate the  Grid data */
    alloc_grid_data(visgrid->grid);

    /*allocate viewpoint */
    visgrid->vp = (Viewpoint *)_MALLOC(sizeof(Viewpoint));

    assert(visgrid->vp);
    copy_viewpoint(visgrid->vp, vp);

    return visgrid;
}

/* ------------------------------------------------------------ */
void free_inmem_visibilitygrid(MemoryVisibilityGrid *visgrid) {
    assert(visgrid);

    if (visgrid->grid) {
        destroy_grid(visgrid->grid);
    }
    if (visgrid->vp) {
        _FREE(visgrid->vp);
    }
    _FREE(visgrid);

    return;
}

/* ------------------------------------------------------------ */
/*set all values of visgrid's Grid to the given value */
void set_inmem_visibilitygrid(MemoryVisibilityGrid *visgrid, float val) {
    assert(visgrid && visgrid->grid && visgrid->grid->hd &&
           visgrid->grid->grid_data);

    dimensionType i, j;

    for (i = 0; i < visgrid->grid->hd->nrows; i++) {
        assert(visgrid->grid->grid_data[i]);
        for (j = 0; j < visgrid->grid->hd->ncols; j++) {
            visgrid->grid->grid_data[i][j] = val;
        }
    }
    return;
}

/* ------------------------------------------------------------ */
/*set the (i,j) value of visgrid's Grid to the given value */
void add_result_to_inmem_visibilitygrid(MemoryVisibilityGrid *visgrid,
                                        dimensionType i, dimensionType j,
                                        float val) {
    assert(visgrid && visgrid->grid && visgrid->grid->hd &&
           visgrid->grid->grid_data);
    assert(i < visgrid->grid->hd->nrows);
    assert(j < visgrid->grid->hd->ncols);
    assert(visgrid->grid->grid_data[i]);

    visgrid->grid->grid_data[i][j] = val;

    return;
}

/* ------------------------------------------------------------ */
/*  The following functions are used to convert the visibility results
   recorded during the viewshed computation into the output grid into
   tehe output required by the user.

   x is assumed to be the visibility value computed for a cell during the
   viewshed computation.

   The value computed during the viewshed is the following:

   x is NODATA if the cell is NODATA;


   x is INVISIBLE if the cell is invisible;

   x is the vertical angle of the cell wrt the viewpoint if the cell is
   visible---the angle is a value in (0,180).
 */
int is_visible(float x) {
    /* if GRASS is on, we cannot guarantee that NODATA is negative; so
       we need to check */
    int isnull = x == NOVALUE;

    if (isnull)
        return 0;
    else
        return (x >= 0);
}
int is_invisible_not_nodata(float x) { return ((int)x == (int)INVISIBLE); }

int is_invisible_nodata(float x) {
    return (!is_visible(x)) && (!is_invisible_not_nodata(x));
}

/* ------------------------------------------------------------ */
/* This function is called when the program runs in
   viewOptions.outputMode == OUTPUT_BOOL. */
float booleanVisibilityOutput(float x) {
    /* NODATA and INVISIBLE are both negative values */
    if (is_visible(x))
        return BOOL_VISIBLE;
    else
        return BOOL_INVISIBLE;
}

/* ------------------------------------------------------------ */
/* This function is called when the program runs in
   viewOptions.outputMode == OUTPUT_ANGLE. In this case x represents
   the right value.  */
float angleVisibilityOutput(float x) { return x; }
