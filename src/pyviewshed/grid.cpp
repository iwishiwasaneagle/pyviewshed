
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

#include "grid.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "defs.h"

/* ------------------------------------------------------------ */
/*copy from b to a */
void copy_header(GridHeader *a, GridHeader b) {
    assert(a);
    a->nrows = b.nrows;
    a->ncols = b.ncols;
    a->xllcorner = b.xllcorner;
    a->yllcorner = b.yllcorner;
    a->ns_res = b.ns_res;
    a->ew_res = b.ew_res;
    a->nodata_value = b.nodata_value;
    return;
}

/* ------------------------------------------------------------ */
/*returns 1 if value is Nodata, 0 if it is not */
int is_nodata(float value) { return value == NOVALUE; }

/* ------------------------------------------------------------ */
/* create an empty grid and return it. The header and the data are set
   to NULL.  */
Grid *create_empty_grid() {
    Grid *ptr_grid = (Grid *)_MALLOC(sizeof(Grid));

    assert(ptr_grid);

    /*initialize structure */
    ptr_grid->hd = NULL;
    ptr_grid->grid_data = NULL;

#ifdef _DEBUG_ON
    printf("**DEBUG: createEmptyGrid \n");
    fflush(stdout);
#endif

    return ptr_grid;
}

/* ------------------------------------------------------------ */
/* allocate memory for grid_data; grid must have a header that gives
   the dimensions */
void alloc_grid_data(Grid *pgrid) {
    assert(pgrid);
    assert(pgrid->hd);

    pgrid->grid_data = (float **)_MALLOC(pgrid->hd->nrows * sizeof(float *));

    assert(pgrid->grid_data);

    dimensionType i;

    for (i = 0; i < pgrid->hd->nrows; i++) {
        pgrid->grid_data[i] =
            (float *)_MALLOC(pgrid->hd->ncols * sizeof(float));

        assert(pgrid->grid_data[i]);
    }

#ifdef _DEBUG_ON
    printf("**DEBUG: allocGridData\n");
    fflush(stdout);
#endif

    return;
}

/* ------------------------------------------------------------ */
/*destroy the structure and reclaim all memory allocated */
void destroy_grid(Grid *grid) {
    assert(grid);

    /*free grid data if its allocated */
    if (grid->grid_data) {
        dimensionType i;

        for (i = 0; i < grid->hd->nrows; i++) {
            if (!grid->grid_data[i]) _FREE((float *)grid->grid_data[i]);
        }

        _FREE((float **)grid->grid_data);
    }

    _FREE(grid->hd);
    _FREE(grid);

#ifdef _DEBUG_ON
    printf("**DEBUG: grid destroyed.\n");
    fflush(stdout);
#endif

    return;
}
