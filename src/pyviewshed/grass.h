
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

#ifndef _GRASS_H
#define _GRASS_H

#include <math.h>

#include <Eigen/Core>

#include "eventlist.h"
#include "grid.h"
#include "visibility.h"

/* calculate ENTER and EXIT event elevation */
surface_type calculate_event_elevation(AEvent e, int nrows, int ncols,
                                       dimensionType vprow, dimensionType vpcol,
                                       Eigen::Matrix<float, 3, -1> raster);

/*  ************************************************************ */
/* input: an array capable to hold the max number of events, a raster
   name, a viewpoint and the viewOptions; action: figure out all events
   in the input file, and write them to the event list. data is
   allocated and initialized with all the cells on the same row as the
   viewpoint. it returns the number of events. initialize and fill
   AEvent* with all the events for the map.  Used when solving in
   memory, so the AEvent* should fit in memory.  */
size_t init_event_list_in_memory(AEvent *eventList, Eigen::MatrixXf,
                                 Viewpoint *vp, GridHeader *hd,
                                 ViewOptions viewOptions, surface_type ***data,
                                 MemoryVisibilityGrid *visgrid);

#endif /*_GRASS_H*/
