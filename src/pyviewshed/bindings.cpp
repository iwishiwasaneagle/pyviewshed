#include <fmt/core.h>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>

#include "defs.h"
#include "grid.h"
#include "logging.h"
#include "viewshed.h"
#include "visibility.h"

Eigen::MatrixXf viewshed_in_to_memory(pybind11::EigenDRef<Eigen::MatrixXf> dem,
                                      Viewpoint *vp, ViewOptions vo,
                                      GridHeader *gh) {
    float **output = viewshed_in_memory(dem, gh, vp, vo);

    _LOG_F(2, "Initializing matrix");
    Eigen::MatrixXf retMatrix;
    retMatrix.resize(gh->nrows, gh->ncols);
    _LOG_F(2, "Initialized matrix with shape (%d,%d)", retMatrix.rows(),
           retMatrix.cols());
    for (int i = 0; i < gh->nrows; i++) {
        _LOG_SCOPE_F(5, "Converting row %d", i);
        for (int j = 0; j < gh->ncols; j++) {
            retMatrix(i, j) = output[i][j];
        }
    };
    return retMatrix;
}

using namespace pybind11::literals;
PYBIND11_MODULE(pyviewshed, m) {
    m.def("viewshed_in_to_memory", &viewshed_in_to_memory);

    pybind11::class_<ViewOptions>(m, "ViewOptions")
        .def(pybind11::init<float, float, float>())
        .def_readwrite("obsElev", &ViewOptions::obsElev)
        .def_readwrite("tgtElev", &ViewOptions::tgtElev)
        .def_readwrite("maxDist", &ViewOptions::maxDist)
        .def("__repr__", [](const ViewOptions &a) {
            return fmt::format(
                "<pyviewshed.ViewOptions(obsElev={:.2f},tgtElev={:."
                "2f},maxDist={:.2f})>",
                a.obsElev, a.tgtElev, a.maxDist);
        });

    pybind11::class_<Cell_head>(m, "Cell_head")
        .def(pybind11::init<int, int, double, double, double, double, double,
                            double>())
        .def_readwrite("rows", &Cell_head::rows)
        .def_readwrite("cols", &Cell_head::cols)
        .def_readwrite("ew_res", &Cell_head::ew_res)
        .def_readwrite("ns_res", &Cell_head::ns_res)
        .def_readwrite("north", &Cell_head::north)
        .def_readwrite("east", &Cell_head::east)
        .def_readwrite("south", &Cell_head::south)
        .def_readwrite("west", &Cell_head::west)
        .def("__repr__", [](const Cell_head &a) {
            return fmt::format(
                "<pyviewshed.Cell_head(rows={},cols={},ew"
                "_res={},ns_res={},N={},E={},S={},W={})>",
                a.rows, a.cols, a.ew_res, a.ns_res, a.north, a.east, a.south,
                a.west);
        });

    pybind11::class_<GridHeader>(m, "GridHeader")
        .def(pybind11::init<dimensionType, dimensionType, double, double,
                            double, double, Cell_head>())
        .def_readwrite("ncols", &GridHeader::ncols)
        .def_readwrite("nrows", &GridHeader::nrows)
        .def_readwrite("xllcorner", &GridHeader::xllcorner)
        .def_readwrite("yllcorner", &GridHeader::yllcorner)
        .def_readwrite("ew_res", &GridHeader::ew_res)
        .def_readwrite("ns_res", &GridHeader::ns_res)
        .def_readwrite("window", &GridHeader::window)
        .def("__repr__", [](const GridHeader &a) {
            return fmt::format(
                "<pyviewshed.GridHeader(ncols={},nrows={},xllcorner="
                "{},yllcorner={},ew_res={},ns_res={})>",
                a.ncols, a.nrows, a.xllcorner, a.yllcorner, a.ew_res, a.ns_res);
        });

    pybind11::class_<Viewpoint>(m, "Viewpoint")
        .def(pybind11::init<dimensionType, dimensionType, float, float>())
        .def_readwrite("row", &Viewpoint::row)
        .def_readwrite("col", &Viewpoint::col)
        .def_readwrite("elev", &Viewpoint::elev)
        .def_readwrite("target_offset", &Viewpoint::target_offset)
        .def("__repr__", [](const Viewpoint &a) {
            return fmt::format(
                "<pyviewshed.Viewpoint(row={},col={},elev={},target_offset={})"
                ">",
                a.row, a.col, a.elev, a.target_offset);
        });
}