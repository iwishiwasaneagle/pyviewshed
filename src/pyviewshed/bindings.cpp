#include <pybind11/eigen.h>

#include "defs.h"
#include "grid.h"
#include "logging.h"
#include "pybind11/pybind11.h"
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
            std::ostringstream stringStream;
            stringStream << "<pyviewshed.ViewOptions(obsElev=" << a.obsElev
                         << ",tgtElev=" << a.tgtElev << ",maxDist=" << a.maxDist
                         << ")>";
            return stringStream.str();
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
            std::ostringstream stringStream;
            stringStream << "<pyviewshed.Cell_head(rows=" << a.rows
                         << ",cols=" << a.cols << ",ew_res=" << a.ew_res
                         << ",ns_res=" << a.ns_res << ",north=" << a.north
                         << ",east=" << a.east << ",south" << a.south
                         << ",west=" << a.west << ")>";
            return stringStream.str();
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
            std::ostringstream stringStream;
            stringStream << "<pyviewshed.GridHeader(ncols=" << a.ncols
                         << ",nrows=" << a.nrows << ",xllcorner=" << a.xllcorner
                         << ",yllcorner=" << a.yllcorner
                         << ",ew_res=" << a.ew_res << ",ns_res=" << a.ns_res
                         << ")>";
            return stringStream.str();
        });

    pybind11::class_<Viewpoint>(m, "Viewpoint")
        .def(pybind11::init<dimensionType, dimensionType, float, float>())
        .def_readwrite("row", &Viewpoint::row)
        .def_readwrite("col", &Viewpoint::col)
        .def_readwrite("elev", &Viewpoint::elev)
        .def_readwrite("target_offset", &Viewpoint::target_offset)
        .def("__repr__", [](const Viewpoint &a) {
            std::ostringstream stringStream;
            stringStream << "<pyviewshed.Viewpoint(row=" << a.row
                         << ",col=" << a.col << ",elev=" << a.elev
                         << ",target_offset" << a.target_offset << ")>";
            return stringStream.str();
        });
}