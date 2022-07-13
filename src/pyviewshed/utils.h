//
// Created by jhewers on 13/07/22.
//

#ifndef PYVIEWSHED_UTILS_H
#define PYVIEWSHED_UTILS_H

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include "logging.h"

template <typename T>
Eigen::Matrix<T, -1, -1> arr_to_eigen(T* array, int nrows, int ncols) {
    static_assert(std::is_floating_point<T>::value,
                  "array type must be a floating point type");

    _LOG_SCOPE_FUNCTION(1);
    _LOG_VERBOSE("Input shape (%d, %d)", nrows, ncols);

    Eigen::Map<Eigen::Matrix<T, -1, -1, Eigen::RowMajor>> output(array, nrows,
                                                                 ncols);

    _LOG_VERBOSE("Output shape (%d, %d)", output.rows(), output.cols());
    return output;
}

static std::string toString(Eigen::MatrixXf mat) {
    std::stringstream ss;
    ss << mat;
    return ss.str();
}

#endif  // PYVIEWSHED_UTILS_H
