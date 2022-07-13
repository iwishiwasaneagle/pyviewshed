#define EIGEN_STACK_ALLOCATION_LIMIT 0
#include <Eigen/Core>

#include "logging.h"
#include "utils.h"
#include "viewshed.h"
#include "visibility.h"

#define N 10000
#define M 5000

int main(int argc, char* argv[]) {
    loguru::init(argc, argv);

    static const Eigen::Matrix<float, N, M> input =
        Eigen::Matrix<float, N, M>::Random();

    _LOG_VERBOSE("Input size (%d, %d)", input.rows(), input.cols());

    Viewpoint vp(3, 3, 1.9, 1.7);
    Cell_head window(input.rows(), input.cols(), 1, 1, input.rows(),
                     input.cols(), 0, 0);
    GridHeader gh(window.cols, window.rows, window.west, window.south,
                  window.ew_res, window.ns_res, window);
    ViewOptions vo(vp.elev, vp.target_offset, 100);

    float** output = viewshed_in_memory(input, &gh, &vp, vo);

    Eigen::MatrixXf output_eig =
        arr_to_eigen<float>(&(output[0][0]), gh.nrows, gh.ncols);
    _LOG_VERBOSE("Output size(%d, %d)", output_eig.rows(), output_eig.cols());

    return 0;
}