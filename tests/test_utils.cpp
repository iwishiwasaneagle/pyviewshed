//
// Created by jhewers on 13/07/22.
//

#include <Eigen/Core>
#include <catch2/catch.hpp>

#include "utils.h"

SCENARIO("2D c-style arrays can be converted to eigen matrices", "[eigen]") {
    GIVEN("A (2,3) shaped float array") {
        float arr[2][3] = {{0.0f, 1.0f, 2.0f}, {3.0f, 4.0f, 5.0f}};
        Eigen::Matrix<float, 2, 3> expected;
        expected << 0, 1, 2, 3, 4, 5;
        CAPTURE(expected, arr);

        WHEN("The array is converted to an eigen matrix") {
            Eigen::MatrixXf output = arr_to_eigen(&(arr[0][0]), 2, 3);
            CAPTURE(output);
            THEN("The resultant data is equivalent") {
                REQUIRE(output == expected);
            }
        }
    }
    GIVEN("A (20,30) shaped double array") {
        Eigen::Array<double, 20, 30> expected =
            Eigen::Array<double, 20, 30>::Random(20, 30);
        CAPTURE(expected.size());
        double data[20][30];
        for (int i = 0; i < 20; i++) {
            for (int j = 0; j < 30; j++) {
                data[i][j] = expected(i, j);
            }
        }

        WHEN("The array is converted to an eigen matrix") {
            Eigen::Matrix<double, 20, 30> output =
                arr_to_eigen(&(data[0][0]), 20, 30);
            CAPTURE(output.size());
            THEN("The resultant data is equivalent") {
                double a, b;
                for (int i = 0; i < 20; i++) {
                    for (int j = 0; j < 30; j++) {
                        a = output(i, j);
                        b = expected(i, j);
                        CAPTURE(a, b);
                        CHECK(a == b);
                    }
                }
            }
        }
    }
    GIVEN("A (5,4) shaped double array") {
        static const short n = 5;
        static const short m = 4;

        Eigen::Array<double, n, m> expected;
        expected << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17,
            18, 19, 20;
        CAPTURE(expected, expected.rows(), expected.cols());
        double data[n][m];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                data[i][j] = expected(i, j);
                CAPTURE(data[i][j]);
            }
        }

        WHEN("The array is converted to an eigen matrix") {
            Eigen::Matrix<double, n, m> output =
                arr_to_eigen(&(data[0][0]), n, m);
            CAPTURE(output);
            THEN("The resultant data is equivalent at (0,0)") {
                REQUIRE(output(0, 0) == 1);
                REQUIRE(expected(0, 0) == 1);
                REQUIRE(data[0][0] == 1);
            }
            THEN("The resultant data is equivalent at (3,4)") {
                REQUIRE(output(4, 3) == 20);
                REQUIRE(expected(4, 3) == 20);
                REQUIRE(data[4][3] == 20);
            }
            THEN("The resultant data is equivalent at (2,2)") {
                REQUIRE(output(2, 2) == 11);
                REQUIRE(expected(2, 2) == 11);
                REQUIRE(data[2][2] == 11);
            }
            THEN("The resultant array is the same") {
                double a, b;
                for (int i = 0; i < n; i++) {
                    for (int j = 0; j < m; j++) {
                        a = output(i, j);
                        b = expected(i, j);
                        CAPTURE(a, b);
                        CHECK(a == b);
                    }
                }
            }
        }
    }
}