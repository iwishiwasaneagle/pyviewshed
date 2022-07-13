find_package(Catch2 REQUIRED)
message(STATUS "Found Catch2")

get_target_property(PROJECT_SOURCES pyviewshed SOURCES)
list(REMOVE_ITEM PROJECT_SOURCES src/bindings.cpp)
list(REMOVE_ITEM PROJECT_SOURCES ${PROJECT_SOURCE_DIR}/src/main.cpp)

set(TEST_SOURCES
        test_main.cpp
        test_utils.cpp)

add_executable(test_pyviewshed ${TEST_SOURCES} ${PROJECT_SOURCES})

target_link_libraries(test_pyviewshed pybind11::module Catch2::Catch2 --coverage )
target_compile_options(test_pyviewshed PRIVATE --coverage)
target_include_directories(test_pyviewshed PRIVATE ${SOURCE_DIR})

include(CTest)
include(Catch)
catch_discover_tests(test_pyviewshed)
