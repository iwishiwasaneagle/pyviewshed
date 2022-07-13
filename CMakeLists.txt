cmake_minimum_required(VERSION 3.20)
project(pyviewshed VERSION 0.0.1)
set(CMAKE_CXX_STANDARD 20)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug)
endif()
set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
option(BUILD_MAIN "Build main instead of python bindings" OFF)
option(BUILD_TESTS "Build tests" ON)
set(SOURCE_DIR src/pyviewshed)

# --- Eigen3 --- #
find_package(Eigen3 3 REQUIRED)
message(STATUS "Found Eigen 3.4 at ${EIGEN3_INCLUDE_DIR}")
include_directories(${EIGEN3_INCLUDE_DIR})


# --- fmt --- #
add_subdirectory(extern/fmt)

# --- PyBind11 --- #
add_subdirectory(extern/pybind11)

# --- What binary to build --- #
if(BUILD_MAIN)
add_executable(pyviewshed ${PROJECT_SOURCE_DIR}/${SOURCE_DIR}/main.cpp)
else(BUILDMAIN)
pybind11_add_module(pyviewshed MODULE ${SOURCE_DIR}/bindings.cpp)
endif(BUILD_MAIN)

# --- Loguru --- #
include_directories(extern/loguru)
target_sources(pyviewshed PRIVATE ${PROJECT_SOURCE_DIR}/extern/loguru/loguru.cpp)
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)
target_link_libraries(pyviewshed PRIVATE Threads::Threads)

# --- Files --- #
include_directories(${SOURCE_DIR})
add_subdirectory(${SOURCE_DIR})

# --- Installation --- #
if(SKBUILD)
    set(INST_DIR .)
else()
    set(INST_DIR ${PROJECT_SOURCE_DIR}/${SOURCE_DIR})
endif()
message(STATUS "Installing into ${INST_DIR} from ${PROJECT_SOURCE_DIR} as SKBUILD=\"${SKBUILD}\"")
install(TARGETS ${MODULE_EXECUTABLE} DESTINATION ${INST_DIR})

# --- Tests --- #
if(BUILD_TESTS)
    include(CTest)
    set_property(GLOBAL PROPERTY USE_FOLDERS ON)
    add_subdirectory(tests)
endif()
