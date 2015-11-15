cmake_minimum_required (VERSION 2.6)
project(SOLO_bksolution)

find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})
set(LIBS ${LIBS} ${GSL_LIBRARIES})

set(LIBS ${LIBS} m)

include_directories(${interp2d_SOURCE_DIR} ${quasimontecarlo_SOURCE_DIR})

add_library(bksol discretize.c solver.cpp)
target_link_libraries(bksol ${LIBS})
