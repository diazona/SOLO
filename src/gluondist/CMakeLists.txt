cmake_minimum_required (VERSION 2.6)
project(SOLO_gluondist)

find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})
set(LIBS ${LIBS} ${GSL_LIBRARIES})

set(LIBS ${LIBS} m)

include_directories(${interp2d_SOURCE_DIR} ${quasimontecarlo_SOURCE_DIR})

add_library(gdist gluondist.cpp)
target_link_libraries(gdist interp2d ${LIBS})

add_executable(gluondisteval
  gluondist_driver.cpp
  gluondist.cpp
  ../mstwpdf.cc
  ../coupling.cpp
  ../factorizationscale.cpp
  ../configuration/context.cpp
  ../utils/utils.cpp)
target_link_libraries(gluondisteval dsspinlo interp2d ${LIBS})
