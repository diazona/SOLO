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
  ${SOLO_SOURCE_DIR}/mstwpdf.cc
  ${SOLO_SOURCE_DIR}/coupling.cpp
  ${SOLO_SOURCE_DIR}/factorizationscale.cpp
  ${SOLO_SOURCE_DIR}/configuration/context.cpp
  ${SOLO_SOURCE_DIR}/configuration/configuration.cpp
  ${SOLO_SOURCE_DIR}/utils/utils.cpp)
target_link_libraries(gluondisteval dsspinlo interp2d ${LIBS})

install(TARGETS gdist
 RUNTIME DESTINATION bin
 LIBRARY DESTINATION lib
 ARCHIVE DESTINATION lib
 COMPONENT libraries
 CONFIGURATIONS Debug Release
)
install(FILES gluondist.h
 DESTINATION include/SOLO
 COMPONENT libraries
 CONFIGURATIONS Debug Release
)
install(TARGETS gluondisteval
 RUNTIME DESTINATION bin
 LIBRARY DESTINATION lib
 ARCHIVE DESTINATION lib
 COMPONENT auxiliary_executables
 CONFIGURATIONS Debug Release
)
