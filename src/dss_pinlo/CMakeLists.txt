cmake_minimum_required (VERSION 2.6)
project(dsspinlo)

find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})
set(LIBS ${LIBS} ${GSL_LIBRARIES})

set(LIBS ${LIBS} m)

include_directories(${interp2d_SOURCE_DIR})

add_library(dsspinlo dss_pinlo.cpp)
target_link_libraries(dsspinlo interp2d ${LIBS})

add_executable(dsspinlotest
  dss_pinlo_test.cpp
  dss_pinlo.cpp)
target_link_libraries(dsspinlotest interp2d ${LIBS})

install(TARGETS dsspinlo
 RUNTIME DESTINATION bin
 LIBRARY DESTINATION lib
 ARCHIVE DESTINATION lib
 COMPONENT libraries
 CONFIGURATIONS Debug Release
)
install(FILES dss_pinlo.h
 DESTINATION include/SOLO
 COMPONENT libraries
 CONFIGURATIONS Debug Release
)
install(TARGETS dsspinlotest
 RUNTIME DESTINATION bin
 LIBRARY DESTINATION lib
 ARCHIVE DESTINATION lib
 COMPONENT testing
 CONFIGURATIONS Debug
)
