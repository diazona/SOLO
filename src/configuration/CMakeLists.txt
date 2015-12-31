cmake_minimum_required (VERSION 2.6)

find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})
set(LIBS ${LIBS} ${GSL_LIBRARIES})

set(LIBS ${LIBS} m)

include_directories(${interp2d_SOURCE_DIR} ${quasimontecarlo_SOURCE_DIR})

add_executable(contextreader
  contextreader.cpp
  context.cpp
  configuration.cpp
  ${SOLO_SOURCE_DIR}/mstwpdf.cc
  ${SOLO_SOURCE_DIR}/utils/utils.cpp
  ${SOLO_SOURCE_DIR}/coupling.cpp
  ${SOLO_SOURCE_DIR}/factorizationscale.cpp)
target_link_libraries(contextreader dsspinlo interp2d gdist ${LIBS})
set_target_properties(contextreader PROPERTIES COMPILE_DEFINITIONS CONTEXT_TEST)

install(TARGETS contextreader
 RUNTIME DESTINATION bin
 LIBRARY DESTINATION lib
 ARCHIVE DESTINATION lib
 COMPONENT auxiliary_executables
 CONFIGURATIONS Debug Release
)