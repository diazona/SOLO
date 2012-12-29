cmake_minimum_required (VERSION 2.6)
project (oneloopcalc)

find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})
set(LIBS ${LIBS} ${GSL_LIBRARIES})

set(LIBS ${LIBS} m)

add_executable(oneloopcalc oneloopcalc.cpp integrator.cpp hardfactors_position.cpp hardfactors_momentum.cpp dss_pinlo.cpp mstwpdf.cc cubature.c gluondist.cpp integrationcontext.cpp integrationtype.cpp context.cpp utils.cpp)
target_link_libraries(oneloopcalc ${CMAKE_SOURCE_DIR}/libinterp2d.a ${LIBS})

add_executable(contextreader context.cpp utils.cpp gluondist.cpp)
target_link_libraries(contextreader ${CMAKE_SOURCE_DIR}/libinterp2d.a ${LIBS})
set_target_properties(contextreader PROPERTIES COMPILE_DEFINITIONS CONTEXT_TEST)

add_executable(utils_test utils.cpp)
set_target_properties(utils_test PROPERTIES COMPILE_DEFINITIONS UTILS_TEST)
