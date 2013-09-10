cmake_minimum_required (VERSION 2.6)
project (oneloopcalc)

find_package(OpenSSL REQUIRED)
include_directories(${OPENSSL_INCLUDE_DIR})
set(LIBS ${LIBS} ${OPENSSL_LIBRARIES})

find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})
set(LIBS ${LIBS} ${GSL_LIBRARIES})

set(LIBS ${LIBS} m)

add_custom_target(git_revision.h
    git log -1 "--format=format:#define GIT_REVISION \"%H\"%n" HEAD > git_revision.h 
    COMMAND git diff-index --quiet HEAD || echo "#define GIT_DIRTY" >> git_revision.h
    WORKING_DIRECTORY ${oneloopcalc_SOURCE_DIR} VERBATIM)

include_directories(${interp2d_SOURCE_DIR} ${quasimontecarlo_SOURCE_DIR})

add_executable(oneloopcalc oneloopcalc.cpp integrator.cpp hardfactors_position.cpp hardfactors_radial.cpp hardfactors_momentum.cpp dss_pinlo.cpp mstwpdf.cc cubature.c gluondist.cpp coupling.cpp factorizationscale.cpp integrationcontext.cpp integrationtype.cpp context.cpp utils.cpp)
target_link_libraries(oneloopcalc interp2d quasimontecarlo ${LIBS})
add_dependencies(oneloopcalc git_revision.h)

add_executable(contextreader context.cpp utils.cpp coupling.cpp gluondist.cpp factorizationscale.cpp)
target_link_libraries(contextreader interp2d ${LIBS})
set_target_properties(contextreader PROPERTIES COMPILE_DEFINITIONS CONTEXT_TEST)

add_executable(utils_test EXCLUDE_FROM_ALL utils.cpp)
set_target_properties(utils_test PROPERTIES COMPILE_DEFINITIONS UTILS_TEST)

add_executable(dss_pinlo_test EXCLUDE_FROM_ALL dss_pinlo_test.cpp dss_pinlo.cpp)
target_link_libraries(dss_pinlo_test interp2d ${LIBS})

add_executable(gluondist gluondist_driver.cpp gluondist.cpp context.cpp utils.cpp coupling.cpp factorizationscale.cpp)
target_link_libraries(gluondist interp2d ${LIBS})

add_executable(kovr kov-position/gauleg.cpp kov-position/interr.cpp kov-position/kovr.cpp)
target_link_libraries(kovr m)

set_target_properties(oneloopcalc contextreader utils_test dss_pinlo_test gluondist kovr PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
