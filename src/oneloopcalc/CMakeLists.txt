cmake_minimum_required (VERSION 2.6)
project (SOLO_oneloopcalc)

find_package(OpenSSL REQUIRED)
include_directories(${OPENSSL_INCLUDE_DIR})
set(LIBS ${LIBS} ${OPENSSL_LIBRARIES})

find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})
set(LIBS ${LIBS} ${GSL_LIBRARIES})

find_package(MuParser REQUIRED)
include_directories(${MUPARSER_INCLUDE_DIRS})
set(LIBS ${LIBS} ${MUPARSER_LIBRARIES})

set(LIBS ${LIBS} m)

find_package(Git)
if(GIT_FOUND)
    add_custom_target(git_revision.h
        ${GIT_EXECUTABLE} log -1 "--format=format:#define GIT_REVISION \"%H\"%n" HEAD > git_revision.h 
        COMMAND git diff-index --quiet HEAD || echo "#define GIT_DIRTY" >> git_revision.h
        WORKING_DIRECTORY ${SOLO_oneloopcalc_BINARY_DIR} VERBATIM)
else()
    add_custom_target(git_revision.h echo >> git_revision.h WORKING_DIRECTORY ${SOLO_oneloopcalc_BINARY_DIR} VERBATIM)
endif()

include_directories(${gslmuparser_SOURCE_DIR} ${interp2d_SOURCE_DIR} ${quasimontecarlo_SOURCE_DIR})

add_executable(oneloopcalc
    oneloopcalc.cpp
    ${SOLO_SOURCE_DIR}/mstwpdf.cc
    ${SOLO_SOURCE_DIR}/coupling.cpp
    ${SOLO_SOURCE_DIR}/factorizationscale.cpp
    ${SOLO_SOURCE_DIR}/configuration/context.cpp
    ${SOLO_SOURCE_DIR}/hardfactors/hardfactor_parsed.cpp
    ${SOLO_SOURCE_DIR}/hardfactors/hardfactor_parser.cpp
    ${SOLO_SOURCE_DIR}/hardfactors/hardfactor_group_parser.cpp
    ${SOLO_SOURCE_DIR}/hardfactors/hardfactors_position.cpp
    ${SOLO_SOURCE_DIR}/hardfactors/hardfactors_radial.cpp
    ${SOLO_SOURCE_DIR}/hardfactors/hardfactors_momentum.cpp
    ${SOLO_SOURCE_DIR}/integration/cubature.c
    ${SOLO_SOURCE_DIR}/integration/integrationcontext.cpp
    ${SOLO_SOURCE_DIR}/integration/integrationtype.cpp
    ${SOLO_SOURCE_DIR}/integration/integrator.cpp
    ${SOLO_SOURCE_DIR}/utils/utils.cpp)
target_link_libraries(oneloopcalc gslmuparser interp2d quasimontecarlo dsspinlo gdist ${LIBS})
add_dependencies(oneloopcalc git_revision.h)
