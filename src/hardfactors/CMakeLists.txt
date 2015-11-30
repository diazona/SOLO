cmake_minimum_required (VERSION 2.6)

find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})
set(LIBS ${LIBS} ${GSL_LIBRARIES})

find_package(MuParser REQUIRED)
include_directories(${MUPARSER_INCLUDE_DIRS})
set(LIBS ${LIBS} ${MUPARSER_LIBRARIES})

set(LIBS ${LIBS} m)

include_directories(${gslmuparser_SOURCE_DIR} ${interp2d_SOURCE_DIR})

add_executable(hfparser
    hfparser.cpp
    hardfactor.cpp
    hardfactor_parser.cpp
    ${SOLO_SOURCE_DIR}/integration/integrationtype.cpp
    ${SOLO_SOURCE_DIR}/integration/integrationcontext.cpp
    ${SOLO_SOURCE_DIR}/mstwpdf.cc
    ${SOLO_SOURCE_DIR}/utils/utils.cpp)
target_link_libraries(hfparser gslmuparser interp2d dsspinlo ${LIBS})

install(TARGETS hfparser
 RUNTIME DESTINATION bin
 LIBRARY DESTINATION lib
 ARCHIVE DESTINATION lib
 COMPONENT auxiliary_executables
 CONFIGURATIONS Debug Release
)
install(FILES hardfactors.cfg
 DESTINATION config
 COMPONENT configuration
 CONFIGURATIONS Debug Release
)
