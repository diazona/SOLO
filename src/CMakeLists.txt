cmake_minimum_required (VERSION 2.6)
project (SOLO)

find_package(OpenSSL REQUIRED)
include_directories(${OPENSSL_INCLUDE_DIR})
set(LIBS ${LIBS} ${OPENSSL_LIBRARIES})

set(LIBS ${LIBS} m)

add_executable(kovr kov-position/gauleg.cpp kov-position/interr.cpp kov-position/kovr.cpp)
target_link_libraries(kovr m)

add_subdirectory(configuration)
add_subdirectory(dss_pinlo)
add_subdirectory(gluondist)
add_subdirectory(hardfactors)
add_subdirectory(oneloopcalc)
add_subdirectory(utils)
