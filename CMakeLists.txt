cmake_minimum_required(VERSION 2.6)
project(SOLO)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

add_subdirectory(gslmuparser EXCLUDE_FROM_ALL)
add_subdirectory(interp2d EXCLUDE_FROM_ALL)
add_subdirectory(quasimontecarlo EXCLUDE_FROM_ALL)
add_subdirectory(src)
