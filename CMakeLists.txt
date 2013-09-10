cmake_minimum_required(VERSION 2.6)
project(SOLO)

add_subdirectory(interp2d EXCLUDE_FROM_ALL)
add_subdirectory(quasimontecarlo EXCLUDE_FROM_ALL)
add_subdirectory(src)
