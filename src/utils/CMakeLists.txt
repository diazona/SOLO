cmake_minimum_required (VERSION 2.6)

add_executable(utils_test EXCLUDE_FROM_ALL utils.cpp)
set_target_properties(utils_test PROPERTIES COMPILE_DEFINITIONS UTILS_TEST)

