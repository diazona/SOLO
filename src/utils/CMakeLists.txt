cmake_minimum_required (VERSION 2.6)

add_executable(utils_test utils.cpp)
set_target_properties(utils_test PROPERTIES COMPILE_DEFINITIONS UTILS_TEST)

install(TARGETS utils_test
 RUNTIME DESTINATION bin
 LIBRARY DESTINATION lib
 ARCHIVE DESTINATION lib
 COMPONENT testing
 CONFIGURATIONS Debug
)
