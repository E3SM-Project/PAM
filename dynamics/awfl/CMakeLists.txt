
add_library(dycore compute_time_step.cpp)
target_include_directories(dycore PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
include (${YAKL_HOME}/yakl_utils.cmake)
yakl_process_target(dycore)
target_link_libraries(dycore pam_core)

