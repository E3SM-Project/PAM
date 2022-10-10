
add_library(dycore compute_time_step.cpp 
                   fct_positivity.cpp
                   riemann/riemann_rho_theta_full.cpp
                   tendencies/rho_theta/tendencies_rho_theta.cpp)
target_include_directories(dycore PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(dycore PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/riemann)
target_include_directories(dycore PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/ader_differential_transforms)
target_include_directories(dycore PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/tendencies/rho_theta)
include (${YAKL_HOME}/yakl_utils.cmake)
yakl_process_target(dycore)
target_link_libraries(dycore pam_core)

