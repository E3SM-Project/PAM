
add_library(dycore compute_time_step.cpp 
                   fct_positivity.cpp
                   riemann/riemann_rho_theta_full_x.cpp
                   riemann/riemann_rho_theta_full_y.cpp
                   riemann/riemann_rho_theta_full_z.cpp
                   riemann/riemann_rho_theta_full_xyz.cpp
                   tendencies/rho_theta/tendencies_rho_theta.cpp
                   tendencies/rho_theta/tendencies_rho_theta_compressible_explicit_ader_x.cpp
                   tendencies/rho_theta/tendencies_rho_theta_compressible_explicit_ader_y.cpp
                   tendencies/rho_theta/tendencies_rho_theta_compressible_explicit_ader_z.cpp
                   tendencies/rho_theta/tendencies_rho_theta_compressible_explicit_semidiscrete_xyz.cpp)
target_include_directories(dycore PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(dycore PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/riemann)
target_include_directories(dycore PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/ader_differential_transforms)
target_include_directories(dycore PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/tendencies/rho_theta)
include (${YAKL_HOME}/yakl_utils.cmake)
yakl_process_target(dycore)
target_link_libraries(dycore pam_core)

