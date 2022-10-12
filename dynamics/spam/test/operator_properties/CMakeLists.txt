
add_executable(exterior_derivative exterior_derivative.cpp)
yakl_process_target(exterior_derivative)
target_link_libraries(exterior_derivative dycore pam_core ${PAM_LINK_FLAGS})

add_test(NAME exterior_derivative_test COMMAND ./exterior_derivative)
