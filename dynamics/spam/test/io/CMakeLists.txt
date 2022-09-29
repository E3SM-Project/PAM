
add_executable(field_test field_test.cpp)
yakl_process_target(field_test)
target_link_libraries(field_test dycore pam_core ${PAM_LINK_FLAGS})

add_test(NAME run_field_test COMMAND ./field_test)
