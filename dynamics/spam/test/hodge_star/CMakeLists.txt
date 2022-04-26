
add_executable(hodge_star hodge_star.cpp)
yakl_process_target(hodge_star)
target_link_libraries(hodge_star dycore pam_core ${NCFLAGS})

add_test(NAME hodge_star_test COMMAND ./hodge_star)
