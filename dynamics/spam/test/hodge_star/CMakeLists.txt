
add_executable(hodge_star hodge_star.cpp)
yakl_process_target(hodge_star)
target_link_libraries(hodge_star dycore pam_core ${PAM_LINK_FLAGS})

add_test(NAME hodge_star_test COMMAND ./hodge_star)

add_executable(hodge_star_extruded hodge_star_extruded.cpp)
yakl_process_target(hodge_star_extruded)
target_link_libraries(hodge_star_extruded dycore pam_core ${PAM_LINK_FLAGS})
target_compile_definitions(hodge_star_extruded PUBLIC _EXTRUDED)

add_test(NAME extruded_hodge_star_test COMMAND ./hodge_star_extruded)
