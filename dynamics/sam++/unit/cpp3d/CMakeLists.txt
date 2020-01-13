
add_executable(cpp3d ../dmdf.F90 ../cpp_driver.F90 ${CPP_SRC})
target_link_libraries(cpp3d ${NCFLAGS})
set_property(TARGET cpp3d APPEND PROPERTY COMPILE_FLAGS ${DEFS3D} )

