
set(PAM_CORE_SRC pam_interface/pam_interface_globals.cpp
                 pam_interface/pam_interface_extern_c.cpp
                 pam_interface/pam_fortran_interface.F90
                 pam_coupler_globals.cpp)

add_library(pam_core ${PAM_CORE_SRC})
target_include_directories(pam_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(pam_core PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
target_include_directories(pam_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/pam_interface)
target_include_directories(pam_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/modules)
set_source_files_properties(pam_interface/pam_fortran_interface.F90 PROPERTIES COMPILE_FLAGS "${YAKL_F90_FLAGS}")
include(${YAKL_HOME}/yakl_utils.cmake)
yakl_process_target(pam_core)

