
set(PAM_CORE_SRC mmf_interface/mmf_interface_globals.cpp
                 mmf_interface/mmf_interface_extern_c.cpp
                 mmf_interface/mmf_fortran_interface.F90
                 pam_coupler_globals.cpp)

add_library(pam_core ${PAM_CORE_SRC})
target_include_directories(pam_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(pam_core PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
target_include_directories(pam_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/mmf_interface)
target_include_directories(pam_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/mmf_routines)
set_source_files_properties(mmf_interface/mmf_fortran_interface.F90 PROPERTIES COMPILE_FLAGS "${YAKL_F90_FLAGS}")
include(${YAKL_HOME}/yakl_utils.cmake)
yakl_process_target(pam_core)

