
set(SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/../../fortran/cloud.F90
            ${CMAKE_CURRENT_SOURCE_DIR}/../../fortran/ice_fall.F90
            ${CMAKE_CURRENT_SOURCE_DIR}/../../fortran/micro_diagnose.F90
            ${CMAKE_CURRENT_SOURCE_DIR}/../../fortran/micro.F90
            ${CMAKE_CURRENT_SOURCE_DIR}/../../fortran/precip_fall.F90
            ${CMAKE_CURRENT_SOURCE_DIR}/../../fortran/precip_init.F90
            ${CMAKE_CURRENT_SOURCE_DIR}/../../fortran/precip_proc.F90
            ${CMAKE_CURRENT_SOURCE_DIR}/../../fortran/sat.F90
            fortran_driver.F90 )


# Main driver
add_executable(fortran_driver ${SOURCES})
include_directories(${YAKL_BIN})
include_directories(${NF_INCLUDE})

include(${YAKL_HOME}/yakl_utils.cmake)
yakl_process_target(fortran_driver)

set_source_files_properties(${SOURCES} PROPERTIES COMPILE_FLAGS "${YAKL_F90_FLAGS}")

target_link_libraries(fortran_driver ${NF_LIBS})


