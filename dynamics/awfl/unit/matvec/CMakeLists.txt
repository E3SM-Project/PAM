
set(MYSRC matvec.cpp)

add_executable(matvec ${MYSRC})
target_compile_features(matvec PUBLIC cxx_std_14)
include_directories(../..)
include_directories(../../../../include)
target_link_libraries(matvec yakl)

include(${CMAKE_CURRENT_SOURCE_DIR}/../../../../externals/YAKL/process_cxx_source_files.cmake)
process_cxx_source_files(${MYSRC})

add_test(NAME matvec_test COMMAND ./matvec)

