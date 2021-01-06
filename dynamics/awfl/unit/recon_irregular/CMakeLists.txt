
set(MYSRC recon_irregular.cpp)

add_executable(recon_irregular ${MYSRC})
target_compile_features(recon_irregular PUBLIC cxx_std_14)
include_directories(../..)
target_link_libraries(recon_irregular yakl)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../../../externals/eigen)

set_source_files_properties(${MYSRC} PROPERTIES COMPILE_FLAGS "${YAKL_CXX_FLAGS}")
if ("${ARCH}" STREQUAL "CUDA")
  set_source_files_properties(${MYSRC} PROPERTIES LANGUAGE CUDA)
  include_directories(../../cub)
endif()
if ("${ARCH}" STREQUAL "HIP")
  include_directories(../../hipCUB/hipcub/include)
  include_directories(../../rocPRIM/rocprim/include)
endif()

add_test(NAME recon_irregular_test COMMAND ./recon_irregular)

