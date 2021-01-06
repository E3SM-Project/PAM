
set(MYSRC recon_regular.cpp)

add_executable(recon_regular ${MYSRC})
target_compile_features(recon_regular PUBLIC cxx_std_14)
include_directories(../..)
target_link_libraries(recon_regular yakl)

set_source_files_properties(${MYSRC} PROPERTIES COMPILE_FLAGS "${YAKL_CXX_FLAGS}")
if ("${ARCH}" STREQUAL "CUDA")
  set_source_files_properties(${MYSRC} PROPERTIES LANGUAGE CUDA)
  include_directories(../../cub)
endif()
if ("${ARCH}" STREQUAL "HIP")
  include_directories(../../hipCUB/hipcub/include)
  include_directories(../../rocPRIM/rocprim/include)
endif()

add_test(NAME recon_regular_test COMMAND ./recon_regular)

