
set(MYSRC anelastic_direct.cpp)

add_executable(anelastic_direct ${MYSRC})
target_compile_features(anelastic_direct PUBLIC cxx_std_14)
include_directories(../..)
target_link_libraries(anelastic_direct yakl)

set_source_files_properties(${MYSRC} PROPERTIES COMPILE_FLAGS "${YAKL_CXX_FLAGS}")
if ("${ARCH}" STREQUAL "CUDA")
  set_source_files_properties(${MYSRC} PROPERTIES LANGUAGE CUDA)
  include_directories(../../cub)
endif()
if ("${ARCH}" STREQUAL "HIP")
  include_directories(../../hipCUB/hipcub/include)
  include_directories(../../rocPRIM/rocprim/include)
endif()

add_test(NAME anelastic_direct_test COMMAND ./anelastic_direct)

