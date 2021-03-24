
set(MYSRC burger_prim.cpp)

add_executable(burger_prim ${MYSRC})
target_compile_features(burger_prim PUBLIC cxx_std_14)
include_directories(../..)
target_link_libraries(burger_prim yakl)

set_source_files_properties(${MYSRC} PROPERTIES COMPILE_FLAGS "${YAKL_CXX_FLAGS}")
if ("${ARCH}" STREQUAL "CUDA")
  set_source_files_properties(${MYSRC} PROPERTIES LANGUAGE CUDA)
  include_directories(../../cub)
endif()
if ("${ARCH}" STREQUAL "HIP")
  include_directories(../../hipCUB/hipcub/include)
  include_directories(../../rocPRIM/rocprim/include)
endif()

# add_test(NAME burger_prim_test COMMAND ./burger_prim)

