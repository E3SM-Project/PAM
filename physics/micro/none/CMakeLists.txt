
add_library(micro INTERFACE)
target_compile_features(micro INTERFACE cxx_std_14)
target_include_directories(micro INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})

