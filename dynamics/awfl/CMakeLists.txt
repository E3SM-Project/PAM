
add_library(dycore INTERFACE)
target_compile_features(dycore INTERFACE cxx_std_14)
target_include_directories(dycore INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})

