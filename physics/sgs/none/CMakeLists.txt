
add_library(sgs INTERFACE)
target_compile_features(sgs INTERFACE cxx_std_14)
target_include_directories(sgs INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})

