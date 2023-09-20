
add_library(dycore INTERFACE)
target_include_directories(dycore INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_definitions(dycore INTERFACE PAMA_DYCORE)

