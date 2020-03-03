file(GLOB SRC *.cpp)

add_library(spam++dycore STATIC ${SRC})

add_executable(topology_test test/topology_test.cpp)
