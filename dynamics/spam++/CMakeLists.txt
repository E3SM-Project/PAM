file(GLOB SRC *.cpp grids/*.cpp models/*.cpp operators/differential/*.cpp operators/reconstructions/*.cpp timesteppers/*.cpp)

add_library(spam++dycore STATIC ${SRC})
