PROJECT(RayTrace)
cmake_minimum_required(VERSION 3.1)


set(CMAKE_CXX_STANDARD 11) # C++11...
set(CMAKE_CXX_STANDARD_REQUIRED ON) #...is required...
set(CMAKE_CXX_EXTENSIONS OFF) #...without compiler extensions like gnu++11
MESSAGE(STATUS "C++11 enabled. ")


#FIND_PACKAGE(wxWidgets REQUIRED)
#include(${wxWidgets_USE_FILE})


SET(SRCS
  Vector.cpp
  test.cpp
  )

add_executable(rayTracer ${SRCS})
#target_link_libraries(rayTracer ${wxWidgets_LIBRARIES})
