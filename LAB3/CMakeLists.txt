cmake_minimum_required(VERSION 3.10)
project(SplineProject)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W3")
endif()

include_directories(headers)

add_executable(SplineTest
    main.cpp
    sources/Point.cpp
    sources/Smoothing_Spline_1D.cpp
)