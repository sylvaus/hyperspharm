cmake_minimum_required (VERSION 2.8)
project (HyperSpharm)

# project version
set(VERSION_MAJOR 0)
set(VERSION_MINOR 1)

if(UNIX)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Werror -std=gnu++11")
endif()

# Add the header folder
include_directories(include)

add_library(libfft SHARED src/fft.cpp)

add_executable(main src/main.cpp)
target_link_libraries(main libfft)
