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
add_library(libutils SHARED src/utils.cpp)

add_executable(main src/main.cpp)
target_link_libraries(main libfft libutils)

if (tests)
  enable_testing()
  find_package(GTest REQUIRED)
  find_package(Threads REQUIRED)
  include_directories(${GTEST_INCLUDE_DIRS})
  
  file(GLOB TESTS_SRC ${PROJECT_SOURCE_DIR}/tests/*.cpp)
  add_executable(tests ${TESTS_SRC})
  target_link_libraries(tests libfft libutils ${GTEST_BOTH_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
  add_test(AllTestsInFoo tests)
endif()
