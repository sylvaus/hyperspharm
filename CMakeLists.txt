cmake_minimum_required (VERSION 2.8)
project (HyperSpharm)

# project version
set(VERSION_MAJOR 0)
set(VERSION_MINOR 1)

if(UNIX)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Werror -std=gnu++11")
endif()

if(COVERAGE)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
endif()

FIND_PACKAGE(OpenMP REQUIRED)
if(OPENMP_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
else()
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unknown-pragmas")
endif()

# Add the header folder
include_directories(include)

add_library(libfft SHARED src/fft.cpp include/fft.h include/types.h)
add_library(libutils SHARED src/utils.cpp include/utils.h include/types.h)
add_library(liblegendre SHARED src/legendre.cpp include/legendre.h)
target_link_libraries(liblegendre libutils)

add_executable(main src/main.cpp)
target_link_libraries(main libfft libutils)


if (TESTS)
  enable_testing()
  find_package(GTest REQUIRED)
  find_package(Threads REQUIRED)
  find_package(GSL REQUIRED) 
  include_directories(${GTEST_INCLUDE_DIRS})

  add_executable(benchmark src/benchmark.cpp)
  target_link_libraries(benchmark liblegendre ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY})
  
  file(GLOB TESTS_SRC ${PROJECT_SOURCE_DIR}/tests/*.cpp)
  add_executable(tests ${TESTS_SRC})
  target_link_libraries(tests libfft libutils liblegendre ${GTEST_BOTH_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY})
  add_test(AllTests tests)
endif()
