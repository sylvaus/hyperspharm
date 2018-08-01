cmake_minimum_required (VERSION 3.2)
project (HyperSpharm)

# project version
set(VERSION_MAJOR 0)
set(VERSION_MINOR 1)

set (CMAKE_CXX_STANDARD 11)

if(UNIX)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Werror")
endif()

if(MAX_OPTIMIZATION)
    MESSAGE("Max optimization on. Setting flags \"-ffast-math -march=native -flto\"")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ffast-math -march=native -flto")
endif()

if(COVERAGE)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
endif()

FIND_PACKAGE(OpenMP REQUIRED)
if(OPENMP_FOUND)
    MESSAGE("OpenMP found. Setting the openmp flags \"${OpenMP_CXX_FLAGS}\"")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
else()
    MESSAGE("No OpenMP found")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unknown-pragmas")
endif()

# Add the header folder
include_directories(include)

add_library(libfft STATIC src/fft.cpp include/fft.h include/types.h)
add_library(libutils STATIC src/utils.cpp include/utils.h include/types.h)
add_library(liblegendre STATIC src/legendre.cpp include/legendre.h)
target_link_libraries(liblegendre libutils)
add_library(libgegenbauer STATIC src/gegenbauer.cpp include/gegenbauer.h)
target_link_libraries(libgegenbauer libutils)

add_library(libspharm STATIC src/spharms.cpp include/spharms.h)
target_link_libraries(libspharm liblegendre libutils)

add_library(libhyperspharm STATIC src/hyperspharm.cpp include/hyperspharm.h)
target_link_libraries(libhyperspharm libgegenbauer liblegendre libutils)

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
    target_link_libraries(tests libspharm libfft libutils liblegendre libgegenbauer ${GTEST_BOTH_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY})
    add_test(AllTests tests)
endif()
