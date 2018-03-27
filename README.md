# HyperSpharm
Hyperspharm library

[HyperSpharm Equations](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4033314/)   
[Spharm Equations](https://arxiv.org/pdf/1202.6522.pdf)

## Satus
[![Build Status](https://travis-ci.org/sylvaus/hyperspharm.svg?branch=master)](https://travis-ci.org/sylvaus/hyperspharm)
[![codecov](https://codecov.io/gh/sylvaus/hyperspharm/branch/master/graph/badge.svg)](https://codecov.io/gh/sylvaus/hyperspharm)

Done:
  - Implemented small FFT library
  - Implement tests for FFT library
  - Implement Factorial function with memoisation
  - Implement Prime Factorization function with memoisation
  - Implement Legendre Polynomial functions

TODO:
  - Implement Spharm functions (and tests)
  - Write the hyperspharm equations 
  - Try to simplify/speed up computation by finding link to FFT
  - Implement the hypersharm functions 
  - Implement tests for the hypersharms functions

## Requirements
  - cmake >= 3.0.2
  - gcc >= 5.4.0

### For testing only
  - google test: 
    ```
    sudo apt-get install libgtest-dev
    cd /usr/src/gtest
    sudo cmake .
    sudo make
    sudo cp *.a /usr/lib
    ```
  - gsl library (2.1):
    ```
    sudo apt-get install libgsl2 gsl-*
    ```
    or source install [source code](ftp://ftp.gnu.org/gnu/gsl/gsl-2.1.tar.gz)
