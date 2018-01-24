# HyperSpharm
Hyperspharm library

[HyperSpharm Equations](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4033314/)   
[Spharm Equations](https://arxiv.org/pdf/1202.6522.pdf)

## Satus

Done:
  - Implemented small FFT library
  - Write tests for FFT library
  - Write Factorial function with memoisation
  - Write Prime Factorization function with memoisation
  

TODO:
  - Write Legendre Polynomial functions with memoisation (Coefficient calculation done)
  - Write the hyperspharm equations 
  - Try to simplify/speed up computation by finding link to FFT (1 week)
  - Write the header files 
  - Implement the hypersharm functions 
  - Write tests for the hypersharms functions

## Requirements
  - cmake >= 3.0.2
  - gcc >= 5.4.0

### For testing only
  - google test: 
    ```
    sudo apt-get install libgtest-dev
    cd /usr/src/gtest
    sudo cmake
    sudo make
    sudo cp *.a /usr/lib
    ```
