# HyperSpharm
Hyperspharm library

[HyperSpharm Equations](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4033314/)   
[Spharm Equations](https://arxiv.org/pdf/1202.6522.pdf)

## Satus

Done:
  - Implemented small FFT library

TODO:
  - Write tests for FFT library
  - Write the hyperspharm equations (3 days)
  - Try to simplify/speed up computation by finding link to FFT (1 week)
  - Write the header files (1-2 weeks)
  - Implement the functions (2-3 weeks)
  - Write the tests (1-2weeks)

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
