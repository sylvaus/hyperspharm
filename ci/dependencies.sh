#!/bin/bash

sudo apt-get update
sudo apt-get install -y libgtest-dev cmake make lcov
cd /usr/src/gtest
sudo cmake .
sudo make
sudo cp ./*.a /usr/lib

# For hypershparm testing
cd $TRAVIS_BUILD_DIR
wget ftp://ftp.gnu.org/gnu/gsl/gsl-2.1.tar.gz
tar -xvf gsl-2.1.tar.gz
cd gsl-2.1/
./configure
make
sudo make install