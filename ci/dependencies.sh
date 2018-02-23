#!/bin/bash

sudo apt-get update
sudo apt-get install -y libgtest-dev cmake make lcov
cd /usr/src/gtest
sudo cmake .
sudo make
sudo cp ./*.a /usr/lib

# For hypershparm librarygcov main_test.cpp
cd $TRAVIS_BUILD_DIR
sudo apt-get install -y libgsl* gsl-*