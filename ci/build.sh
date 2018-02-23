#!/bin/bash

cd $TRAVIS_BUILD_DIR
mkdir build
cd build
cmake ..
make