#!/bin/bash

cd $TRAVIS_BUILD_DIR
mkdir build
cd build
cmake .. -DTESTS:=1 -DCOVERAGE:=1
make && make tests
