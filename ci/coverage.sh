#!/bin/bash

cd $TRAVIS_BUILD_DIR
# capture coverage info
lcov --directory . --capture --output-file coverage.info
# filter out system
lcov --remove coverage.info '/usr/*' --output-file coverage.info
# filter out tests
lcov --remove coverage.info 'tests/*' --output-file coverage.info
lcov --list coverage.info #debug info
