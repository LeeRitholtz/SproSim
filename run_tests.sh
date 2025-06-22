#!/bin/bash

# Script to build the project and run tests from the build directory

set -e

echo "Building project..."
cmake --build build --parallel

echo "Running tests..."
(cd build && ctest --output-on-failure)
