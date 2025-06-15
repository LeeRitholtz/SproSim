#!/bin/bash

# SproSim Demo Runner Script
# This script builds and runs the coffee brewing simulation demo

set -e  # Exit on any error

echo "🔧 Building SproSim..."
echo "======================="

# Create build directory if it doesn't exist
if [ ! -d "build" ]; then
    echo "Creating build directory..."
    mkdir build
fi

# Navigate to build directory
cd build

# Configure and build
echo "Configuring with CMake..."
cmake .. -DCMAKE_BUILD_TYPE=Release

echo "Building project..."
cmake --build . --config Release

echo ""
echo "✅ Build complete!"
echo ""

# Run the demo
echo "☕ Running Coffee Brewing Simulation Demo..."
echo "============================================="
echo ""

./sprosim_demo

echo ""
echo "🎉 Demo completed successfully!"
echo ""
echo "💡 You can also run:"
echo "   ./sprosim_tests  - to run the test suite"
echo "   ctest            - to run tests via CTest"
