#!/bin/bash

# SproSim Demo Runner Script
# This script builds and runs the coffee brewing simulation demo

set -e  # Exit on any error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Function to print colored output
print_color() {
    color=$1
    shift
    echo -e "${color}$@${NC}"
}

# Function to show usage
show_usage() {
    print_color $BLUE "SproSim Demo Runner"
    echo "==================="
    echo ""
    echo "Usage: $0 [options] -- [demo arguments]"
    echo ""
    echo "Runner Options:"
    echo "  --help, -h        Show this help message"
    echo "  --skip-build      Skip building (use existing binary)"
    echo "  --clean           Clean build before compiling"
    echo "  --debug           Build in debug mode"
    echo "  --scenario NAME   Run a predefined scenario (see below)"
    echo ""
    echo "Demo Arguments (after --):"
    echo "  --dose GRAMS          Coffee dose (default: 18.0)"
    echo "  --particles N         Number of particles (default: 350)"
    echo "  --tamping-force LBS   Tamping force in lbs (default: 20.0)"
    echo "  --grind TYPE          Grind size: coarse|medium|fine (default: medium)"
    echo "  --grind-range MIN MAX Custom grind range in microns"
    echo "  --pressure BAR        Brewing pressure (default: 9.0)"
    echo "  --temperature C       Water temperature (default: 95.0)"
    echo "  --time SECONDS        Brewing time (default: 30.0)"
    echo "  --yield GRAMS         Target output yield (default: 36.0)"
    echo "  --3d                  Use 3D particles"
    echo "  --vtk                 Export VTK files"
    echo "  --help                Show demo help"
    echo ""
    echo "Predefined Scenarios:"
    echo "  espresso       Traditional espresso (18g in, 36g out, 30s, 9 bar)"
    echo "  ristretto      Short shot (18g in, 20g out, 20s, 9 bar, fine grind)"
    echo "  lungo          Long shot (18g in, 54g out, 45s, 9 bar, coarse grind)"
    echo "  turbo          Turbo shot (20g in, 40g out, 15s, 6 bar, coarse grind)"
    echo "  filter         Filter-style (15g in, 250g out, 180s, 1 bar, coarse grind)"
    echo "  3d-demo        3D particle demonstration with VTK export"
    echo "  performance    Performance test (1000 particles)"
    echo ""
    echo "Examples:"
    echo "  $0                                    # Run with defaults"
    echo "  $0 -- --dose 20 --grind fine         # Custom parameters"
    echo "  $0 --scenario espresso                # Run espresso scenario"
    echo "  $0 --skip-build -- --3d --vtk         # Skip build, use 3D with VTK"
}

# Default values
SKIP_BUILD=false
CLEAN_BUILD=false
BUILD_TYPE="Release"
SCENARIO=""
DEMO_ARGS=""

# Parse runner options
while [[ $# -gt 0 ]]; do
    case $1 in
        --help|-h)
            show_usage
            exit 0
            ;;
        --skip-build)
            SKIP_BUILD=true
            shift
            ;;
        --clean)
            CLEAN_BUILD=true
            shift
            ;;
        --debug)
            BUILD_TYPE="Debug"
            shift
            ;;
        --scenario)
            SCENARIO="$2"
            shift 2
            ;;
        --)
            shift
            DEMO_ARGS="$@"
            break
            ;;
        *)
            DEMO_ARGS="$@"
            break
            ;;
    esac
done

# Set demo arguments based on scenario
case $SCENARIO in
    espresso)
        print_color $YELLOW "Running Espresso scenario..."
        DEMO_ARGS="--dose 18 --yield 36 --time 30 --pressure 9 --grind medium --tamping-force 30"
        ;;
    ristretto)
        print_color $YELLOW "Running Ristretto scenario..."
        DEMO_ARGS="--dose 18 --yield 20 --time 20 --pressure 9 --grind fine --tamping-force 35"
        ;;
    lungo)
        print_color $YELLOW "Running Lungo scenario..."
        DEMO_ARGS="--dose 18 --yield 54 --time 45 --pressure 9 --grind coarse --tamping-force 25"
        ;;
    turbo)
        print_color $YELLOW "Running Turbo Shot scenario..."
        DEMO_ARGS="--dose 20 --yield 40 --time 15 --pressure 6 --grind coarse --tamping-force 20"
        ;;
    filter)
        print_color $YELLOW "Running Filter-style scenario..."
        DEMO_ARGS="--dose 15 --yield 250 --time 180 --pressure 1 --grind coarse --tamping-force 5"
        ;;
    3d-demo)
        print_color $YELLOW "Running 3D Demonstration..."
        DEMO_ARGS="--3d --vtk --particles 500 --dose 18"
        ;;
    performance)
        print_color $YELLOW "Running Performance Test..."
        DEMO_ARGS="--particles 1000 --dose 20 --time 60"
        ;;
    "")
        # No scenario, use provided args or defaults
        ;;
    *)
        print_color $RED "Unknown scenario: $SCENARIO"
        exit 1
        ;;
esac

# Build if not skipping
if [ "$SKIP_BUILD" = false ]; then
    print_color $BLUE "🔧 Building SproSim..."
    echo "======================="

    # Clean if requested
    if [ "$CLEAN_BUILD" = true ] && [ -d "build" ]; then
        print_color $YELLOW "Cleaning previous build..."
        rm -rf build
    fi

    # Create build directory if it doesn't exist
    if [ ! -d "build" ]; then
        echo "Creating build directory..."
        mkdir build
    fi

    # Navigate to build directory
    cd build

    # Configure and build
    echo "Configuring with CMake ($BUILD_TYPE mode)..."
    cmake .. -DCMAKE_BUILD_TYPE=$BUILD_TYPE

    echo "Building project..."
    cmake --build . --config $BUILD_TYPE

    echo ""
    print_color $GREEN "✅ Build complete!"
    echo ""

    # Return to root directory
    cd ..
else
    print_color $YELLOW "⏭️  Skipping build (using existing binary)"
    echo ""
fi

# Check if demo binary exists
if [ ! -f "build/sprosim_demo" ]; then
    print_color $RED "Error: Demo binary not found. Please build first."
    exit 1
fi

# Run the demo
print_color $BLUE "☕ Running Coffee Brewing Simulation Demo..."
echo "============================================="

if [ -n "$DEMO_ARGS" ]; then
    print_color $YELLOW "Parameters: $DEMO_ARGS"
fi
echo ""

# Execute demo with arguments
./build/sprosim_demo $DEMO_ARGS

echo ""
print_color $GREEN "🎉 Demo completed successfully!"
echo ""

# Show additional options
print_color $BLUE "💡 Additional options:"
echo "   ./build/sprosim_tests  - Run the test suite"
echo "   ./run_tests.sh         - Run tests with coverage"
echo "   ./run_demo.sh --help   - See all demo options"

# If VTK files were created, notify user
if [[ "$DEMO_ARGS" == *"--vtk"* ]]; then
    echo ""
    print_color $YELLOW "📊 VTK files exported to current directory"
    echo "   Open *.vtp files in ParaView to visualize the simulation"
fi
