# SproSim - Coffee Brewing Simulation Library

SproSim is a C++ library for simulating specialty coffee brewing processes, with a focus on espresso extraction physics. The library models coffee particle extraction, water flow dynamics, and bed compaction to provide insights into brewing quality and optimization.

## Features

- **Coffee Bed Modeling**: Realistic simulation of coffee grounds with particle size distribution
- **Water Flow Physics**: Darcy flow through porous media with pressure-driven extraction
- **Extraction Kinetics**: Temperature-dependent extraction with saturation effects
- **Particle Interactions**: Drag forces and bed compaction modeling
- **Real-time Simulation**: Time-stepped physics solver for dynamic brewing analysis
- **ParaView Visualization**: Export simulation data for 3D visualization and analysis

## Requirements

- CMake 3.20 or higher
- C++20 compatible compiler (GCC 10+, Clang 12+, or MSVC 2019+)
- Catch2 v3 (for testing) - will be automatically found or can be installed via package manager
- ParaView (optional, for visualization) - download from https://www.paraview.org/

### Installing Dependencies

**Ubuntu/Debian:**
```bash
sudo apt update
sudo apt install cmake build-essential libcatch2-dev
```

**macOS (with Homebrew):**
```bash
brew install cmake catch2
# Optional: Install ParaView
brew install --cask paraview
```

**Windows:**
- Install Visual Studio 2019+ with C++ development tools
- Install CMake from https://cmake.org/download/
- Catch2 will be handled by CMake's FetchContent if not found

## Building the Project

Clean and build the project from scratch with the following steps:

1. **Clean build directory** (optional but recommended for fresh start):
   ```bash
   rm -rf build
   ```

2. **Configure and build the project** with CMake:
   ```bash
   mkdir -p build
   cd build
   cmake .. -DCMAKE_BUILD_TYPE=Release
   cmake --build . --config Release
   ```

3. **Run tests** using one of the methods below.

You can automate these steps using the IDE tasks:
- **Clean Build** — removes the build directory
- **Build Project** — configures and builds the project in Release mode

### Prerequisites Check

Before building, ensure you have all required dependencies:

```bash
# Check CMake version (must be 3.20+)
cmake --version

# Check C++ compiler
g++ --version    # or clang++ --version

# For testing (optional but recommended)
# Ubuntu/Debian: apt list --installed | grep catch2
# macOS: brew list | grep catch2
```

### Step-by-Step Build Process

(Section removed for conciseness; see above summary.)

1. **Clone and navigate to the project:**
   ```bash
   git clone <repository-url>
   cd SproSim
   ```

2. **Create and enter build directory:**
   ```bash
   mkdir build
   cd build
   ```

3. **Configure the build with CMake:**
   ```bash
   # Debug build (default, includes debugging symbols)
   cmake ..
   
   # Or Release build (optimized for performance)
   cmake .. -DCMAKE_BUILD_TYPE=Release
   
   # Or specify custom install prefix
   cmake .. -DCMAKE_INSTALL_PREFIX=/usr/local
   ```

4. **Build the project:**
   ```bash
   # Build all targets (library, demo, tests)
   cmake --build .
   
   # Or build with multiple cores for faster compilation
   cmake --build . --parallel 4
   
   # Or build specific targets only
   cmake --build . --target sprosim          # Library only
   cmake --build . --target sprosim_demo     # Demo only
   cmake --build . --target sprosim_tests    # Tests only
   ```

5. **Install the library (optional):**
   ```bash
   # Install to system directories (may require sudo)
   cmake --build . --target install
   
   # Or install to custom location
   cmake --install . --prefix /path/to/install
   ```

### Build Verification

After building, verify your build with tests as described in the "Running Tests" section.

After building, you should see these executables in your build directory:
- `sprosim_demo` - Demo application
- `sprosim_tests` - Test suite
- `libsprosim.a` - Static library file

### Troubleshooting Build Issues

**CMake version too old:**
```bash
# Ubuntu/Debian
sudo apt update && sudo apt install cmake

# macOS
brew install cmake

# Or download from https://cmake.org/download/
```

**Missing C++20 support:**
```bash
# Ubuntu/Debian - install newer GCC
sudo apt install gcc-10 g++-10
export CC=gcc-10 CXX=g++-10

# macOS - update Xcode command line tools
xcode-select --install
```

**Catch2 not found:**
```bash
# Ubuntu/Debian
sudo apt install libcatch2-dev

# macOS
brew install catch2

# Or let CMake download it automatically (no action needed)
```

## Running the Project

### Quick Start with Convenience Script

The fastest way to build and run the demo:

```bash
cd SproSim
./run_demo.sh
```

This script will:
1. Create build directory if needed
2. Configure with CMake
3. Build the project with optimized settings
4. Run the demo application
5. Show additional options

### Running the Demo Application

#### Method 1: Direct Execution
```bash
cd SproSim/build
./sprosim_demo
```

#### Method 2: From project root
```bash
cd SproSim
./build/sprosim_demo
```

#### What the Demo Does
**What the Demo Does

The demo runs a complete espresso brewing simulation featuring:

- **Coffee Setup**: 18g dose in 58mm portafilter with 350 particles
- **Grind Distribution**: Realistic particle sizes from 200-800 microns  
- **Brewing Parameters**: 9 bar pressure, 95°C temperature
- **Simulation Duration**: 30 seconds with 10ms time steps
- **Real-time Monitoring**: Progress, extraction yield, and TDS tracking
- **VTK Export**: Automatic export of visualization data for ParaView

#### Demo Output Explained

**Initial Setup:**
```
Creating coffee bed...
  Dose: 18.00g                    # Coffee mass
  Portafilter diameter: 58mm      # Basket size
  Total particles: 350            # Simulated coffee grounds
  Initial bed height: 12.50mm     # Uncompressed bed depth
  Initial porosity: 0.420         # Void fraction (42% empty space)
```

**Real-time Progress:**
```
Running brewing simulation...
   0.0% complete | Time:  0.0s | Extraction:  0.00% | TDS:  0.000%
  10.0% complete | Time:  3.0s | Extraction:  8.45% | TDS:  1.240%
  50.0% complete | Time: 15.0s | Extraction: 18.20% | TDS:  1.890%
 100.0% complete | Time: 30.0s | Extraction: 21.30% | TDS:  2.180%
```

- **Extraction**: Percentage of extractable coffee compounds removed
- **TDS**: Total Dissolved Solids concentration in the liquid

**Final Analysis:**
```
=== Final Brewing Results ===
Coffee Bed Analysis:
  Final extraction yield: 21.30%     # Total extraction (18-24% is ideal)
  Total dissolved solids: 2.18%      # Strength (1.5-2.5% typical for espresso)
  Bed compaction: 1.090               # Compression factor
  Final porosity: 0.367               # Reduced void space from pressure

Brewing Quality Assessment:
  Result: WELL-EXTRACTED (balanced, sweet)    # Quality evaluation
  Uniformity: EXCELLENT (even extraction)     # Extraction consistency

=== ParaView Visualization Instructions ===
1. Install ParaView from https://www.paraview.org/
2. Open ParaView and load: ./demo_output/brewing_simulation.pvd
3. Click 'Apply' to load the time series data
4. Use the play button to animate the extraction process
5. Color particles by 'extraction_state' to see extraction progress
6. Add velocity vectors to visualize water flow
```

## Running Tests

Two ways to run tests after building:

### Method 1: Using CTest (Recommended)
Run all tests via CTest from the build directory for a quick pass/fail summary and test discovery:

```bash
cd build
ctest --output-on-failure
```

You can also run this via the IDE task **Run CTest Runner**.

### Method 2: Direct Test Execution
Run the Catch2 test executable directly for detailed output and debug information:

```bash
cd build
./tests/sprosim_tests
```

This task is available as **Run Catch2 Tests** in your IDE.

### Understanding Test Output
- CTest output summarizes overall test results.
- Running the Catch2 executable directly provides full test case details and assertion failures.

Ensure you run all tests after building to verify project correctness.

#### Method 1: Using CTest (Recommended)
```bash
cd SproSim/build
ctest

# For verbose output
ctest --verbose

# Run specific tests
ctest -R physics  # Run tests matching "physics"

# Run tests in parallel
ctest --parallel 4
```

#### Method 2: Direct Test Execution
```bash
cd SproSim/build
./sprosim_tests

# Run specific test cases
./sprosim_tests "[physics]"           # Run physics tests only
./sprosim_tests "PhysicsSolver*"      # Run PhysicsSolver tests
./sprosim_tests --list-tests          # List all available tests
```

#### Understanding Test Output

**Successful Test Run:**
```
Test project /path/to/SproSim/build
    Start 1: sprosim_tests
1/1 Test #1: sprosim_tests ....................   Passed    0.25 sec

100% tests passed, 0 tests failed out of 1

Total Test time (real) =   0.25 sec
```

**Individual Test Details:**
```bash
./sprosim_tests --success  # Show successful test details
```

#### Test Categories

The test suite includes:
- **Physics Tests**: Solver accuracy and flow field validation
- **Extraction Tests**: Particle extraction kinetics
- **Integration Tests**: Complete brewing simulation workflows
- **Boundary Tests**: Edge cases and error conditions

### Performance Monitoring

For performance analysis during longer simulations:

```bash
# Time the demo execution
time ./sprosim_demo

# Monitor system resources (macOS)
top -pid `pgrep sprosim_demo`

# Monitor system resources (Linux)
htop -p `pgrep sprosim_demo`
```

## Visualization with ParaView

The simulation automatically exports VTK data files that can be visualized in ParaView for detailed analysis.

### Quick Visualization Setup

1. **Run the demo** to generate visualization data:
   ```bash
   cd SproSim/build
   ./sprosim_demo
   ```

2. **Open ParaView** and load the time series:
   - File → Open → Select `demo_output/brewing_simulation.pvd`
   - Click "Apply" to load all timesteps

3. **Visualize particles** with extraction colors:
   - Filters → Glyph (to make particles visible as spheres)
   - Set Glyph Type to "Sphere" 
   - Set Scale Factor to 5.0
   - Color by "extraction_state" (blue = unextracted, brown = extracted)

4. **Animate the extraction**:
   - Use VCR controls to play through the 30-second brewing process
   - Watch particles change color as they extract

### Advanced Visualization

**Flow Field Analysis:**
- Load individual flow files: `timestep_flow_*.vts`
- Color by "velocity_magnitude" to see water flow patterns
- Add arrow glyphs to show flow direction

**Cross-Section Views:**
- Use Slice filter to cut through the coffee bed
- Analyze extraction patterns at different bed depths
- Compare extraction uniformity across the bed

**Automated Visualization:**
```python
# In ParaView Python shell:
exec(open('visualization/paraview/scripts/visualize_espresso.py').read())
quick_espresso_viz('./demo_output/brewing_simulation.pvd')
```

**Export Options:**
- Save screenshots of key moments
- Export animations as MP4/AVI movies
- Generate publication-quality figures

### Visualization Data Structure

The demo creates these files in `demo_output/`:
- `brewing_simulation.pvd` - Main collection file (load this)
- `timestep_particles_*.vtu` - Particle data for each timestep
- `timestep_flow_*.vts` - Flow field data for each timestep

Each particle contains:
- Position (x, y coordinates)
- `extraction_state` (0.0 = fresh, 1.0 = fully extracted)
- `concentration` (local dissolved coffee concentration)
- `particle_size` (diameter in meters)

## Library Usage

### Basic Example

```cpp
#include "sprosim/CoffeeBed.h"
#include "sprosim/CoffeeParticle.h"
#include "sprosim/WaterFlow.h"
#include "sprosim/PhysicsSolver.h"

using namespace sprosim;

// Create coffee bed (18g dose, 58mm portafilter)
auto bed = std::make_shared<CoffeeBed>(18.0, 58.0);

// Add coffee particles
auto particle = std::make_shared<CoffeeParticle>(0.005, 0.01, 0.0005);
bed->add_particle(particle);

// Create water flow field (58x30 grid)
auto flow = std::make_shared<WaterFlow>(58, 30);

// Setup physics solver
PhysicsSolver::Parameters params{
    .permeability = 1e-12,
    .fluid_viscosity = 1e-3,
    .extraction_rate = 0.05,
    .temperature = 368.15,        // 95°C
    .inlet_pressure = 9e5,        // 9 bar
    .outlet_pressure = 1e5,       // 1 bar
    .particle_drag = 0.1,
    .flow_resistance = 0.2,
    .saturation_concentration = 0.15,
    .temperature_factor = 0.008
};

PhysicsSolver solver(bed, flow, params);

// Run simulation
for (int i = 0; i < 3000; i++) {
    solver.simulate_step(0.01);  // 10ms time steps
}

// Get results
double extraction_yield = bed->get_extraction_yield();
double tds = bed->get_total_dissolved_solids();
```

### Key Classes

- **`CoffeeBed`**: Manages collection of coffee particles and bed properties
- **`CoffeeParticle`**: Individual particle with extraction state and physics
- **`WaterFlow`**: 2D flow field with velocity and concentration tracking
- **`PhysicsSolver`**: Main simulation engine that couples all physics

## Project Structure

```
SproSim/
├── include/sprosim/          # Public headers
│   ├── interfaces/           # Abstract interfaces
│   ├── CoffeeBed.h
│   ├── CoffeeParticle.h
│   ├── WaterFlow.h
│   └── PhysicsSolver.h
├── src/                      # Implementation files
├── tests/                    # Unit tests with Catch2
├── demo/                     # Demo application
├── docs/                     # Documentation
└── CMakeLists.txt
```

## Physics Model

SproSim implements a coupled physics model including:

1. **Darcy Flow**: Pressure-driven flow through porous coffee bed
2. **Mass Transfer**: Concentration-dependent extraction kinetics
3. **Bed Mechanics**: Particle packing and compaction under pressure
4. **Thermal Effects**: Temperature-dependent extraction rates

The model is suitable for:
- Espresso brewing optimization
- Pour-over coffee analysis
- Grind size impact studies
- Equipment design validation

## Development

### Development Workflow

1. **Set up development environment:**
   ```bash
   # Build in debug mode for development
   cd SproSim/build
   cmake .. -DCMAKE_BUILD_TYPE=Debug
   cmake --build .
   ```

2. **Make incremental builds during development:**
   ```bash
   # Only rebuild changed files
   cmake --build .
   
   # Or use make directly (if using Unix Makefiles generator)
   make -j4
   ```

3. **Run tests frequently:**
   ```bash
   # Quick test run
   ctest
   
   # Run tests with debugging info
   ctest --verbose --output-on-failure
   ```

### Adding New Features

1. **Implement new classes** following the existing interface patterns in `include/sprosim/interfaces/`
2. **Add source files** to `src/` directory and update `CMakeLists.txt`:
   ```cmake
   add_library(sprosim
       src/CoffeeBed.cpp
       src/CoffeeParticle.cpp
       src/WaterFlow.cpp
       src/PhysicsSolver.cpp
       src/YourNewClass.cpp  # Add here
   )
   ```
3. **Write comprehensive unit tests** in `tests/` directory
4. **Update the demo** in `demo/main.cpp` to showcase new functionality
5. **Test thoroughly:**
   ```bash
   cmake --build . && ctest && ./sprosim_demo
   ```

### Debugging

#### Debug Builds
```bash
# Configure for debugging
cmake .. -DCMAKE_BUILD_TYPE=Debug -DCMAKE_CXX_FLAGS="-g -O0"

# Build with debug symbols
cmake --build .
```

#### Using Debuggers
```bash
# GDB (Linux/macOS)
gdb ./sprosim_demo
gdb ./sprosim_tests

# LLDB (macOS)
lldb ./sprosim_demo
lldb ./sprosim_tests
```

#### Memory Debugging
```bash
# Valgrind (Linux)
valgrind --leak-check=full ./sprosim_demo
valgrind --leak-check=full ./sprosim_tests

# AddressSanitizer (GCC/Clang)
cmake .. -DCMAKE_CXX_FLAGS="-fsanitize=address -g"
cmake --build .
./sprosim_demo  # Will automatically detect memory issues
```

### Code Quality

#### Static Analysis
```bash
# Clang Static Analyzer
scan-build cmake --build .

# Cppcheck
cppcheck --enable=all --std=c++20 src/ include/
```

#### Code Formatting

The project uses clang-format with LLVM-based style for consistent code formatting. A `.clang-format` configuration file is provided in the project root.

**Install clang-format:**
```bash
# Ubuntu/Debian
sudo apt install clang-format

# macOS
brew install clang-format

# Windows (with Visual Studio)
# clang-format is included with Visual Studio 2019+
```

**Format code:**
```bash
# Format all source files in-place
find src include -name "*.cpp" -o -name "*.h" | xargs clang-format -i

# Check formatting without modifying files
find src include -name "*.cpp" -o -name "*.h" | xargs clang-format --dry-run --Werror
```

### Contributing

1. **Fork the repository** and clone your fork
2. **Create a feature branch:**
   ```bash
   git checkout -b feature/your-feature-name
   ```
3. **Make changes with appropriate tests:**
   ```bash
   # Edit code, add tests
   cmake --build . && ctest  # Verify everything works
   ```
4. **Commit with descriptive messages:**
   ```bash
   git add .
   git commit -m "Add feature: detailed description of changes"
   ```
5. **Push and submit a pull request:**
   ```bash
   git push origin feature/your-feature-name
   ```

### Build System Customization

#### Custom CMake Options
```bash
# Disable tests
cmake .. -DBUILD_TESTING=OFF

# Enable additional warnings
cmake .. -DCMAKE_CXX_FLAGS="-Wall -Wextra -Wpedantic"

# Use different compiler
cmake .. -DCMAKE_CXX_COMPILER=clang++

# Set custom installation directory
cmake .. -DCMAKE_INSTALL_PREFIX=/opt/sprosim
```

#### Integration with IDEs

**Visual Studio Code:**
1. Install CMake Tools extension
2. Open project folder
3. Select kit (compiler)
4. Configure and build through command palette

**CLion:**
1. Open project (CLion auto-detects CMake)
2. Configure CMake settings if needed
3. Build and run through IDE

**Qt Creator:**
1. Open CMakeLists.txt as project
2. Configure build settings
3. Build and debug through IDE

## License

GPL-3.0-or-later - See LICENSE file for details.

## Authors

- Lee Ritholtz <lee@sprosim.com> - SproSim Labs

## Support

For questions or issues, please create an issue on the project repository.