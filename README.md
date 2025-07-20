# SproSim - A "Vibe Coded" Espresso Extraction Simulation Library

SproSim is a C++ library for espresso simulation, enabling data-driven understanding of brewing parameter impacts and systematic coffee research. It models complete physics - from puck preparation through extraction - to predict shot outcomes and systematically explore recipe parameter spaces. This provides a hackable platform for testing coffee theories with quantitative data.

### *Disclaimer & Thoughts*:
This project is very experimental and well over 85% of the code (and docs) has been developed by an LLM! The code base is basically a disaster right now, but it is slowly improving.

The bottom line is that, right now, building a significant project with an AI coding agent/assistant isn't as simple as writing "build me an espresso extraction physics simulation library." I believe that it would be very hard to be productive on a real project if you don't have experience with the domain you're building in as well as the languages and tools being used. You also do need to know your code base, even if most of it is generated, because it can very quickly become an unscalable and unmaintainable mess (have it generate UML diagrams!). A clean and clear architecture approach definitely helps the LLM perform better because it provides natural constraints. It would be interesting to experiment with what combination of system design patterns, languages, and tools lend themselves best to LLM co-development. Overall, I'd say building this project with an LLM has been an awesome experience so far. It's certainly a very powerful productivity multiplier and it is definitely a skill that improves over time with increased reps.

## Features

- **Parameter Impact Analysis**: Understand how dose, pressure, grind size, time, and temperature choices affect extraction outcomes
- **Parameter Space Exploration**: Systematically test brewing parameter combinations and their interactions
- **Puck Preparation Analysis**: Model how tamping, distribution, and bed compaction choices influence flow and extraction
- **Extraction Prediction**: Forecast yield, strength, and uniformity patterns from brewing parameters
- **Quantitative Brewing Research**: Generate data-driven insights into brewing cause-and-effect relationships
- **Brewing Theory Testing**: Validate coffee science hypotheses with controlled simulation experiments
- **Complete Physics Modeling**: From particle-level extraction to flow dynamics and bed mechanics
- **Research Visualization**: 3D analysis of flow patterns, extraction gradients, and brewing dynamics

## Physics Architecture

- **Flow Models**: Pluggable flow physics (current: Darcy, future: Forchheimer, Hagen-Poiseuille)
- **Transport Models**: Swappable extraction theories (current: LinearExtraction, future: multi-compound, etc.)
- **Permeability Models**: How bed structure affects flow resistance (current: Constant, KozenyCarman, future: Ergun equation, Burke-Plummer, Happel-Brenner)
- **Particle Models**: Individual particle behavior and physics (current: 2D/3D particles, future: non-spherical shapes, non-uniform material properties, particle interactions)

## Requirements

- CMake 3.20+, C++20 compiler
- Catch2 (auto-installed for testing)
- ParaView (optional, for visualization)

## Quick Start

```bash
./run_demo.sh                    # Build and run with defaults
./run_demo.sh --scenario espresso # Run predefined espresso shot
./build/tests/sprosim_tests      # Run test suite
```

## Building

```bash
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build .
```

## Library Usage

### Basic Example

```cpp
#include "sprosim/SimulationManager.h"
using namespace sprosim;

SimulationManager::Configuration config;
config.dose = 18.0;           // grams
config.target_yield = 36.0;   // grams
config.brewing_time = 30.0;   // seconds
config.pressure_bar = 9.0;
config.temperature_celsius = 95.0;

SimulationManager sim(config);
auto result = sim.run_simulation();
```

### Configuration Options

```cpp
config.dose = 18.0;                    // Coffee dose (grams)
config.particle_count = 350;           // Number of particles
config.pressure_bar = 9.0;             // Brewing pressure
config.temperature_celsius = 95.0;     // Water temperature
config.brewing_time = 30.0;            // Max brew time (seconds)
config.target_yield = 36.0;            // Target output (grams)
config.use_3d = false;                 // Use 3D particles
config.enable_vtk_export = false;      // Export for ParaView
```

### Advanced Usage with Transport Models

```cpp
#include "sprosim/PhysicsSolver.h"
#include "sprosim/models/transport/LinearExtractionModel.h"

// Create components
auto bed = std::make_shared<CoffeeBed>(18.0, 58.0);
auto flow = std::make_shared<WaterFlow>(58, 30);
auto transport_model = std::make_shared<LinearExtractionModel>();
auto flow_model = std::make_shared<DarcyFlowModel>(/* permeability model */);

PhysicsSolver::Parameters params{/* physics parameters */};
PhysicsSolver solver(bed, flow, params, flow_model, transport_model);

// Custom transport model
class MyTransportModel : public ITransportModel {
    void update_extraction(std::shared_ptr<IWaterFlow> water_flow,
                          std::shared_ptr<CoffeeBed> bed,
                          const Parameters& params, double dt) override {
        // Your extraction physics here
    }
};
```

## Project Structure

```
include/sprosim/          # Public API
├── SimulationManager.h   # High-level simulation control
├── PhysicsSolver.h       # Low-level physics engine
├── interfaces/           # Abstract interfaces
└── models/              # Physics models (transport, flow, permeability)
src/                     # Implementation
tests/                   # Unit tests
demo/                    # Example usage
visualization/           # ParaView integration
```

## Visualization

Export VTK files for ParaView analysis:

```bash
./build/sprosim_demo --vtk
paraview *.pvd  # Load time series data
```

## Running Tests

```bash
cd build && ctest                # Using CTest
./run_tests.sh                  # Using convenience script
```

## License

GPL-3.0-or-later

## Authors

Lee Ritholtz <lee@sprosim.com> - SproSim Labs
