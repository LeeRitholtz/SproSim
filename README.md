# SproSim - A "Vibe Coded" Espresso Extraction Simulation Library

SproSim is a C++ library for espresso simulation, enabling data-driven understanding of brewing parameter impacts and systematic coffee research. It models complete physics - from puck preparation through extraction - to predict shot outcomes and systematically explore both recipe parameter spaces and physics regimes. This provides a hackable platform for testing coffee theories with quantitative data.

## Features

- **Parameter Impact Analysis**: Understand how dose, pressure, grind size, time, and temperature choices affect extraction outcomes
- **Parameter Space Exploration**: Systematically test brewing parameter combinations AND physics model choices to understand both recipe effects and model accuracy
- **Puck Preparation Analysis**: Model how tamping, distribution, and bed compaction choices influence flow and extraction
- **Extraction Prediction**: Forecast yield, strength, and uniformity patterns from brewing parameters
- **Quantitative Brewing Research**: Generate data-driven insights into brewing cause-and-effect relationships
- **Brewing Theory Testing**: Validate coffee science hypotheses with controlled simulation experiments
- **Physics Model Validation**: Map how different brewing parameter ranges activate different physics regimes (flow, transport, permeability) to understand which models apply under which conditions
- **Complete Physics Modeling**: From particle-level extraction to flow dynamics and bed mechanics
- **Research Visualization**: 3D analysis of flow patterns, extraction gradients, and brewing dynamics

## Physics Architecture

- **Flow Models**: Pluggable flow physics (Darcy)
- **Transport Models**: Swappable extraction theories (LinearExtraction)
- **Permeability Models**: How bed structure affects flow resistance (Constant, KozenyCarman)
- **Particle Models**: Individual particle behavior and physics (2D and 3D particles)

See [PROJECT_STRUCTURE.md](PROJECT_STRUCTURE.md) for the full directory layout and architecture diagram.

## Requirements

- CMake 3.20+, C++20 compiler
- Catch2 (auto-installed for testing)
- ParaView (optional, for visualization)

## Quick Start

```bash
# Build
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build .

# Run demo
./run_demo.sh                    # Build and run with defaults
./run_demo.sh --scenario espresso # Run predefined espresso shot

# Run tests
./run_tests.sh                   # Using convenience script
cd build && ctest                # Or using CTest directly
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

## Visualization

Export VTK files for ParaView analysis:

```bash
./build/sprosim_demo --vtk
paraview *.pvd  # Load time series data
```

See [visualization/README.md](visualization/README.md) for details.

## License

GPL-3.0-or-later

## Authors

Lee Ritholtz <lee@sprosim.com> - SproSim Labs
