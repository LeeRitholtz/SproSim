# SproSim - A "Vibe Coded" Espresso Extraction Simulation Library

SproSim is a C++ library that simulates espresso extraction at the particle level. It models water flow through a coffee bed, solute transport out of individual particles, and the permeability of the puck — then reports extraction yield, TDS, and other brewing outcomes. The physics models are pluggable, so you can swap in different theories for flow, transport, and permeability to see how they change the results under different brewing conditions.

## Disclaimer:
This library is very experimental and well over 80% of the code (and docs) has been developed in collaboration with an LLM.

For my thoughts on this process, as well as general updates on the project, please subscribe to the [SproSim Substack](https://sprosim.substack.com).

## Features

- **Pluggable Physics**: Interface-driven architecture — swap flow, transport, and permeability models independently
- **Particle-Level Simulation**: Models individual coffee particles with size distribution, extraction state, and concentration tracking
- **Configurable Brewing Parameters**: Dose, pressure, temperature, grind size, brew time, and target yield
- **2D and 3D Modes**: Choose between 2D (fast) and 3D (realistic) particle simulations
- **VTK Visualization**: Export time series data for ParaView analysis of flow patterns and extraction gradients
- **CLI Demo**: Pre-built scenarios for espresso, ristretto, lungo, and more

## Physics Models

| Category | Interface | Current Implementation |
|----------|-----------|----------------------|
| Flow | `IFlowModel` | Darcy flow |
| Transport | `ITransportModel` | Linear first-order extraction kinetics |
| Permeability | `IPermeabilityModel` | Constant, Kozeny-Carman |
| Particles | `ICoffeeParticle` | 2D, 3D with size distribution |

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
