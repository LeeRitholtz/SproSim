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

## Examples

See the [examples/](examples/) directory for compilable usage examples:

- **[basic_simulation.cpp](examples/basic_simulation.cpp)** — Configure and run a simulation using `SimulationManager`
- **[custom_transport_model.cpp](examples/custom_transport_model.cpp)** — Implement a custom `ITransportModel` and plug it into `PhysicsSolver`

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
