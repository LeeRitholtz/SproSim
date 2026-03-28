# SproSim Visualization

VTK export is **implemented and working**. The `VTKExporter` class writes VTK XML files directly (no VTK library dependency) for loading in ParaView.

## Components

- **VTKExporter** (`paraview/exporters/`): Exports particle data (`.vtu`), flow fields (`.vts`), and time series collections (`.pvd`)
- **ParaView Script** (`paraview/scripts/visualize_espresso.py`): Python automation for visualization
- **State File** (`paraview/states/espresso_basic.pvsm`): Pre-configured ParaView session
- **Examples** (`paraview/examples/README.md`): Workflow documentation

## Usage

```cpp
#include "VTKExporter.h"

VTKExporter exporter;
exporter.set_output_directory("output/");
exporter.start_time_series("brewing.pvd", "timestep");

// During simulation loop:
exporter.add_timestep_to_series(bed, flow, time, timestep);

// After simulation:
exporter.finalize_time_series();
```

Then open the `.pvd` file in ParaView to view the full time series.

## Roadmap

- [x] Implement VTK export functionality (particle, flow, time series)
- [x] Create basic ParaView workflow and example state file
- [ ] Document visualization best practices more thoroughly
- [ ] Add more example ParaView state files for different analysis types
- [ ] Web-based visualization interface (future, see `web/README.md` for notes)

## Dependencies

- ParaView (for viewing exported files — not required at build time)