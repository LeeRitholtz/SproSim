# ParaView Visualization Examples

## Quick Start

### 1. Run the Simulation
```bash
cd SproSim/build
./sprosim_demo --vtk
```

This creates VTK files in `./demo_output/`:
- `brewing_simulation.pvd` - Main collection file (open this)
- `timestep_particles_*.vtu` - Particle data for each timestep
- `timestep_flow_*.vts` - Flow field data for each timestep

### 2. Load in ParaView

1. Open ParaView
2. File → Open → Select `brewing_simulation.pvd`
3. Click "Apply" in the Properties panel
4. The timeline will show all exported timesteps

## Particle Visualization

1. **Create Glyphs**: Filters → Glyph → Sphere, Scale Factor ~5.0
2. **Color by Extraction**: Set Coloring to `extraction_state` (blue = fresh, brown = extracted)
3. **Animate**: Use the VCR controls to play the time series

### Available Data Arrays

| Array | File Type | Description |
|-------|-----------|-------------|
| `extraction_state` | `.vtu` | 0.0 (fresh) to 1.0 (fully extracted) |
| `concentration` | `.vtu` | Dissolved solids concentration |
| `particle_size` | `.vtu` | Particle diameter in meters |
| `velocity` | `.vts` | 3-component water velocity field |
| `pressure` | `.vts` | Water pressure field |

## Flow Field Visualization

1. Load a `timestep_flow_*.vts` file
2. Color by `velocity` magnitude or `pressure`
3. Filters → Glyph → Arrow to show velocity vectors

## State File

Load `../states/espresso_basic.pvsm` via File → Load State for a pre-configured particle visualization.

For more on ParaView: https://www.paraview.org/tutorials/