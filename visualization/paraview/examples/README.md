# ParaView Visualization Examples

This directory contains examples and workflows for visualizing SproSim espresso extraction simulations in ParaView.

## Quick Start

### 1. Run the Simulation
```bash
cd SproSim/build
./sprosim_demo
```

This creates VTK files in `./demo_output/`:
- `brewing_simulation.pvd` - Main collection file
- `timestep_particles_*.vtu` - Particle data for each timestep
- `timestep_flow_*.vts` - Flow field data for each timestep

### 2. Load in ParaView

#### Option A: Manual Loading
1. Open ParaView
2. File → Open → Select `brewing_simulation.pvd`
3. Click "Apply" in the Properties panel
4. The timeline will show 101 timesteps

#### Option B: Automated Script
1. Open ParaView
2. Tools → Python Shell
3. Navigate to the script: `cd /path/to/SproSim/visualization/paraview/scripts`
4. Run: `exec(open('visualize_espresso.py').read())`
5. Execute: `quick_espresso_viz('./demo_output/brewing_simulation.pvd')`

## Visualization Workflows

### Basic Particle Visualization

1. **Load Data**: Open `brewing_simulation.pvd`
2. **Create Glyphs**: 
   - Filters → Glyph
   - Glyph Type: Sphere
   - Scale Mode: Scalar
   - Scale Factor: 5.0 (to make particles visible)
3. **Color by Extraction**:
   - In Properties, set Coloring to "extraction_state"
   - Adjust color map: blue (unextracted) → brown (extracted)
4. **Add Color Bar**: View → Color Map Editor → Show Color Legend

### Advanced Particle Analysis

1. **Size by Particle Diameter**:
   - In Glyph filter, set Scale Mode to "scalar"
   - Set Scalars to "particle_size"
   - Adjust Scale Factor as needed

2. **Filter by Extraction Level**:
   - Filters → Threshold
   - Set Array to "extraction_state"
   - Adjust minimum/maximum values to show only partially extracted particles

3. **Slice Through Coffee Bed**:
   - Filters → Slice
   - Set Origin to middle of bed
   - Set Normal to (0, 1, 0) for horizontal slice

### Flow Field Visualization

The flow field data is in separate `.vts` files. To visualize:

1. **Load Flow Data**: Open individual `timestep_flow_*.vts` files
2. **Vector Field**:
   - Color by "velocity_magnitude"
   - Filters → Glyph → Arrow to show velocity vectors
3. **Streamlines**:
   - Filters → Stream Tracer
   - Set Vectors to "velocity"
   - Add seed points at inlet

### Animation and Time Series

1. **Time Controls**: Use the VCR controls in toolbar
2. **Animation Settings**:
   - View → Animation View
   - Set Mode to "Snap To TimeSteps"
   - Adjust frame rate for smooth playback
3. **Export Animation**:
   - File → Save Animation
   - Choose format (AVI, MP4, image sequence)

## Common Visualization Recipes

### Recipe 1: Extraction Progress Over Time
```
Goal: Show how coffee particles extract over the 30-second brewing time

Steps:
1. Load brewing_simulation.pvd
2. Create sphere glyphs (Scale Factor: 3.0)
3. Color by "extraction_state"
4. Set color map: Blue → Orange → Brown
5. Play animation to see extraction progress
6. Save as movie file
```

### Recipe 2: Particle Size Distribution Analysis
```
Goal: Analyze extraction differences between particle sizes

Steps:
1. Load data and create glyphs
2. Color by "particle_size" to see size distribution
3. Add Threshold filter:
   - Array: particle_size
   - Show only large particles (> 0.0006m)
4. Compare extraction rates of different sizes
5. Use Calculator filter to create size/extraction ratio
```

### Recipe 3: Cross-Section Analysis
```
Goal: See internal structure of coffee bed during extraction

Steps:
1. Load data with glyphs
2. Add Slice filter at bed center (Y = 0.015)
3. Color slice by extraction_state
4. Add second slice at different height
5. Compare extraction patterns at different bed depths
```

### Recipe 4: Flow Pattern Analysis
```
Goal: Visualize water flow through coffee bed

Steps:
1. Load flow field: timestep_flow_000000.vts
2. Color by "velocity_magnitude"
3. Add Arrow glyphs:
   - Glyph Type: Arrow
   - Scale by velocity_magnitude
   - Scale Factor: 0.01
4. Add Streamlines:
   - Seed Type: Point Cloud
   - Place seeds at top of domain
5. Animate through time to see flow evolution
```

## Useful ParaView Features

### Data Analysis
- **Calculator**: Create custom fields (e.g., extraction rate)
- **Plot Over Time**: Track values at specific particles
- **Histogram**: Analyze distribution of extraction values
- **Probe Location**: Get exact values at specific points

### Visualization Enhancements
- **Lighting**: Adjust for better particle visibility
- **Opacity**: Make particles semi-transparent to see internal structure  
- **Annotations**: Add text showing brewing parameters
- **Multiple Views**: Compare different timesteps side-by-side

### Export Options
- **Screenshots**: High-resolution images for publications
- **Animations**: MP4/AVI movies of extraction process
- **Data**: Export processed data back to CSV/VTK
- **Interactive**: HTML export for web sharing

## Tips and Tricks

### Performance
- Use fewer timesteps for initial exploration
- Reduce sphere resolution for faster rendering
- Use Level of Detail (LOD) for large datasets

### Visual Quality
- Adjust lighting for better particle definition
- Use gradient backgrounds for professional look
- Enable ambient occlusion for depth perception
- Set appropriate point/line sizes for your display

### Scientific Analysis
- Always include scale bars and legends
- Document color map ranges for reproducibility
- Save state files (.pvsm) for consistent visualization
- Export quantitative data alongside visualizations

## Troubleshooting

### Common Issues

**Files won't load**:
- Check file paths in .pvd file
- Ensure all .vtu/.vts files are in same directory
- Verify file permissions

**Particles too small to see**:
- Increase Scale Factor in Glyph filter
- Try different glyph types (Box, Cone)
- Adjust point size in representation properties

**Animation runs too fast/slow**:
- Adjust frame rate in Animation View
- Change timestep intervals in data export
- Use Animation Settings to control playback

**Colors don't match expected values**:
- Check data range in Information panel
- Rescale color map to data range
- Verify correct array is selected for coloring

## Example State Files

This directory contains pre-configured ParaView state files:

- `espresso_basic.pvsm` - Basic particle visualization
- `extraction_analysis.pvsm` - Advanced extraction analysis setup
- `flow_visualization.pvsm` - Flow field analysis setup

To use: File → Load State → Select .pvsm file

## Next Steps

- Explore the Python scripting examples in `../scripts/`
- Try combining particle and flow visualizations
- Create custom filters for specific analysis needs
- Export results for further analysis in other tools

---

*For more ParaView tutorials, visit: https://www.paraview.org/tutorials/*
*SproSim documentation: See main project README.md*