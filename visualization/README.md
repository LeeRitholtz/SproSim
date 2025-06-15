# SproSim Visualization

This directory contains visualization tools and interfaces for the SproSim espresso simulation library.

## Architecture Overview

We use a **two-phase visualization approach**:

1. **Phase 1 (Current)**: ParaView/VTK for research and development
2. **Phase 2 (Future)**: Web-based interface for end users

## Directory Structure

```
visualization/
├── README.md              # This file
├── paraview/              # ParaView/VTK visualization (Phase 1)
│   ├── exporters/         # VTK export classes
│   ├── scripts/           # ParaView Python scripts
│   ├── states/            # Saved ParaView state files
│   └── examples/          # Example visualization workflows
└── web/                   # Web-based visualization (Phase 2)
    ├── backend/           # REST API server
    ├── frontend/          # JavaScript/WebGL interface
    └── shared/            # Common data formats
```

## Phase 1: ParaView/VTK Visualization

### Purpose
- **Research & Development**: Detailed analysis of simulation physics
- **Publication Quality**: High-resolution figures for papers
- **Debugging**: Visual verification of simulation correctness
- **Collaboration**: Share VTK files with other researchers

### Components
- **VTK Exporters**: C++ classes to export simulation data to VTK format
- **ParaView Scripts**: Python automation for common visualization tasks  
- **State Files**: Pre-configured ParaView sessions for different analysis types
- **Documentation**: Workflows for typical visualization tasks

### Usage
```cpp
// In your simulation code:
VTKExporter exporter;
exporter.export_particles(coffee_bed.particles(), "particles_001.vtu");
exporter.export_flow_field(water_flow, "flow_001.vtu");
exporter.export_time_series("simulation.pvd");
```

Then open the `.vtu` or `.pvd` files in ParaView for visualization.

## Phase 2: Web Interface (Future)

### Purpose
- **End User Access**: Browser-based simulation for coffee professionals
- **Real-time Interaction**: Parameter adjustment during simulation
- **Demonstrations**: Easy sharing and presentation of results
- **Education**: Interactive learning tool for espresso physics

### Planned Components
- **Backend API**: REST endpoints for simulation control and data streaming
- **3D Renderer**: WebGL-based particle and flow visualization
- **Dashboard**: Real-time plots and metrics
- **Parameter Controls**: Interactive sliders and input fields

### Planned Usage
```javascript
// Future web interface:
const simulator = new EspressoSimulator();
simulator.setParameters({
    grindSize: 0.5,
    pressure: 9.0,
    temperature: 93.0
});
simulator.start();
simulator.onUpdate((data) => {
    visualizer.updateParticles(data.particles);
    dashboard.updateMetrics(data.metrics);
});
```

## Development Roadmap

### Current Focus (v0.2.0)
- [ ] Implement VTK export functionality
- [ ] Create basic ParaView workflow
- [ ] Document visualization best practices
- [ ] Provide example ParaView state files

### Future Development (v0.4.0+)
- [ ] Design REST API for real-time data
- [ ] Implement web-based 3D visualization
- [ ] Create interactive parameter controls
- [ ] Add comparative visualization tools

## Getting Started

### For Developers (Phase 1)
1. Install ParaView from https://www.paraview.org/
2. Build SproSim with VTK export enabled
3. Run simulation with VTK export
4. Load generated files in ParaView
5. Use provided state files for common visualizations

### For End Users (Phase 2 - Coming Later)
1. Open web browser
2. Navigate to simulation interface
3. Adjust espresso parameters
4. Run simulation and view results in real-time

## Contributing

When adding new visualization features:

1. **ParaView/VTK**: Add exporters to `paraview/exporters/`
2. **Web Interface**: Follow REST API conventions in `web/backend/`
3. **Documentation**: Update this README and add examples
4. **Testing**: Include visualization tests in main test suite

## Dependencies

### Phase 1 (Current)
- VTK library (for data export)
- ParaView (for visualization)

### Phase 2 (Future)
- Web server framework (crow, httplib, etc.)
- JavaScript libraries (Three.js, Chart.js)
- WebGL-capable browser

---

*Last updated: 2024*
*See main TODO.md for detailed development priorities*