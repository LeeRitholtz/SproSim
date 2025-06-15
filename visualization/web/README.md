# Web Interface (Phase 2)

This directory will contain the web-based visualization interface for SproSim, planned for a future release.

## Overview

The web interface will provide real-time, interactive visualization of espresso extraction simulations, designed for end users including coffee professionals, equipment manufacturers, and educators.

## Planned Architecture

### Backend (`backend/`)
- **REST API Server**: Real-time simulation control and data streaming
- **WebSocket Support**: Live data updates during simulation
- **Simulation Management**: Start/stop/pause simulation controls
- **Parameter Validation**: Ensure physical parameter constraints

### Frontend (`frontend/`)
- **3D Visualization**: WebGL-based particle and flow field rendering
- **Interactive Dashboard**: Real-time plots and extraction metrics
- **Parameter Controls**: Sliders and inputs for brewing parameters
- **Responsive Design**: Works on desktop, tablet, and mobile

### Shared (`shared/`)
- **Data Formats**: JSON schemas for API communication
- **Type Definitions**: Common data structures
- **Validation Rules**: Parameter bounds and constraints

## Planned Features

### Real-Time Visualization
- Live 3D particle rendering with extraction state colors
- Flow field visualization with velocity vectors and streamlines
- Interactive cross-sections through the coffee bed
- Concentration heat maps and iso-surfaces

### Interactive Controls
- **Brewing Parameters**: Grind size, dose, pressure, temperature
- **Simulation Controls**: Play/pause/speed adjustment
- **View Controls**: Camera rotation, zoom, cross-section slicing
- **Comparison Mode**: Side-by-side parameter comparisons

### Dashboard & Analytics
- **Real-time Plots**: Extraction yield, pressure, flow rate over time
- **Key Metrics**: Current yield percentage, brew ratio, extraction time
- **Predictions**: Estimated final yield and optimal stop time
- **Export Options**: Download results as CSV or images

## Technology Stack (Planned)

### Backend
- **HTTP Server**: cpp-httplib or Crow framework
- **JSON Processing**: nlohmann/json
- **WebSocket**: For real-time data streaming
- **Threading**: For non-blocking simulation execution

### Frontend
- **3D Graphics**: Three.js with WebGL
- **UI Framework**: React or Vue.js
- **Plotting**: Chart.js or D3.js
- **Styling**: CSS3 with responsive design

### Communication
- **REST API**: JSON over HTTP for configuration
- **WebSocket**: Binary or JSON for real-time data
- **CORS Support**: For development and deployment flexibility

## Development Timeline

This web interface is planned for **v0.4.0** or later, after the core physics simulation and ParaView visualization are complete.

### Prerequisites
1. Stable C++ simulation core
2. Validated physics models
3. ParaView workflow established
4. Performance optimization completed

### Development Phases
1. **API Design**: Define REST endpoints and data formats
2. **Backend Implementation**: HTTP server with simulation integration
3. **Frontend Prototype**: Basic 3D visualization and controls
4. **Dashboard Development**: Real-time plots and metrics
5. **Polish & Testing**: Responsive design and cross-browser testing

## File Structure (Planned)

```
web/
├── README.md              # This file
├── backend/
│   ├── src/
│   │   ├── api_server.cpp # HTTP/WebSocket server
│   │   ├── simulation_manager.cpp # Simulation lifecycle
│   │   └── data_serializer.cpp # JSON conversion
│   ├── include/
│   └── CMakeLists.txt
├── frontend/
│   ├── src/
│   │   ├── components/    # UI components
│   │   ├── visualization/ # 3D rendering
│   │   ├── dashboard/     # Plots and metrics
│   │   └── utils/         # Helpers
│   ├── public/
│   ├── package.json
│   └── webpack.config.js
└── shared/
    ├── schemas/           # JSON API schemas
    ├── types/             # TypeScript definitions
    └── validation/        # Parameter validation rules
```

## Getting Started (Future)

When this interface is available:

### For Users
1. Open web browser
2. Navigate to simulation URL
3. Adjust espresso parameters using sliders
4. Click "Start Brewing" to begin simulation
5. Watch real-time extraction visualization
6. Export results when complete

### For Developers
1. Build and start backend server
2. Install frontend dependencies (`npm install`)
3. Start development server (`npm run dev`)
4. Open browser to `localhost:3000`
5. Backend API available at `localhost:8080/api`

## Contributing (Future)

When ready for development:
- Follow REST API conventions
- Use TypeScript for type safety
- Write unit tests for both backend and frontend
- Follow responsive design principles
- Ensure cross-browser compatibility

---

*This is a placeholder for future development. See main TODO.md and visualization/README.md for current priorities.*

*Current focus: ParaView/VTK visualization (Phase 1)*