# Visualization CMakeLists.txt
cmake_minimum_required(VERSION 3.20)

# Build visualization components (no VTK libraries required - just XML export)
add_library(sprosim_visualization
    paraview/exporters/VTKExporter.cpp
)

# Set include directories
target_include_directories(sprosim_visualization
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
)

# Link with core library only
target_link_libraries(sprosim_visualization
    PUBLIC
        sprosim
)

# Install visualization library
install(TARGETS sprosim_visualization
    EXPORT sprosim-targets
    INCLUDES DESTINATION include
)

# Install visualization headers
install(DIRECTORY paraview/exporters/
    DESTINATION include/sprosim/visualization
    FILES_MATCHING PATTERN "*.h"
)

message(STATUS "Building SproSim visualization (VTK XML export)")

# Always install ParaView scripts and examples (they don't need compilation)
install(DIRECTORY paraview/scripts/
    DESTINATION share/sprosim/paraview/scripts
    FILES_MATCHING PATTERN "*.py"
)

install(DIRECTORY paraview/states/
    DESTINATION share/sprosim/paraview/states
    FILES_MATCHING PATTERN "*.pvsm"
)

install(DIRECTORY paraview/examples/
    DESTINATION share/sprosim/paraview/examples
)

# Web interface placeholder (future development)
# add_subdirectory(web) # Uncomment when web interface is ready
