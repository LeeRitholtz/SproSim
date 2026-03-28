// Basic SproSim example: configure and run an espresso simulation.
//
// Build: cmake --build build
// Run:   ./build/basic_simulation

#include <iomanip>
#include <iostream>

#include "sprosim/SimulationManager.h"

int main() {
    using namespace sprosim;

    // Configure a standard espresso shot
    SimulationManager::Configuration config;
    config.dose = 18.0;                // grams of coffee
    config.target_yield = 36.0;        // grams of liquid output
    config.brewing_time = 30.0;        // seconds
    config.pressure_bar = 9.0;         // bar
    config.temperature_celsius = 93.0; // degrees C
    config.particle_count = 350;
    config.grind_size_min = 300.0;     // microns
    config.grind_size_max = 600.0;     // microns
    config.extraction_rate = 0.01;

    // Optional: enable VTK export for ParaView visualization
    // config.export_vtk = true;
    // config.vtk_output_dir = "./vtk_output";

    // Optional: use 3D particles (slower, more realistic)
    // config.use_3d_particles = true;

    SimulationManager sim(config);

    // Optional: report progress during simulation
    sim.set_progress_callback(
        [](double time, double progress, double yield) {
            std::cout << std::fixed << std::setprecision(1)
                      << "t=" << time << "s  "
                      << progress << "%  "
                      << "yield=" << std::setprecision(3) << yield * 100.0 << "%\n";
        });

    auto result = sim.run_simulation();

    // Print results
    std::cout << "\n--- Results ---\n"
              << std::fixed << std::setprecision(2)
              << "Extraction yield: " << result.final_extraction_yield * 100.0 << "%\n"
              << "TDS:              " << result.total_dissolved_solids << "%\n"
              << "Porosity:         " << result.final_porosity << "\n"
              << "Sim time:         " << result.simulation_time << "s\n"
              << "Steps:            " << result.total_steps << "\n"
              << "Stopped because:  " << result.termination_reason << "\n";

    return 0;
}