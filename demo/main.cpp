#include <cstring>
#include <iostream>
#include <memory>

#include "sprosim/SimulationManager.h"

using namespace sprosim;

void print_usage(const char* program_name) {
    std::cout << "Usage: " << program_name << " [options]\n"
              << "Options:\n"
              << "  --dose GRAMS          Coffee dose (default: 18.0)\n"
              << "  --particles N         Number of particles (default: 350)\n"
              << "  --tamping-force LBS   Tamping force in lbs (default: 20.0)\n"
              << "  --grind-range MIN MAX Grind range in microns (default: 300-600)\n"
              << "  --pressure BAR        Brewing pressure (default: 9.0)\n"
              << "  --temperature C       Water temperature (default: 95.0)\n"
              << "  --time SECONDS        Brewing time (default: 30.0)\n"
              << "  --yield GRAMS         Target output yield (default: 36.0)\n"
              << "  --extraction-rate VAL Extraction rate constant (default: 0.01)\n"
              << "  --3d                  Use 3D particles\n"
              << "  --vtk                 Export VTK files\n"
              << "  --help                Show this help\n";
}

int main(int argc, char* argv[]) {
    try {
        // Default configuration
        SimulationManager::Configuration config;

        // Parse command line arguments
        for (int i = 1; i < argc; i++) {
            if (strcmp(argv[i], "--help") == 0) {
                print_usage(argv[0]);
                return 0;
            } else if (strcmp(argv[i], "--dose") == 0 && i + 1 < argc) {
                config.dose = std::stod(argv[++i]);
            } else if (strcmp(argv[i], "--particles") == 0 && i + 1 < argc) {
                config.particle_count = std::stoi(argv[++i]);
            } else if (strcmp(argv[i], "--tamping-force") == 0 && i + 1 < argc) {
                config.tamping_force = std::stod(argv[++i]);
            } else if (strcmp(argv[i], "--grind-range") == 0 && i + 2 < argc) {
                config.grind_size_min = std::stod(argv[++i]);
                config.grind_size_max = std::stod(argv[++i]);
            } else if (strcmp(argv[i], "--pressure") == 0 && i + 1 < argc) {
                config.pressure_bar = std::stod(argv[++i]);
            } else if (strcmp(argv[i], "--temperature") == 0 && i + 1 < argc) {
                config.temperature_celsius = std::stod(argv[++i]);
            } else if (strcmp(argv[i], "--time") == 0 && i + 1 < argc) {
                config.brewing_time = std::stod(argv[++i]);
            } else if (strcmp(argv[i], "--yield") == 0 && i + 1 < argc) {
                config.target_yield = std::stod(argv[++i]);
                config.enable_yield_termination = true;
            } else if (strcmp(argv[i], "--extraction-rate") == 0 && i + 1 < argc) {
                config.extraction_rate = std::stod(argv[++i]);
            } else if (strcmp(argv[i], "--3d") == 0) {
                config.use_3d_particles = true;
            } else if (strcmp(argv[i], "--vtk") == 0) {
                config.export_vtk = true;
            }
        }

        // Print configuration
        std::cout << "=== SproSim Brewing Demo ===\n";
        std::cout << "Dose: " << config.dose << "g | ";
        std::cout << "Particles: " << config.particle_count << " | ";
        std::cout << "Grind: " << config.grind_size_min << "-" << config.grind_size_max << " μm | ";
        std::cout << "Pressure: " << config.pressure_bar << " bar\n";
        std::cout << "Temperature: " << config.temperature_celsius << "°C | ";
        std::cout << "Time: " << config.brewing_time << "s | ";
        std::cout << "Yield: " << config.target_yield << "g | ";
        std::cout << "Type: " << (config.use_3d_particles ? "3D" : "2D") << "\n\n";

        // Create simulation manager
        SimulationManager manager(config);

        // Set progress callback
        manager.set_progress_callback([](double time, double progress, double extraction) {
            std::cout << "Progress: " << std::fixed << std::setprecision(1) << progress << "% | ";
            std::cout << "Extraction: " << extraction * 100 << "%\n";
        });

        // Run simulation
        std::cout << "Running simulation...\n";
        auto results = manager.run_simulation();

        // Display results
        std::cout << "\n=== Brewing Results ===\n";
        std::cout << "Extraction Yield: " << results.final_extraction_yield * 100 << "%\n";
        std::cout << "Total Dissolved Solids: " << results.total_dissolved_solids * 100 << "%\n";
        std::cout << "Final Porosity: " << results.final_porosity << "\n";
        std::cout << "Simulation Time: " << results.simulation_time << "s\n";
        std::cout << "Termination: " << results.termination_reason << "\n";

        // Quality assessment
        std::cout << "\nQuality: ";
        if (results.is_under_extracted) {
            std::cout << "UNDER-EXTRACTED (increase time/temperature/fineness)\n";
        } else if (results.is_over_extracted) {
            std::cout << "OVER-EXTRACTED (decrease time/temperature/fineness)\n";
        } else {
            std::cout << "WELL-EXTRACTED (balanced)\n";
        }

        if (config.export_vtk) {
            std::cout << "\nVisualization: Open " << config.vtk_output_dir
                      << "/brewing.pvd in ParaView\n";
        }

        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}
