#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <vector>

#include "sprosim/CoffeeBed.h"
#include "sprosim/CoffeeParticle.h"
#include "sprosim/CoffeeParticle3D.h"
#include "sprosim/PhysicsSolver.h"
#include "sprosim/WaterFlow.h"

#ifdef SPROSIM_HAS_VTK
#include "paraview/exporters/VTKExporter.h"
#endif

using namespace sprosim;
using namespace sprosim::particle_3d_utils;

class RealisticEspressoDemo {
  public:
    RealisticEspressoDemo() {
        std::cout << "=== SproSim Realistic Coffee Brewing Simulation ===" << std::endl;
        std::cout << "Simulating espresso extraction with realistic random particle distribution"
                  << std::endl;
        std::cout << std::endl;
    }

    void run() {
        // Create realistic coffee bed with proper random distribution
        auto coffee_bed = create_realistic_coffee_bed();

        // Create enhanced water flow field
        auto water_flow = create_enhanced_water_flow();

        // Setup physics solver with realistic espresso parameters
        auto solver = create_physics_solver(coffee_bed, water_flow);

        // Run simulation with realistic visualization
        run_realistic_simulation(coffee_bed, water_flow, solver);

        // Display results with depth analysis
        display_realistic_results(coffee_bed, water_flow);
    }

  private:
    std::shared_ptr<CoffeeBed> create_realistic_coffee_bed() {
        std::cout << "Creating realistic coffee bed with proper particle distribution..."
                  << std::endl;

        // Standard espresso parameters
        double dose_grams = 18.0;
        double portafilter_diameter_mm = 58.0;
        double bed_depth_mm = 12.0; // Typical tamped bed depth

        auto bed = std::make_shared<CoffeeBed>(dose_grams, portafilter_diameter_mm);

        // Use enhanced 3D particle distribution with normal X,Y positioning
        std::cout
            << "Generating realistic 3D particle distribution with normal spatial distribution..."
            << std::endl;

        double bed_radius = portafilter_diameter_mm / 2000.0; // Convert to meters
        double bed_depth = bed_depth_mm / 1000.0;             // Convert to meters

        int total_particles = 500; // More particles for realistic density
        std::pair<double, double> size_range = {200e-6, 800e-6}; // 200-800 microns

        // Generate 3D particles with enhanced normal distribution for X,Y coordinates
        auto particles_3d =
            generate_random_3d_distribution(total_particles, bed_radius, bed_depth, size_range);

        std::cout << "Generated " << particles_3d.size()
                  << " particles with realistic spatial distribution..." << std::endl;

        // Convert 3D particles to 2D particles for the coffee bed
        for (const auto& particle_3d : particles_3d) {
            // Get 3D position and size from the enhanced distribution
            auto [x, y, z] = particle_3d.get_position_3d();
            double particle_size = particle_3d.get_size();

            // Use depth (z) as Y coordinate for 2D visualization
            double y_2d = z;

            // Create 2D particle from 3D data (X,Y,radius)
            auto particle = std::make_shared<CoffeeParticle>(x, y_2d, particle_size / 2.0);
            bed->add_particle(particle);
        }

        std::cout << "  Dose: " << dose_grams << "g" << std::endl;
        std::cout << "  Portafilter diameter: " << portafilter_diameter_mm << "mm" << std::endl;
        std::cout << "  Bed depth: " << bed_depth_mm << "mm" << std::endl;
        std::cout << "  Total particles placed: " << bed->get_particles().size() << std::endl;
        std::cout << "  Particle size range: 200-800 microns" << std::endl;
        std::cout << "  Distribution: Log-normal size + Normal X,Y positioning" << std::endl;
        std::cout << "  Placement: Enhanced 3D->2D with realistic clustering" << std::endl;
        std::cout << "  Size stratification: Built into enhanced distribution" << std::endl;
        std::cout << "  Initial bed height: " << std::fixed << std::setprecision(2)
                  << bed->get_bed_height() * 1000 << "mm" << std::endl;
        std::cout << "  Initial porosity: " << std::fixed << std::setprecision(3)
                  << bed->get_porosity() << std::endl;
        std::cout << std::endl;

        return bed;
    }

    std::shared_ptr<WaterFlow> create_enhanced_water_flow() {
        std::cout << "Creating enhanced water flow field..." << std::endl;

        // Create enhanced grid for better resolution
        size_t width = 58; // X direction (mm) - across portafilter
        size_t depth = 15; // Depth through coffee bed (mm)

        auto flow = std::make_shared<WaterFlow>(width, depth);

        std::cout << "  Grid dimensions: " << width << "x" << depth << " cells" << std::endl;
        std::cout << "  Cell size: " << flow->get_cell_size() * 1000 << "mm" << std::endl;
        std::cout << "  Flow direction: Through coffee bed depth (Y direction)" << std::endl;
        std::cout << "  Enhanced resolution for realistic flow patterns" << std::endl;
        std::cout << std::endl;

        return flow;
    }

#include "sprosim/Parameters.h"

    std::shared_ptr<PhysicsSolver> create_physics_solver(std::shared_ptr<CoffeeBed> bed,
                                                         std::shared_ptr<WaterFlow> flow) {

        std::cout << "Setting up enhanced physics solver..." << std::endl;

        // Realistic espresso brewing parameters
        Parameters params{
            .permeability = 8e-13,            // Slightly lower for finer grind
            .fluid_viscosity = 1e-3,          // Water viscosity at 95°C [Pa·s]
            .extraction_rate = 0.08,          // Higher rate for finer particles
            .temperature = 368.15,            // 95°C in Kelvin
            .inlet_pressure = 9e5,            // 9 bar inlet pressure [Pa]
            .outlet_pressure = 1e5,           // 1 bar outlet pressure [Pa]
            .particle_drag = 0.15,            // Higher drag for realistic flow
            .flow_resistance = 0.25,          // More resistance in packed bed
            .saturation_concentration = 0.18, // Higher saturation
            .temperature_factor = 0.010       // Stronger temperature dependence
        };

        auto solver = std::make_shared<PhysicsSolver>(bed, flow, params);

        std::cout << "  Brewing temperature: " << params.temperature - 273.15 << "°C" << std::endl;
        std::cout << "  Brewing pressure: "
                  << (params.inlet_pressure - params.outlet_pressure) / 1e5 << " bar" << std::endl;
        std::cout << "  Enhanced extraction rate: " << params.extraction_rate << " s⁻¹"
                  << std::endl;
        std::cout << "  Enhanced flow resistance: " << params.flow_resistance << std::endl;
        std::cout << std::endl;

        return solver;
    }

    void run_realistic_simulation(std::shared_ptr<CoffeeBed> bed, std::shared_ptr<WaterFlow> flow,
                                  std::shared_ptr<PhysicsSolver> solver) {

        std::cout << "Running realistic brewing simulation..." << std::endl;
        std::cout << "Simulating 30 seconds of realistic espresso extraction..." << std::endl;

#ifdef SPROSIM_HAS_VTK
        std::cout << "VTK export enabled - saving realistic visualization data..." << std::endl;
        sprosim::visualization::VTKExporter vtk_exporter;
        vtk_exporter.set_output_directory("./demo_realistic_output");
        vtk_exporter.start_time_series("brewing_realistic_simulation.pvd", "timestep_realistic");
#else
        std::cout << "VTK export disabled - install VTK for realistic ParaView visualization"
                  << std::endl;
#endif
        std::cout << std::endl;

        const double total_time = 30.0; // 30 seconds
        const double dt = 0.01;         // 10ms time steps
        const int total_steps = static_cast<int>(total_time / dt);
        const int report_interval = total_steps / 15;      // Report every ~6.7%
        const int vtk_export_interval = total_steps / 120; // Export every 0.83%

        double current_time = 0.0;

        // Track extraction progress for analysis
        std::vector<double> extraction_history;
        std::vector<double> time_history;

        for (int step = 0; step < total_steps; step++) {
            solver->simulate_step(dt);
            current_time += dt;

            // Record extraction progress
            if (step % (total_steps / 100) == 0) {
                extraction_history.push_back(bed->get_extraction_yield());
                time_history.push_back(current_time);
            }

#ifdef SPROSIM_HAS_VTK
            // Export VTK data more frequently for smooth animation
            if (step % vtk_export_interval == 0 || step == total_steps - 1) {
                vtk_exporter.add_timestep_to_series(bed, flow, current_time, step);
            }
#endif

            // Report progress
            if (step % report_interval == 0 || step == total_steps - 1) {
                double progress = static_cast<double>(step) / total_steps * 100.0;
                double extraction_yield = bed->get_extraction_yield();
                double tds = bed->get_total_dissolved_solids();

                std::cout << std::fixed << std::setprecision(1);
                std::cout << "  " << std::setw(5) << progress << "% complete";
                std::cout << " | Time: " << std::setw(4) << current_time << "s";
                std::cout << " | Extraction: " << std::setprecision(2) << std::setw(5)
                          << extraction_yield << "%";
                std::cout << " | TDS: " << std::setprecision(3) << std::setw(6) << tds << "%";
                std::cout << " | Rate: " << std::setprecision(4) << std::setw(6)
                          << (extraction_yield / current_time) << "%/s" << std::endl;
            }
        }

#ifdef SPROSIM_HAS_VTK
        vtk_exporter.finalize_time_series();
        std::cout << std::endl;
        std::cout << "Realistic VTK files saved to ./demo_realistic_output/" << std::endl;
        std::cout << "Open brewing_realistic_simulation.pvd in ParaView for realistic visualization"
                  << std::endl;
#endif
        std::cout << std::endl;
    }

    void display_realistic_results(std::shared_ptr<CoffeeBed> bed,
                                   std::shared_ptr<WaterFlow> flow) {

        std::cout << "=== Realistic Brewing Results ===" << std::endl;
        std::cout << std::endl;

        // Coffee bed analysis
        std::cout << "Realistic Coffee Bed Analysis:" << std::endl;
        std::cout << "  Final extraction yield: " << std::fixed << std::setprecision(2)
                  << bed->get_extraction_yield() << "%" << std::endl;
        std::cout << "  Total dissolved solids: " << std::fixed << std::setprecision(3)
                  << bed->get_total_dissolved_solids() << "%" << std::endl;
        std::cout << "  Bed compaction: " << std::fixed << std::setprecision(3)
                  << bed->get_compaction() << std::endl;
        std::cout << "  Final porosity: " << std::fixed << std::setprecision(3)
                  << bed->get_porosity() << std::endl;
        std::cout << "  Final bed height: " << std::fixed << std::setprecision(2)
                  << bed->get_bed_height() * 1000 << "mm" << std::endl;
        std::cout << std::endl;

        // Enhanced Flow field analysis
        auto [width, height] = flow->get_grid_dimensions();
        double avg_velocity = 0.0;
        double max_velocity = 0.0;

        for (size_t i = 0; i < width; i++) {
            for (size_t j = 0; j < height; j++) {
                auto [vx, vy] = flow->get_velocity(i, j);
                double velocity_magnitude = sqrt(vx * vx + vy * vy);
                avg_velocity += velocity_magnitude;
                max_velocity = std::max(max_velocity, velocity_magnitude);
            }
        }
        avg_velocity /= (width * height);

        std::cout << "Enhanced Flow Field Analysis:" << std::endl;
        std::cout << "  Average velocity: " << std::fixed << std::setprecision(4)
                  << avg_velocity * 1000 << " mm/s" << std::endl;
        std::cout << "  Maximum velocity: " << std::fixed << std::setprecision(4)
                  << max_velocity * 1000 << " mm/s" << std::endl;
        std::cout << "  Flow uniformity: " << std::fixed << std::setprecision(3)
                  << (avg_velocity / max_velocity) << " (closer to 1.0 = more uniform)"
                  << std::endl;
        std::cout << std::endl;

        // Particle analysis by depth layers
        const auto& particles = bed->get_particles();

        // Analyze particles by depth layers (using Y coordinate as depth)
        std::vector<std::vector<double>> layer_extractions(3); // Top, middle, bottom layers
        std::vector<std::vector<double>> layer_sizes(3);       // Track sizes too

        // Find min/max Y to determine layer boundaries
        double min_y = 1e6, max_y = -1e6;
        for (const auto& particle : particles) {
            auto [x, y] = particle->get_position();
            min_y = std::min(min_y, y);
            max_y = std::max(max_y, y);
        }

        double layer_thickness = (max_y - min_y) / 3.0;

        for (const auto& particle : particles) {
            auto [x, y] = particle->get_position();
            double extraction = particle->get_extraction_state();
            double size = particle->get_size();

            // Determine layer based on Y position (depth)
            int layer;
            if (y < min_y + layer_thickness) {
                layer = 0; // Top layer (shallow)
            } else if (y < min_y + 2 * layer_thickness) {
                layer = 1; // Middle layer
            } else {
                layer = 2; // Bottom layer (deep)
            }

            layer_extractions[layer].push_back(extraction);
            layer_sizes[layer].push_back(size);
        }

        std::cout << "Particle Extraction by Depth Layers:" << std::endl;

        std::vector<std::string> layer_names = {"Top (inlet)", "Middle", "Bottom (outlet)"};
        for (size_t layer = 0; layer < 3; layer++) {
            if (!layer_extractions[layer].empty()) {
                double sum_ext = 0.0, sum_size = 0.0;
                double min_ext = 1.0, max_ext = 0.0;
                double min_size = 1e6, max_size = 0.0;

                for (size_t i = 0; i < layer_extractions[layer].size(); i++) {
                    double ext = layer_extractions[layer][i];
                    double size = layer_sizes[layer][i];

                    sum_ext += ext;
                    sum_size += size;
                    min_ext = std::min(min_ext, ext);
                    max_ext = std::max(max_ext, ext);
                    min_size = std::min(min_size, size);
                    max_size = std::max(max_size, size);
                }

                double avg_ext = sum_ext / layer_extractions[layer].size();
                double avg_size = sum_size / layer_sizes[layer].size();

                std::cout << "  " << layer_names[layer] << " layer ("
                          << layer_extractions[layer].size() << " particles):" << std::endl;
                std::cout << "    Average extraction: " << std::fixed << std::setprecision(2)
                          << avg_ext * 100 << "%" << std::endl;
                std::cout << "    Extraction range: " << std::fixed << std::setprecision(2)
                          << min_ext * 100 << "% - " << max_ext * 100 << "%" << std::endl;
                std::cout << "    Average particle size: " << std::fixed << std::setprecision(0)
                          << avg_size * 1e6 << " microns" << std::endl;
                std::cout << "    Size range: " << std::fixed << std::setprecision(0)
                          << min_size * 1e6 << " - " << max_size * 1e6 << " microns" << std::endl;
            }
        }
        std::cout << std::endl;

        // Overall extraction uniformity
        double min_extraction = 1.0, max_extraction = 0.0, avg_extraction = 0.0;
        double min_size = 1e6, max_size = 0.0, avg_size = 0.0;

        for (const auto& particle : particles) {
            double extraction = particle->get_extraction_state();
            double size = particle->get_size();

            min_extraction = std::min(min_extraction, extraction);
            max_extraction = std::max(max_extraction, extraction);
            avg_extraction += extraction;

            min_size = std::min(min_size, size);
            max_size = std::max(max_size, size);
            avg_size += size;
        }
        avg_extraction /= particles.size();
        avg_size /= particles.size();

        std::cout << "Overall Realistic Extraction Distribution:" << std::endl;
        std::cout << "  Minimum extraction: " << std::fixed << std::setprecision(2)
                  << min_extraction * 100 << "%" << std::endl;
        std::cout << "  Maximum extraction: " << std::fixed << std::setprecision(2)
                  << max_extraction * 100 << "%" << std::endl;
        std::cout << "  Average extraction: " << std::fixed << std::setprecision(2)
                  << avg_extraction * 100 << "%" << std::endl;
        std::cout << "  Extraction uniformity: " << std::fixed << std::setprecision(3)
                  << (1.0 - (max_extraction - min_extraction)) << std::endl;
        std::cout << std::endl;

        std::cout << "Realistic Particle Size Distribution:" << std::endl;
        std::cout << "  Minimum size: " << std::fixed << std::setprecision(0) << min_size * 1e6
                  << " microns" << std::endl;
        std::cout << "  Maximum size: " << std::fixed << std::setprecision(0) << max_size * 1e6
                  << " microns" << std::endl;
        std::cout << "  Average size: " << std::fixed << std::setprecision(0) << avg_size * 1e6
                  << " microns" << std::endl;
        std::cout << std::endl;

        // Brewing quality assessment
        std::cout << "Realistic Brewing Quality Assessment:" << std::endl;
        double extraction_yield = bed->get_extraction_yield();
        if (extraction_yield < 18.0) {
            std::cout << "  Result: UNDER-EXTRACTED (sour, weak)" << std::endl;
            std::cout << "  Recommendation: Grind finer or increase contact time" << std::endl;
        } else if (extraction_yield > 24.0) {
            std::cout << "  Result: OVER-EXTRACTED (bitter, harsh)" << std::endl;
            std::cout << "  Recommendation: Grind coarser or reduce contact time" << std::endl;
        } else {
            std::cout << "  Result: WELL-EXTRACTED (balanced, sweet)" << std::endl;
            std::cout << "  Recommendation: Maintain current parameters" << std::endl;
        }

        double extraction_uniformity = 1.0 - (max_extraction - min_extraction);
        if (extraction_uniformity > 0.8) {
            std::cout << "  Uniformity: EXCELLENT (even extraction throughout bed)" << std::endl;
        } else if (extraction_uniformity > 0.6) {
            std::cout << "  Uniformity: GOOD (mostly even with some variation)" << std::endl;
        } else {
            std::cout << "  Uniformity: POOR (significant channeling or uneven flow)" << std::endl;
            std::cout << "  Recommendation: Check grind distribution and tamping technique"
                      << std::endl;
        }

        std::cout << std::endl;

#ifdef SPROSIM_HAS_VTK
        std::cout << "=== Realistic ParaView Visualization Instructions ===" << std::endl;
        std::cout << "1. Install ParaView from https://www.paraview.org/" << std::endl;
        std::cout
            << "2. Open ParaView and load: ./demo_realistic_output/brewing_realistic_simulation.pvd"
            << std::endl;
        std::cout << "3. Click 'Apply' to load the realistic time series data" << std::endl;
        std::cout << "4. Create sphere glyphs to visualize particles" << std::endl;
        std::cout << "5. Color particles by 'extraction_state' to see extraction progress"
                  << std::endl;
        std::cout << "6. Color particles by 'particle_size' to see grind distribution" << std::endl;
        std::cout << "7. Use cross-sections (Filters → Slice) to see depth layers" << std::endl;
        std::cout << "8. Compare with spiral demo - notice the realistic random distribution!"
                  << std::endl;
        std::cout << "9. Play animation to watch 30 seconds of realistic brewing" << std::endl;
        std::cout << std::endl;
        std::cout << "Key Differences from Spiral Demo:" << std::endl;
        std::cout << "- RANDOM particle distribution (no artificial spiral)" << std::endl;
        std::cout << "- Log-normal size distribution (realistic grind)" << std::endl;
        std::cout << "- Size stratification (larger particles settle deeper)" << std::endl;
        std::cout << "- Y-axis represents depth through coffee bed" << std::endl;
        std::cout << "- More realistic particle clustering and packing" << std::endl;
        std::cout << std::endl;
#endif

        std::cout << "=== Realistic Simulation Complete ===" << std::endl;
        std::cout << "Realistic particle distribution with " << particles.size() << " particles"
                  << std::endl;
        std::cout << "Log-normal size distribution (200-800 microns)" << std::endl;
        std::cout << "Random spatial distribution with size stratification" << std::endl;
        std::cout << "Enhanced physics for improved brewing dynamics" << std::endl;
        std::cout << "NO ARTIFICIAL SPIRAL - truly random particle placement!" << std::endl;
    }
};

int main() {
    try {
        RealisticEspressoDemo demo;
        demo.run();
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}
