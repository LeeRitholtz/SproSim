#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <random>
#include <vector>

#include "sprosim/CoffeeBed.h"
#include "sprosim/CoffeeParticle.h"
#include "sprosim/PhysicsSolver.h"
#include "sprosim/WaterFlow.h"

#ifdef SPROSIM_HAS_VTK
#include "paraview/exporters/VTKExporter.h"
#endif

using namespace sprosim;

class Realistic3DEspressoDemo {
  public:
    Realistic3DEspressoDemo() {
        std::cout << "=== SproSim 3D Realistic Coffee Brewing Simulation ===" << std::endl;
        std::cout << "Simulating espresso extraction with realistic 3D particle distribution"
                  << std::endl;
        std::cout << std::endl;
    }

    void run() {
        // Create realistic 3D coffee bed
        auto coffee_bed = create_realistic_3d_coffee_bed();

        // Create 3D water flow field
        auto water_flow = create_3d_water_flow();

        // Setup physics solver with realistic espresso parameters
        auto solver = create_physics_solver(coffee_bed, water_flow);

        // Run simulation with 3D visualization
        run_3d_simulation(coffee_bed, water_flow, solver);

        // Display results with 3D analysis
        display_3d_results(coffee_bed, water_flow);
    }

  private:
    std::shared_ptr<CoffeeBed> create_realistic_3d_coffee_bed() {
        std::cout << "Creating realistic 3D coffee bed..." << std::endl;

        // Standard espresso parameters
        double dose_grams = 18.0;
        double portafilter_diameter_mm = 58.0;
        double bed_depth_mm = 12.0; // Typical tamped bed depth

        auto bed = std::make_shared<CoffeeBed>(dose_grams, portafilter_diameter_mm);

        // Realistic particle size distribution (log-normal)
        std::cout << "Generating realistic particle size distribution..." << std::endl;

        // Setup random number generation
        std::random_device rd;
        std::mt19937 gen(42); // Fixed seed for reproducibility

        // Log-normal distribution for particle sizes (200-800 microns, peak at ~400)
        std::lognormal_distribution<double> size_dist(std::log(400e-6), 0.3);

        // Uniform distributions for realistic positioning
        double bed_radius = portafilter_diameter_mm / 2000.0; // Convert to meters
        double bed_depth = bed_depth_mm / 1000.0;             // Convert to meters

        std::uniform_real_distribution<double> angle_dist(0.0, 2.0 * M_PI);
        std::uniform_real_distribution<double> radius_frac_dist(0.0, 1.0);
        std::uniform_real_distribution<double> depth_dist(0.0, bed_depth);

        int total_particles = 500; // More particles for realistic density

        std::cout << "Placing " << total_particles << " particles with realistic distribution..."
                  << std::endl;

        for (int i = 0; i < total_particles; i++) {
            // Generate realistic particle size (200-800 microns)
            double particle_size = size_dist(gen);
            particle_size = std::clamp(particle_size, 200e-6, 800e-6); // Clamp to realistic range

            // Generate random position (avoiding spiral)
            double r_fraction = std::sqrt(radius_frac_dist(gen)); // Square root for uniform area
            double r = bed_radius * r_fraction;
            double theta = angle_dist(gen);
            double z = depth_dist(gen);

            // Convert to Cartesian coordinates
            double x = r * std::cos(theta);
            double y = z; // Use Z (depth) as Y coordinate for 2D visualization

            // Add some size-based stratification
            double size_factor = (particle_size - 200e-6) / (800e-6 - 200e-6); // 0-1
            y = y * (0.7 + 0.6 * size_factor); // Larger particles appear "deeper"

            // Simple overlap avoidance
            bool valid_position = true;
            for (const auto& existing_particle : bed->particles()) {
                auto [ex, ey] = existing_particle->get_position();
                double distance = std::sqrt((x - ex) * (x - ex) + (y - ey) * (y - ey));
                if (distance < particle_size * 1.5) { // Allow some overlap for denser packing
                    valid_position = false;
                    break;
                }
            }

            if (valid_position) {
                auto particle = std::make_shared<CoffeeParticle>(x, y, particle_size / 2.0);
                bed->add_particle(particle);
            }
        }

        std::cout << "  Dose: " << dose_grams << "g" << std::endl;
        std::cout << "  Portafilter diameter: " << portafilter_diameter_mm << "mm" << std::endl;
        std::cout << "  Bed depth: " << bed_depth_mm << "mm" << std::endl;
        std::cout << "  Total particles placed: " << bed->particles().size() << std::endl;
        std::cout << "  Particle size range: 200-800 microns" << std::endl;
        std::cout << "  Distribution: Log-normal (realistic grind)" << std::endl;
        std::cout << "  Placement: Random distribution with size stratification" << std::endl;
        std::cout << "  Initial bed height: " << std::fixed << std::setprecision(2)
                  << bed->get_bed_height() * 1000 << "mm" << std::endl;
        std::cout << "  Initial porosity: " << std::fixed << std::setprecision(3)
                  << bed->get_porosity() << std::endl;
        std::cout << std::endl;

        return bed;
    }

    std::shared_ptr<WaterFlow> create_3d_water_flow() {
        std::cout << "Creating enhanced water flow field..." << std::endl;

        // Create enhanced grid for better resolution
        size_t width = 58; // X direction (mm)
        size_t depth = 15; // Depth through coffee bed (mm)

        auto flow = std::make_shared<WaterFlow>(width, depth);

        std::cout << "  Grid dimensions: " << width << "x" << depth << " cells" << std::endl;
        std::cout << "  Cell size: " << flow->get_cell_size() * 1000 << "mm" << std::endl;
        std::cout << "  Flow direction: Through coffee bed depth" << std::endl;
        std::cout << "  Enhanced resolution for realistic flow patterns" << std::endl;
        std::cout << std::endl;

        return flow;
    }

    std::shared_ptr<PhysicsSolver> create_physics_solver(std::shared_ptr<CoffeeBed> bed,
                                                         std::shared_ptr<WaterFlow> flow) {

        std::cout << "Setting up enhanced physics solver..." << std::endl;

        // Realistic espresso brewing parameters
        PhysicsSolver::Parameters params{
            .permeability = 8e-13,            // Slightly lower for finer grind
            .fluid_viscosity = 1e-3,          // Water viscosity at 95°C [Pa·s]
            .extraction_rate = 0.08,          // Higher rate for finer particles
            .temperature = 368.15,            // 95°C in Kelvin
            .inlet_pressure = 9e5,            // 9 bar inlet pressure [Pa]
            .outlet_pressure = 1e5,           // 1 bar outlet pressure [Pa]
            .particle_drag = 0.15,            // Higher drag for 3D
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

    void run_3d_simulation(std::shared_ptr<CoffeeBed> bed, std::shared_ptr<WaterFlow> flow,
                           std::shared_ptr<PhysicsSolver> solver) {

        std::cout << "Running realistic brewing simulation..." << std::endl;
        std::cout << "Simulating 30 seconds of enhanced espresso extraction..." << std::endl;

#ifdef SPROSIM_HAS_VTK
        std::cout << "VTK export enabled - saving realistic visualization data..." << std::endl;
        sprosim::visualization::VTKExporter vtk_exporter;
        vtk_exporter.set_output_directory("./demo_3d_output");
        vtk_exporter.start_time_series("brewing_3d_simulation.pvd", "timestep_3d");
#else
        std::cout << "VTK export disabled - install VTK for enhanced ParaView visualization"
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
        std::cout << "Enhanced VTK files saved to ./demo_3d_output/" << std::endl;
        std::cout << "Open brewing_3d_simulation.pvd in ParaView for realistic visualization"
                  << std::endl;
#endif
        std::cout << std::endl;
    }

    void display_3d_results(std::shared_ptr<CoffeeBed> bed, std::shared_ptr<WaterFlow> flow) {

        std::cout << "=== Realistic Brewing Results ===" << std::endl;
        std::cout << std::endl;

        // Coffee bed analysis with realistic insights
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
        const auto& particles = bed->particles();

        // Analyze particles by depth layers
        std::vector<std::vector<double>> layer_extractions(3); // Top, middle, bottom layers

        for (const auto& particle : particles) {
            auto [x, y] = particle->get_position();
            double extraction = particle->get_extraction_state();

            // Determine layer based on Y position (depth proxy)
            double normalized_depth = y / bed->get_bed_height();

            if (normalized_depth < 0.33) {
                layer_extractions[0].push_back(extraction); // Top layer
            } else if (normalized_depth < 0.67) {
                layer_extractions[1].push_back(extraction); // Middle layer
            } else {
                layer_extractions[2].push_back(extraction); // Bottom layer
            }
        }

        std::cout << "Particle Extraction by Depth Layers:" << std::endl;

        std::vector<std::string> layer_names = {"Top (inlet)", "Middle", "Bottom (outlet)"};
        for (size_t layer = 0; layer < 3; layer++) {
            if (!layer_extractions[layer].empty()) {
                double sum = 0.0;
                double min_ext = 1.0, max_ext = 0.0;
                for (double ext : layer_extractions[layer]) {
                    sum += ext;
                    min_ext = std::min(min_ext, ext);
                    max_ext = std::max(max_ext, ext);
                }
                double avg_ext = sum / layer_extractions[layer].size();

                std::cout << "  " << layer_names[layer] << " layer ("
                          << layer_extractions[layer].size() << " particles):" << std::endl;
                std::cout << "    Average extraction: " << std::fixed << std::setprecision(2)
                          << avg_ext * 100 << "%" << std::endl;
                std::cout << "    Range: " << std::fixed << std::setprecision(2) << min_ext * 100
                          << "% - " << max_ext * 100 << "%" << std::endl;
            }
        }
        std::cout << std::endl;

        // Overall extraction uniformity
        double min_extraction = 1.0, max_extraction = 0.0, avg_extraction = 0.0;
        for (const auto& particle : particles) {
            double extraction = particle->get_extraction_state();
            min_extraction = std::min(min_extraction, extraction);
            max_extraction = std::max(max_extraction, extraction);
            avg_extraction += extraction;
        }
        avg_extraction /= particles.size();

        std::cout << "Overall Extraction Distribution:" << std::endl;
        std::cout << "  Minimum extraction: " << std::fixed << std::setprecision(2)
                  << min_extraction * 100 << "%" << std::endl;
        std::cout << "  Maximum extraction: " << std::fixed << std::setprecision(2)
                  << max_extraction * 100 << "%" << std::endl;
        std::cout << "  Average extraction: " << std::fixed << std::setprecision(2)
                  << avg_extraction * 100 << "%" << std::endl;
        std::cout << "  Extraction uniformity: " << std::fixed << std::setprecision(3)
                  << (1.0 - (max_extraction - min_extraction)) << std::endl;
        std::cout << std::endl;

        // Brewing quality assessment with realistic considerations
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
        std::cout << "=== Enhanced ParaView Visualization Instructions ===" << std::endl;
        std::cout << "1. Install ParaView from https://www.paraview.org/" << std::endl;
        std::cout << "2. Open ParaView and load: ./demo_3d_output/brewing_3d_simulation.pvd"
                  << std::endl;
        std::cout << "3. Click 'Apply' to load the enhanced time series data" << std::endl;
        std::cout << "4. Create sphere glyphs to visualize particles with realistic distribution"
                  << std::endl;
        std::cout << "5. Color particles by 'extraction_state' to see extraction progress"
                  << std::endl;
        std::cout << "6. Use cross-sections (Filters → Slice) to see different depth layers"
                  << std::endl;
        std::cout << "7. Rotate view to see the realistic coffee bed distribution" << std::endl;
        std::cout << "8. Play animation to watch 30 seconds of realistic brewing" << std::endl;
        std::cout << std::endl;
        std::cout << "Advanced Analysis:" << std::endl;
        std::cout << "- Slice at different Y levels to see depth stratification" << std::endl;
        std::cout << "- Color by 'particle_size' to see grind distribution effects" << std::endl;
        std::cout << "- Use volume rendering for continuous extraction visualization" << std::endl;
        std::cout << "- Create streamlines to visualize realistic water flow paths" << std::endl;
        std::cout << std::endl;
#endif

        std::cout << "=== Realistic Simulation Complete ===" << std::endl;
        std::cout << "Realistic particle distribution with " << particles.size() << " particles"
                  << std::endl;
        std::cout << "Log-normal size distribution (200-800 microns)" << std::endl;
        std::cout << "Realistic spatial distribution with size stratification" << std::endl;
        std::cout << "Enhanced physics for improved brewing dynamics" << std::endl;
    }
};

int main() {
    try {
        Realistic3DEspressoDemo demo;
        demo.run();
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}