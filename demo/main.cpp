#include <iostream>
#include <memory>
#include <iomanip>
#include <vector>
#include <cmath>

#include "sprosim/CoffeeBed.h"
#include "sprosim/CoffeeParticle.h"
#include "sprosim/WaterFlow.h"
#include "sprosim/PhysicsSolver.h"

#ifdef SPROSIM_HAS_VTK
#include "paraview/exporters/VTKExporter.h"
#endif

using namespace sprosim;

class EspressoBrewingDemo {
public:
    EspressoBrewingDemo() {
        std::cout << "=== SproSim Coffee Brewing Simulation Demo ===" << std::endl;
        std::cout << "Simulating espresso extraction with realistic parameters" << std::endl;
        std::cout << std::endl;
    }

    void run() {
        // Create coffee bed (18g dose, 58mm portafilter)
        auto coffee_bed = create_coffee_bed();
        
        // Create water flow field
        auto water_flow = create_water_flow();
        
        // Setup physics solver with realistic espresso parameters
        auto solver = create_physics_solver(coffee_bed, water_flow);
        
        // Run simulation
        run_simulation(coffee_bed, water_flow, solver);
        
        // Display final results
        display_results(coffee_bed, water_flow);
    }

private:
    std::shared_ptr<CoffeeBed> create_coffee_bed() {
        std::cout << "Creating coffee bed..." << std::endl;
        
        // Standard espresso parameters
        double dose_grams = 18.0;
        double portafilter_diameter_mm = 58.0;
        
        auto bed = std::make_shared<CoffeeBed>(dose_grams, portafilter_diameter_mm);
        
        // Add coffee particles with realistic size distribution
        std::cout << "Adding coffee particles with realistic grind distribution..." << std::endl;
        
        // Simulate a medium-fine grind with particle sizes from 200-800 microns
        std::vector<double> particle_sizes = {
            0.0002, 0.0003, 0.0004, 0.0005, 0.0006, 0.0007, 0.0008  // 200-800 microns
        };
        
        // Create particles distributed across the bed
        double bed_radius = portafilter_diameter_mm / 2000.0; // Convert to meters
        int particles_per_size = 50; // Total: 350 particles
        
        for (double size : particle_sizes) {
            for (int i = 0; i < particles_per_size; i++) {
                // Distribute particles randomly across the bed
                double angle = 2.0 * M_PI * i / particles_per_size;
                double radius_fraction = sqrt(static_cast<double>(i) / particles_per_size);
                double x = bed_radius * radius_fraction * cos(angle);
                double y = bed_radius * radius_fraction * sin(angle);
                
                auto particle = std::make_shared<CoffeeParticle>(x, y, size / 2.0); // radius = diameter/2
                bed->add_particle(particle);
            }
        }
        
        std::cout << "  Dose: " << dose_grams << "g" << std::endl;
        std::cout << "  Portafilter diameter: " << portafilter_diameter_mm << "mm" << std::endl;
        std::cout << "  Total particles: " << bed->get_particles().size() << std::endl;
        std::cout << "  Initial bed height: " << std::fixed << std::setprecision(2) 
                  << bed->get_bed_height() * 1000 << "mm" << std::endl;
        std::cout << "  Initial porosity: " << std::fixed << std::setprecision(3) 
                  << bed->get_porosity() << std::endl;
        std::cout << std::endl;
        
        return bed;
    }
    
    std::shared_ptr<WaterFlow> create_water_flow() {
        std::cout << "Creating water flow field..." << std::endl;
        
        // Create a 58x30 grid (58mm wide, 30mm tall)
        size_t width = 58;   // 1mm per cell
        size_t height = 30;  // 1mm per cell
        
        auto flow = std::make_shared<WaterFlow>(width, height);
        
        std::cout << "  Grid dimensions: " << width << "x" << height << " cells" << std::endl;
        std::cout << "  Cell size: " << flow->get_cell_size() * 1000 << "mm" << std::endl;
        std::cout << std::endl;
        
        return flow;
    }
    
    std::shared_ptr<PhysicsSolver> create_physics_solver(
        std::shared_ptr<CoffeeBed> bed, 
        std::shared_ptr<WaterFlow> flow) {
        
        std::cout << "Setting up physics solver..." << std::endl;
        
        // Realistic espresso brewing parameters
        PhysicsSolver::Parameters params{
            .permeability = 1e-12,           // Typical coffee bed permeability [m²]
            .fluid_viscosity = 1e-3,         // Water viscosity at 95°C [Pa·s]
            .extraction_rate = 0.05,         // Extraction rate constant [1/s]
            .temperature = 368.15,           // 95°C in Kelvin
            .inlet_pressure = 9e5,           // 9 bar inlet pressure [Pa]
            .outlet_pressure = 1e5,          // 1 bar outlet pressure [Pa]
            .particle_drag = 0.1,            // Particle drag coefficient
            .flow_resistance = 0.2,          // Flow resistance factor
            .saturation_concentration = 0.15, // Max concentration [kg/m³]
            .temperature_factor = 0.008      // Temperature dependence [K⁻¹]
        };
        
        auto solver = std::make_shared<PhysicsSolver>(bed, flow, params);
        
        std::cout << "  Brewing temperature: " << params.temperature - 273.15 << "°C" << std::endl;
        std::cout << "  Brewing pressure: " << (params.inlet_pressure - params.outlet_pressure) / 1e5 << " bar" << std::endl;
        std::cout << "  Extraction rate: " << params.extraction_rate << " s⁻¹" << std::endl;
        std::cout << std::endl;
        
        return solver;
    }
    
    void run_simulation(
        std::shared_ptr<CoffeeBed> bed,
        std::shared_ptr<WaterFlow> flow,
        std::shared_ptr<PhysicsSolver> solver) {
        
        std::cout << "Running brewing simulation..." << std::endl;
        std::cout << "Simulating 30 seconds of extraction..." << std::endl;
        
#ifdef SPROSIM_HAS_VTK
        std::cout << "VTK export enabled - saving visualization data..." << std::endl;
        sprosim::visualization::VTKExporter vtk_exporter;
        vtk_exporter.set_output_directory("./demo_output");
        vtk_exporter.start_time_series("brewing_simulation.pvd", "timestep");
#else
        std::cout << "VTK export disabled - install VTK for ParaView visualization" << std::endl;
#endif
        std::cout << std::endl;
        
        const double total_time = 30.0;    // 30 seconds
        const double dt = 0.01;            // 10ms time steps
        const int total_steps = static_cast<int>(total_time / dt);
        const int report_interval = total_steps / 10; // Report every 10%
        const int vtk_export_interval = total_steps / 100; // Export every 1%
        
        double current_time = 0.0;
        
        for (int step = 0; step < total_steps; step++) {
            solver->simulate_step(dt);
            current_time += dt;
            
#ifdef SPROSIM_HAS_VTK
            // Export VTK data every 1% of simulation
            if (step % vtk_export_interval == 0 || step == total_steps - 1) {
                vtk_exporter.add_timestep_to_series(bed, flow, current_time, step);
            }
#endif
            
            // Report progress every 10%
            if (step % report_interval == 0 || step == total_steps - 1) {
                double progress = static_cast<double>(step) / total_steps * 100.0;
                double extraction_yield = bed->get_extraction_yield();
                double tds = bed->get_total_dissolved_solids();
                
                std::cout << std::fixed << std::setprecision(1);
                std::cout << "  " << std::setw(5) << progress << "% complete";
                std::cout << " | Time: " << std::setw(4) << current_time << "s";
                std::cout << " | Extraction: " << std::setprecision(2) << std::setw(5) << extraction_yield * 100 << "%";
                std::cout << " | TDS: " << std::setprecision(3) << std::setw(6) << tds * 100 << "%";
                std::cout << std::endl;
            }
        }
        
#ifdef SPROSIM_HAS_VTK
        vtk_exporter.finalize_time_series();
        std::cout << "VTK files saved to ./demo_output/" << std::endl;
        std::cout << "Open brewing_simulation.pvd in ParaView for visualization" << std::endl;
#endif
        std::cout << std::endl;
    }
    
    void display_results(
        std::shared_ptr<CoffeeBed> bed,
        std::shared_ptr<WaterFlow> flow) {
        
        std::cout << "=== Final Brewing Results ===" << std::endl;
        std::cout << std::endl;
        
        // Coffee bed analysis
        std::cout << "Coffee Bed Analysis:" << std::endl;
        std::cout << "  Final extraction yield: " << std::fixed << std::setprecision(2) 
                  << bed->get_extraction_yield() * 100 << "%" << std::endl;
        std::cout << "  Total dissolved solids: " << std::fixed << std::setprecision(3) 
                  << bed->get_total_dissolved_solids() * 100 << "%" << std::endl;
        std::cout << "  Bed compaction: " << std::fixed << std::setprecision(3) 
                  << bed->get_compaction() << std::endl;
        std::cout << "  Final porosity: " << std::fixed << std::setprecision(3) 
                  << bed->get_porosity() << std::endl;
        std::cout << "  Final bed height: " << std::fixed << std::setprecision(2) 
                  << bed->get_bed_height() * 1000 << "mm" << std::endl;
        std::cout << std::endl;
        
        // Flow field analysis
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
        
        std::cout << "Flow Field Analysis:" << std::endl;
        std::cout << "  Average velocity: " << std::fixed << std::setprecision(4) 
                  << avg_velocity * 1000 << " mm/s" << std::endl;
        std::cout << "  Maximum velocity: " << std::fixed << std::setprecision(4) 
                  << max_velocity * 1000 << " mm/s" << std::endl;
        std::cout << std::endl;
        
        // Particle analysis
        const auto& particles = bed->get_particles();
        double min_extraction = 1.0, max_extraction = 0.0, avg_extraction = 0.0;
        
        for (const auto& particle : particles) {
            double extraction = particle->get_extraction_state();
            min_extraction = std::min(min_extraction, extraction);
            max_extraction = std::max(max_extraction, extraction);
            avg_extraction += extraction;
        }
        avg_extraction /= particles.size();
        
        std::cout << "Particle Extraction Distribution:" << std::endl;
        std::cout << "  Minimum extraction: " << std::fixed << std::setprecision(2) 
                  << min_extraction * 100 << "%" << std::endl;
        std::cout << "  Maximum extraction: " << std::fixed << std::setprecision(2) 
                  << max_extraction * 100 << "%" << std::endl;
        std::cout << "  Average extraction: " << std::fixed << std::setprecision(2) 
                  << avg_extraction * 100 << "%" << std::endl;
        std::cout << std::endl;
        
        // Brewing quality assessment
        std::cout << "Brewing Quality Assessment:" << std::endl;
        double extraction_yield = bed->get_extraction_yield();
        if (extraction_yield < 0.18) {
            std::cout << "  Result: UNDER-EXTRACTED (sour, weak)" << std::endl;
        } else if (extraction_yield > 0.24) {
            std::cout << "  Result: OVER-EXTRACTED (bitter, harsh)" << std::endl;
        } else {
            std::cout << "  Result: WELL-EXTRACTED (balanced, sweet)" << std::endl;
        }
        
        double extraction_uniformity = 1.0 - (max_extraction - min_extraction);
        if (extraction_uniformity > 0.8) {
            std::cout << "  Uniformity: EXCELLENT (even extraction)" << std::endl;
        } else if (extraction_uniformity > 0.6) {
            std::cout << "  Uniformity: GOOD (mostly even)" << std::endl;
        } else {
            std::cout << "  Uniformity: POOR (uneven extraction)" << std::endl;
        }
        
        std::cout << std::endl;
        
#ifdef SPROSIM_HAS_VTK
        std::cout << "=== ParaView Visualization Instructions ===" << std::endl;
        std::cout << "1. Install ParaView from https://www.paraview.org/" << std::endl;
        std::cout << "2. Open ParaView and load: ./demo_output/brewing_simulation.pvd" << std::endl;
        std::cout << "3. Click 'Apply' to load the time series data" << std::endl;
        std::cout << "4. Use the play button to animate the extraction process" << std::endl;
        std::cout << "5. Color particles by 'extraction_state' to see extraction progress" << std::endl;
        std::cout << "6. Add velocity vectors to visualize water flow" << std::endl;
        std::cout << std::endl;
#endif
        
        std::cout << "=== Simulation Complete ===" << std::endl;
    }
};

int main() {
    try {
        EspressoBrewingDemo demo;
        demo.run();
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}