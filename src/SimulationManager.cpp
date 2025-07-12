#include "sprosim/SimulationManager.h"
#include "sprosim/CoffeeParticle.h"
#include "sprosim/CoffeeParticle3D.h"
#include <cmath>
#include <iomanip>
#include <iostream>

namespace sprosim {

SimulationManager::SimulationManager(const Configuration& config, const Parameters& physics_params)
    : config_(config), physics_params_(physics_params) {
    // Constructor just stores parameters, actual initialization happens in run_simulation()
}

SimulationManager::SimulationManager(const Configuration& config)
    : config_(config), physics_params_(convert_to_physics_params(config)) {
    // Constructor with automatic physics parameter conversion
}

SimulationManager::Results SimulationManager::run_simulation() {
    initialize_simulation();

    const int total_steps = static_cast<int>(config_.brewing_time / config_.timestep);
    const int report_interval = std::max(1, total_steps / 10); // Report 10 times during sim

    double current_time = 0.0;
    double current_yield_grams = 0.0; // Track output weight in grams
    std::string termination_reason = "time_limit";
    int step = 0;

    setup_vtk_export();

    for (step = 0; step < total_steps; step++) {
        physics_solver_->simulate_step(config_.timestep);
        current_time += config_.timestep;

        // Calculate current output yield in grams
        // TODO: This is simplified - need proper liquid flow tracking
        double extraction_fraction = coffee_bed_->get_extraction_yield();
        current_yield_grams =
            config_.dose * extraction_fraction * 2.0; // Rough 1:2 ratio approximation

        // Check termination conditions
        if (check_termination_conditions(current_time, current_yield_grams)) {
            if (config_.enable_yield_termination && current_yield_grams >= config_.target_yield) {
                termination_reason = "target_yield_reached";
            } else if (config_.enable_convergence_termination) {
                termination_reason = "converged";
            }
            break;
        }

        // Progress reporting
        if (step % report_interval == 0 || step == total_steps - 1) {
            double progress = 100.0 * step / total_steps;
            report_progress(current_time, progress, extraction_fraction);
        }

        // VTK export
        export_vtk_timestep(current_time, step);
    }

    finalize_vtk_export();
    return generate_results(current_time, step, termination_reason);
}

void SimulationManager::set_progress_callback(const ProgressCallback& callback) {
    progress_callback_ = callback;
}

const SimulationManager::Configuration& SimulationManager::get_configuration() const {
    return config_;
}

void SimulationManager::initialize_simulation() {
    coffee_bed_ = create_coffee_bed();
    water_flow_ = create_water_flow();
    physics_solver_ = create_physics_solver();
}

std::shared_ptr<CoffeeBed> SimulationManager::create_coffee_bed() {
    auto bed = std::make_shared<CoffeeBed>(config_.dose, config_.portafilter_diameter);

    // Apply tamping force (convert lbs to Pa for compaction)
    double tamping_pa = config_.tamping_force * 6895.0; // lbs to Pa
    bed->update_compaction(tamping_pa);

    // Create particles with grind distribution
    double bed_radius = config_.portafilter_diameter / 2000.0; // mm to m

    for (int i = 0; i < config_.particle_count; i++) {
        // Distribute particles in circular pattern
        double angle = 2.0 * M_PI * i / config_.particle_count;
        double radius_fraction = std::sqrt(static_cast<double>(i) / config_.particle_count);
        double x = bed_radius * radius_fraction * std::cos(angle);
        double y = bed_radius * radius_fraction * std::sin(angle);

        // Particle size distribution within grind range
        double size_microns =
            config_.grind_size_min + (config_.grind_size_max - config_.grind_size_min) *
                                         (static_cast<double>(i) / config_.particle_count);
        double radius = size_microns * 1e-6 / 2.0; // microns to meters, diameter to radius

        if (config_.use_3d_particles) {
            double z = 0.0; // Start at bed bottom
            auto particle = std::make_shared<CoffeeParticle3D>(x, y, z, radius);
            bed->add_particle(particle);
        } else {
            auto particle = std::make_shared<CoffeeParticle>(x, y, radius);
            bed->add_particle(particle);
        }
    }

    return bed;
}

std::shared_ptr<WaterFlow> SimulationManager::create_water_flow() {
    return std::make_shared<WaterFlow>(config_.grid_width, config_.grid_height);
}

std::shared_ptr<PhysicsSolver> SimulationManager::create_physics_solver() {
    return std::make_shared<PhysicsSolver>(coffee_bed_, water_flow_, physics_params_);
}

Parameters SimulationManager::convert_to_physics_params(const Configuration& config) {
    Parameters physics_params;

    // Convert user-friendly brewing parameters to physics parameters
    physics_params.permeability = 1e-12;                     // Default bed permeability [m²]
    physics_params.fluid_viscosity = 1e-3;                   // Water viscosity [Pa·s]
    physics_params.extraction_rate = config.extraction_rate; // [1/s]
    physics_params.temperature = config.temperature_celsius + 273.15; // °C to K
    physics_params.inlet_pressure = config.pressure_bar * 1e5 + 1e5;  // bar to Pa + atmospheric
    physics_params.outlet_pressure = 1e5;                             // Atmospheric pressure [Pa]
    physics_params.particle_drag = 0.1;             // Default particle drag coefficient
    physics_params.flow_resistance = 0.2;           // Default flow resistance
    physics_params.saturation_concentration = 0.15; // Default saturation [kg/m³]
    physics_params.temperature_factor = 0.008;      // Default Arrhenius factor [K⁻¹]

    return physics_params;
}

bool SimulationManager::check_termination_conditions(double current_time,
                                                     double current_yield_grams) {
    if (config_.enable_yield_termination && current_yield_grams >= config_.target_yield) {
        return true;
    }

    // TODO: Add convergence termination logic
    if (config_.enable_convergence_termination) {
        // Check if extraction rate has fallen below threshold
        // This would require tracking extraction rate over time
    }

    return false;
}

void SimulationManager::report_progress(double current_time, double progress,
                                        double extraction_yield) {
    if (progress_callback_) {
        progress_callback_(current_time, progress, extraction_yield);
    }
}

void SimulationManager::setup_vtk_export() {
#ifdef SPROSIM_HAS_VTK
    if (config_.export_vtk) {
        vtk_exporter_ = std::make_unique<sprosim::visualization::VTKExporter>();
        vtk_exporter_->set_output_directory(config_.vtk_output_dir);
        vtk_exporter_->start_time_series("brewing.pvd", "timestep");
    }
#endif
}

void SimulationManager::export_vtk_timestep(double current_time, int step) {
#ifdef SPROSIM_HAS_VTK
    if (config_.export_vtk && vtk_exporter_) {
        // Export every 50th step to avoid too many files
        const int total_steps = static_cast<int>(config_.brewing_time / config_.timestep);
        if (step % (total_steps / 50 + 1) == 0) {
            vtk_exporter_->add_timestep_to_series(coffee_bed_, water_flow_, current_time, step);
        }
    }
#endif
}

void SimulationManager::finalize_vtk_export() {
#ifdef SPROSIM_HAS_VTK
    if (config_.export_vtk && vtk_exporter_) {
        vtk_exporter_->finalize_time_series();
    }
#endif
}

SimulationManager::Results
SimulationManager::generate_results(double simulation_time, int total_steps,
                                    const std::string& termination_reason) {
    Results results;

    results.final_extraction_yield = coffee_bed_->get_extraction_yield();
    results.total_dissolved_solids = coffee_bed_->get_total_dissolved_solids();
    results.final_porosity = coffee_bed_->get_porosity();
    results.simulation_time = simulation_time;
    results.total_steps = total_steps;
    results.termination_reason = termination_reason;

    assess_extraction_quality(results);

    return results;
}

void SimulationManager::assess_extraction_quality(Results& results) {
    double extraction_percentage = results.final_extraction_yield * 100.0;

    if (extraction_percentage < 18.0) {
        results.is_under_extracted = true;
        results.is_over_extracted = false;
        results.is_well_extracted = false;
    } else if (extraction_percentage > 24.0) {
        results.is_under_extracted = false;
        results.is_over_extracted = true;
        results.is_well_extracted = false;
    } else {
        results.is_under_extracted = false;
        results.is_over_extracted = false;
        results.is_well_extracted = true;
    }
}

} // namespace sprosim
