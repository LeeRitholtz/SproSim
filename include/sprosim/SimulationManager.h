#pragma once

#include "sprosim/CoffeeBed.h"
#include "sprosim/Parameters.h"
#include "sprosim/PhysicsSolver.h"
#include "sprosim/WaterFlow.h"

#ifdef SPROSIM_HAS_VTK
#include "../visualization/paraview/exporters/VTKExporter.h"
#endif

#include <functional>
#include <memory>
#include <string>

namespace sprosim {

/**
 * @brief High-level simulation controller for coffee brewing processes
 *
 * SimulationManager orchestrates complete brewing simulations by managing
 * the simulation lifecycle, termination conditions, progress reporting,
 * and result aggregation. It uses PhysicsSolver for individual timestep
 * physics calculations while handling the broader simulation management.
 *
 * Supports multiple termination conditions:
 * - Time-based: Run for specified duration
 * - Yield-based: Stop when target extraction yield reached
 * - Convergence-based: Stop when extraction rate falls below threshold
 *
 * @example
 * ```cpp
 * SimulationManager::Configuration config{
 *     .dose = 18.0,              // grams
 *     .particle_count = 350,
 *     .brewing_time = 30.0,      // seconds
 *     .target_yield = 0.20,      // 20% extraction
 *     .timestep = 0.01,          // seconds
 *     .export_vtk = true
 * };
 *
 * SimulationManager manager(config);
 * auto results = manager.run_simulation();
 * ```
 */
class SimulationManager {
  public:
    /**
     * @brief Simulation configuration parameters
     */
    struct Configuration {
        // Core recipe parameters
        double dose = 18.0;         // Coffee dose in grams
        double target_yield = 36.0; // Target output weight in grams
        double brewing_time = 30.0; // Maximum brewing time in seconds

        // Coffee bed parameters
        double portafilter_diameter = 58.0; // Portafilter diameter in mm
        double tamping_force = 30.0;        // Tamping force in lbs
        double grind_size_min = 300.0;      // Minimum grind size in microns
        double grind_size_max = 600.0;      // Maximum grind size in microns

        // Brewing conditions
        double pressure_bar = 9.0;         // Brewing pressure in bar
        double temperature_celsius = 95.0; // Water temperature in °C
        double extraction_rate = 0.01;     // Extraction rate constant

        // Grid parameters
        int grid_width = 58;  // Flow field width
        int grid_height = 30; // Flow field height
        // TODO: particle_count should be randomized between a reasonable range
        int particle_count = 350; // Number of particles to simulate
        double timestep = 0.01;   // Simulation timestep in seconds

        // Output parameters
        bool export_vtk = false;                            // Enable VTK export
        std::string vtk_output_dir = "./extraction_output"; // VTK output directory
        bool use_3d_particles = false;                      // Use 3D particles

        // Termination conditions
        bool enable_yield_termination = false;       // Stop at target yield
        bool enable_convergence_termination = false; // Stop when extraction rate low
        double convergence_threshold = 0.001;        // Convergence threshold
    };

    /**
     * @brief Simulation results and statistics
     */
    struct Results {
        double final_extraction_yield = 0.0; // Final extraction yield (0-1)
        double total_dissolved_solids = 0.0; // TDS percentage
        double final_porosity = 0.0;         // Final bed porosity
        double simulation_time = 0.0;        // Total simulation time
        int total_steps = 0;                 // Number of simulation steps
        std::string termination_reason;      // Why simulation ended

        // Quality assessment
        bool is_under_extracted = false; // Extraction < 18%
        bool is_over_extracted = false;  // Extraction > 24%
        bool is_well_extracted = false;  // 18-24% range
    };

    /**
     * @brief Progress reporting callback function type
     *
     * Called periodically during simulation with:
     * - current_time: Current simulation time
     * - progress: Progress percentage (0-100)
     * - extraction_yield: Current extraction yield (0-1)
     */
    using ProgressCallback =
        std::function<void(double current_time, double progress, double extraction_yield)>;

    /**
     * @brief Constructor with simulation configuration and physics parameters
     * @param config Simulation configuration parameters
     * @param physics_params Physics parameters for simulation
     */
    SimulationManager(const Configuration& config, const Parameters& physics_params);

    /**
     * @brief Constructor with automatic physics parameter conversion
     * @param config Simulation configuration parameters
     *
     * Automatically converts brewing parameters to physics parameters.
     */
    explicit SimulationManager(const Configuration& config);

    /**
     * @brief Run complete brewing simulation
     * @return Simulation results and statistics
     */
    Results run_simulation();

    /**
     * @brief Set progress reporting callback
     * @param callback Function to call for progress updates
     */
    void set_progress_callback(const ProgressCallback& callback);

    /**
     * @brief Get current simulation configuration
     * @return Configuration reference
     */
    const Configuration& get_configuration() const;

  private:
    Configuration config_;
    Parameters physics_params_;
    ProgressCallback progress_callback_;

    // Simulation components
    std::shared_ptr<CoffeeBed> coffee_bed_;
    std::shared_ptr<WaterFlow> water_flow_;
    std::shared_ptr<PhysicsSolver> physics_solver_;

#ifdef SPROSIM_HAS_VTK
    std::unique_ptr<sprosim::visualization::VTKExporter> vtk_exporter_;
#endif

    // Setup methods
    void initialize_simulation();
    std::shared_ptr<CoffeeBed> create_coffee_bed();
    std::shared_ptr<WaterFlow> create_water_flow();
    std::shared_ptr<PhysicsSolver> create_physics_solver();

    // Parameter conversion
    Parameters convert_to_physics_params(const Configuration& config);

    // Simulation control
    bool check_termination_conditions(double current_time, double extraction_yield);
    void report_progress(double current_time, double progress, double extraction_yield);
    void setup_vtk_export();
    void export_vtk_timestep(double current_time, int step);
    void finalize_vtk_export();

    // Results processing
    Results generate_results(double simulation_time, int total_steps,
                             const std::string& termination_reason);
    void assess_extraction_quality(Results& results);
};

} // namespace sprosim
