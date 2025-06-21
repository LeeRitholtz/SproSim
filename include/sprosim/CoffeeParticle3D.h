#pragma once
#include "sprosim/interfaces/IParticle.h"
#include <array>
#include <tuple>
#include <utility>
#include <vector>

namespace sprosim {

/**
 * @brief 3D Coffee particle implementation with full spatial positioning
 *
 * Extends the basic coffee particle to support 3D positioning with Z coordinate,
 * enabling realistic coffee bed depth simulation and proper cross-sectional analysis.
 */
class CoffeeParticle3D : public ICoffeeParticle {
  public:
    /**
     * @brief Construct a 3D coffee particle
     * @param x X coordinate in meters (horizontal position)
     * @param y Y coordinate in meters (horizontal position)
     * @param z Z coordinate in meters (depth through coffee bed)
     * @param radius Particle radius in meters
     */
    CoffeeParticle3D(double x, double y, double z, double radius);

    // IParticle interface (2D compatibility)
    std::pair<double, double> get_position() const override;
    double get_size() const override;

    // ICoffeeParticle interface
    double get_extraction_state() const override;
    double get_concentration() const override;
    void update_extraction(double delta_conc, double dt) override;
    void apply_force(double fx, double fy) override;

    // 3D-specific methods
    /**
     * @brief Get full 3D position
     * @return Tuple containing (x, y, z) coordinates in meters
     */
    std::tuple<double, double, double> get_position_3d() const;

    /**
     * @brief Get Z coordinate (depth through coffee bed)
     * @return Z coordinate in meters (0 = top of bed, positive = deeper)
     */
    double get_depth() const;

    /**
     * @brief Apply 3D force to particle
     * @param fx Force in X direction [N]
     * @param fy Force in Y direction [N]
     * @param fz Force in Z direction [N] (vertical through bed)
     */
    void apply_force_3d(double fx, double fy, double fz);

    /**
     * @brief Update 3D position directly
     * @param dx Change in X position [m]
     * @param dy Change in Y position [m]
     * @param dz Change in Z position [m]
     */
    void update_position_3d(double dx, double dy, double dz);

    /**
     * @brief Set absolute 3D position
     * @param x New X coordinate [m]
     * @param y New Y coordinate [m]
     * @param z New Z coordinate [m]
     */
    void set_position_3d(double x, double y, double z);

    /**
     * @brief Get particle volume (for 3D packing calculations)
     * @return Volume in cubic meters
     */
    double get_volume() const;

    /**
     * @brief Get distance to another 3D particle
     * @param other Another 3D particle
     * @return 3D Euclidean distance in meters
     */
    double distance_to(const CoffeeParticle3D& other) const;

    /**
     * @brief Check if particle overlaps with another in 3D space
     * @param other Another 3D particle
     * @return True if particles overlap
     */
    bool overlaps_with(const CoffeeParticle3D& other) const;

    /**
     * @brief Get local packing density around this particle
     * @param neighbors Vector of nearby particles
     * @param search_radius Radius to search for neighbors [m]
     * @return Local packing fraction (0-1)
     */
    double get_local_packing_density(const std::vector<CoffeeParticle3D>& neighbors,
                                     double search_radius) const;

    /**
     * @brief Get particle mass
     * @return Particle mass in kg
     */
    double get_mass() const;

  private:
    std::array<double, 3> position; // [x, y, z] coordinates in meters
    double radius;                  // Particle radius in meters
    double extraction_state;        // Extraction progress (0-1)
    double concentration;           // Local dissolved concentration [kg/m³]
    std::array<double, 3> velocity; // Particle velocity [m/s] for physics

    // Physical constants
    static constexpr double max_extractable_ = 0.22;    // 22% extractable content
    static constexpr double particle_density_ = 1300.0; // kg/m³ (coffee density)
    static constexpr double pi_ = 3.14159265358979323846;

    // Simulation constants
    static constexpr double timestep_dt_ = 0.001;  // Integration timestep [s]
    static constexpr double max_bed_depth_ = 0.02; // Maximum bed depth [m] (20mm)
    static constexpr double max_portafilter_radius_ =
        0.029;                                        // Portafilter radius [m] (58mm diameter / 2)
    static constexpr double velocity_damping_ = 0.95; // Velocity damping factor for friction

    // Helper methods
    double calculate_mass() const;
    void apply_boundary_constraints();
    void update_physics_state(double dt);
};

/**
 * @brief Utility functions for 3D particle operations
 */
namespace particle_3d_utils {

/**
 * @brief Generate random 3D particle distribution
 * @param num_particles Number of particles to generate
 * @param bed_radius Portafilter radius [m]
 * @param bed_depth Coffee bed depth [m]
 * @param size_distribution Particle size distribution parameters
 * @return Vector of 3D positioned particles
 */
std::vector<CoffeeParticle3D>
generate_random_3d_distribution(int num_particles, double bed_radius, double bed_depth,
                                const std::pair<double, double>& size_range = {200e-6, 800e-6});

/**
 * @brief Simulate particle settling under gravity
 * @param particles Vector of particles to settle
 * @param gravity_acceleration Gravitational acceleration [m/s²]
 * @param settling_time Time to simulate settling [s]
 * @param time_step Integration time step [s]
 */
void simulate_gravity_settling(std::vector<CoffeeParticle3D>& particles,
                               double gravity_acceleration = 9.81, double settling_time = 1.0,
                               double time_step = 0.001);

/**
 * @brief Detect and resolve particle overlaps
 * @param particles Vector of particles to check
 * @param max_iterations Maximum collision resolution iterations
 * @return Number of overlaps resolved
 */
int resolve_particle_overlaps(std::vector<CoffeeParticle3D>& particles, int max_iterations = 100);

/**
 * @brief Calculate bed porosity in 3D space
 * @param particles Vector of particles
 * @param bed_volume Total bed volume [m³]
 * @return Porosity fraction (0-1)
 */
double calculate_bed_porosity(const std::vector<CoffeeParticle3D>& particles, double bed_volume);

/**
 * @brief Create particle size stratification (larger particles settle lower)
 * @param particles Vector of particles to stratify
 * @param stratification_strength Strength of size segregation (0-1)
 */
void apply_size_stratification(std::vector<CoffeeParticle3D>& particles,
                               double stratification_strength = 0.3);
} // namespace particle_3d_utils

} // namespace sprosim