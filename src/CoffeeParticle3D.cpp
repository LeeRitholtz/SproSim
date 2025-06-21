#include "sprosim/CoffeeParticle3D.h"
#include <algorithm>
#include <cmath>
#include <random>

namespace sprosim {

CoffeeParticle3D::CoffeeParticle3D(double x, double y, double z, double radius)
    : position{x, y, z},
      radius(radius),
      extraction_state(0.0),
      concentration(0.0),
      velocity{0.0, 0.0, 0.0}
{}

std::pair<double, double> CoffeeParticle3D::get_position() const {
    return {position[0], position[1]};
}

std::tuple<double, double, double> CoffeeParticle3D::get_position_3d() const {
    return {position[0], position[1], position[2]};
}

double CoffeeParticle3D::get_depth() const {
    return position[2];
}

double CoffeeParticle3D::get_size() const {
    return radius * 2.0;
}

double CoffeeParticle3D::get_extraction_state() const {
    return extraction_state;
}

double CoffeeParticle3D::get_concentration() const {
    return concentration;
}

void CoffeeParticle3D::update_extraction(double delta_conc, double dt) {
    concentration += delta_conc;
    extraction_state = std::min(1.0, extraction_state + delta_conc / max_extractable_);
}

void CoffeeParticle3D::apply_force(double fx, double fy) {
    apply_force_3d(fx, fy, 0.0);
}

void CoffeeParticle3D::apply_force_3d(double fx, double fy, double fz) {
    const double dt = timestep_dt_;  // 1ms timestep
    double mass = calculate_mass();
    
    // F = ma -> a = F/m
    double ax = fx / mass;
    double ay = fy / mass;
    double az = fz / mass;
    
    // Update velocity
    velocity[0] += ax * dt;
    velocity[1] += ay * dt;
    velocity[2] += az * dt;
    
    // Update position using Verlet integration
    position[0] += velocity[0] * dt + 0.5 * ax * dt * dt;
    position[1] += velocity[1] * dt + 0.5 * ay * dt * dt;
    position[2] += velocity[2] * dt + 0.5 * az * dt * dt;
    
    // Apply boundary constraints
    apply_boundary_constraints();
}

void CoffeeParticle3D::update_position_3d(double dx, double dy, double dz) {
    position[0] += dx;
    position[1] += dy;
    position[2] += dz;
    apply_boundary_constraints();
}

void CoffeeParticle3D::set_position_3d(double x, double y, double z) {
    position[0] = x;
    position[1] = y;
    position[2] = z;
    apply_boundary_constraints();
}

double CoffeeParticle3D::get_volume() const {
    return (4.0 / 3.0) * pi_ * radius * radius * radius;
}

double CoffeeParticle3D::distance_to(const CoffeeParticle3D& other) const {
    auto [ox, oy, oz] = other.get_position_3d();
    double dx = position[0] - ox;
    double dy = position[1] - oy;
    double dz = position[2] - oz;
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

bool CoffeeParticle3D::overlaps_with(const CoffeeParticle3D& other) const {
    double distance = distance_to(other);
    double combined_radius = radius + other.radius;
    return distance < combined_radius;
}

double CoffeeParticle3D::get_local_packing_density(
    const std::vector<CoffeeParticle3D>& neighbors, 
    double search_radius) const {
    
    double search_volume = (4.0 / 3.0) * pi_ * search_radius * search_radius * search_radius;
    double particle_volume_sum = 0.0;
    
    for (const auto& neighbor : neighbors) {
        if (distance_to(neighbor) <= search_radius) {
            particle_volume_sum += neighbor.get_volume();
        }
    }
    
    return particle_volume_sum / search_volume;
}

double CoffeeParticle3D::get_mass() const {
    return particle_density_ * get_volume();
}

double CoffeeParticle3D::calculate_mass() const {
    return get_mass();
}

void CoffeeParticle3D::apply_boundary_constraints() {
    // Ensure particle stays within reasonable bounds
    // Z should be non-negative (above bed bottom)
    position[2] = std::max(0.0, position[2]);
    
    // Limit maximum depth (prevent particles from falling through bed)
    position[2] = std::min(max_bed_depth_, position[2]); // 20mm max depth
    
    // Radial constraint (stay within portafilter)
    double r = std::sqrt(position[0] * position[0] + position[1] * position[1]);
    const double max_radius = max_portafilter_radius_; // 58mm diameter / 2
    if (r > max_radius) {
        double scale = max_radius / r;
        position[0] *= scale;
        position[1] *= scale;
    }
}

void CoffeeParticle3D::update_physics_state(double dt) {
    // Apply damping to velocity (friction)
    const double damping = velocity_damping_;
    velocity[0] *= damping;
    velocity[1] *= damping;
    velocity[2] *= damping;
}

namespace particle_3d_utils {

std::vector<CoffeeParticle3D> generate_random_3d_distribution(
    int num_particles,
    double bed_radius,
    double bed_depth,
    const std::pair<double, double>& size_range) {
    
    std::vector<CoffeeParticle3D> particles;
    particles.reserve(num_particles);
    
    // Random number generation
    std::random_device rd;
    std::mt19937 gen(42); // Fixed seed for reproducibility
    
    // Size distribution (log-normal)
    double mean_size = (size_range.first + size_range.second) / 2.0;
    std::lognormal_distribution<double> size_dist(std::log(mean_size), 0.3);
    
    // Position distributions - using normal distribution for X,Y coordinates
    std::normal_distribution<double> x_dist(0.0, bed_radius * 0.3);
    std::normal_distribution<double> y_dist(0.0, bed_radius * 0.3);
    std::uniform_real_distribution<double> depth_dist(0.0, bed_depth);
    
    for (int i = 0; i < num_particles; i++) {
        // Generate particle size
        double particle_size = size_dist(gen);
        particle_size = std::clamp(particle_size, size_range.first, size_range.second);
        double radius = particle_size / 2.0;
        
        // Generate 3D position with normal distribution for X,Y
        double x, y, r;
        int attempts = 0;
        do {
            x = x_dist(gen);
            y = y_dist(gen);
            r = std::sqrt(x * x + y * y);
            attempts++;
        } while (r > bed_radius && attempts < 100); // Keep within portafilter bounds
        
        double z = depth_dist(gen);
        
        // Size stratification: larger particles tend to settle deeper
        double size_factor = (particle_size - size_range.first) / (size_range.second - size_range.first);
        z *= (0.5 + 0.5 * size_factor); // Larger particles go deeper
        
        // Create particle
        CoffeeParticle3D particle(x, y, z, radius);
        
        // Simple overlap avoidance
        bool valid_position = true;
        for (const auto& existing : particles) {
            if (particle.distance_to(existing) < (radius + existing.get_size() / 2.0)) {
                valid_position = false;
                break;
            }
        }
        
        if (valid_position) {
            particles.push_back(particle);
        }
    }
    
    return particles;
}

void simulate_gravity_settling(
    std::vector<CoffeeParticle3D>& particles,
    double gravity_acceleration,
    double settling_time,
    double time_step) {
    
    int num_steps = static_cast<int>(settling_time / time_step);
    
    for (int step = 0; step < num_steps; step++) {
        for (auto& particle : particles) {
            // Apply gravity force
            double mass = particle.get_mass();
            double gravity_force = mass * gravity_acceleration;
            
            // Apply force in positive Z direction (downward)
            particle.apply_force_3d(0.0, 0.0, gravity_force);
        }
        
        // Simple collision resolution
        resolve_particle_overlaps(particles, 5);
    }
}

int resolve_particle_overlaps(
    std::vector<CoffeeParticle3D>& particles,
    int max_iterations) {
    
    int overlaps_resolved = 0;
    
    for (int iter = 0; iter < max_iterations; iter++) {
        bool found_overlap = false;
        
        for (size_t i = 0; i < particles.size(); i++) {
            for (size_t j = i + 1; j < particles.size(); j++) {
                if (particles[i].overlaps_with(particles[j])) {
                    found_overlap = true;
                    overlaps_resolved++;
                    
                    // Simple separation: move particles apart
                    auto [x1, y1, z1] = particles[i].get_position_3d();
                    auto [x2, y2, z2] = particles[j].get_position_3d();
                    
                    double dx = x2 - x1;
                    double dy = y2 - y1;
                    double dz = z2 - z1;
                    double distance = std::sqrt(dx * dx + dy * dy + dz * dz);
                    
                    if (distance > 0) {
                        double required_separation = (particles[i].get_size() + particles[j].get_size()) / 2.0;
                        double overlap = required_separation - distance;
                        
                        // Normalize direction
                        dx /= distance;
                        dy /= distance;
                        dz /= distance;
                        
                        // Move particles apart by half the overlap each
                        double move_distance = overlap / 2.0;
                        
                        particles[i].update_position_3d(-dx * move_distance, -dy * move_distance, -dz * move_distance);
                        particles[j].update_position_3d(dx * move_distance, dy * move_distance, dz * move_distance);
                    }
                }
            }
        }
        
        if (!found_overlap) break;
    }
    
    return overlaps_resolved;
}

double calculate_bed_porosity(
    const std::vector<CoffeeParticle3D>& particles,
    double bed_volume) {
    
    double total_particle_volume = 0.0;
    for (const auto& particle : particles) {
        total_particle_volume += particle.get_volume();
    }
    
    return 1.0 - (total_particle_volume / bed_volume);
}

void apply_size_stratification(
    std::vector<CoffeeParticle3D>& particles,
    double stratification_strength) {
    
    for (auto& particle : particles) {
        auto [x, y, z] = particle.get_position_3d();
        double size = particle.get_size();
        
        // Larger particles settle deeper
        double size_factor = (size - 200e-6) / (800e-6 - 200e-6); // Normalize to 0-1
        double depth_adjustment = stratification_strength * size_factor * 0.005; // Max 5mm adjustment
        
        particle.set_position_3d(x, y, z + depth_adjustment);
    }
}

} // namespace particle_3d_utils

} // namespace sprosim