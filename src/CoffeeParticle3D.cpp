#include "sprosim/CoffeeParticle3D.h"
#include <algorithm>
#include <cmath>
#include <random>

namespace sprosim {

CoffeeParticle3D::CoffeeParticle3D(double x, double y, double z, double radius)
    : position_{x, y, z},
      radius_(radius),
      extraction_state_(0.0),
      concentration_(0.0),
      velocity_{0.0, 0.0, 0.0}
{}

std::pair<double, double> CoffeeParticle3D::get_position() const {
    return {position_[0], position_[1]};
}

std::tuple<double, double, double> CoffeeParticle3D::get_position_3d() const {
    return {position_[0], position_[1], position_[2]};
}

double CoffeeParticle3D::get_depth() const {
    return position_[2];
}

double CoffeeParticle3D::get_size() const {
    return radius_ * 2.0;
}

double CoffeeParticle3D::get_extraction_state() const {
    return extraction_state_;
}

double CoffeeParticle3D::get_concentration() const {
    return concentration_;
}

void CoffeeParticle3D::update_extraction(double delta_conc, double dt) {
    concentration_ += delta_conc;
    extraction_state_ = std::min(1.0, extraction_state_ + delta_conc / max_extractable_);
}

void CoffeeParticle3D::apply_force(double fx, double fy) {
    apply_force_3d(fx, fy, 0.0);
}

void CoffeeParticle3D::apply_force_3d(double fx, double fy, double fz) {
    const double dt = 0.001;  // 1ms timestep
    double mass = calculate_mass();
    
    // F = ma -> a = F/m
    double ax = fx / mass;
    double ay = fy / mass;
    double az = fz / mass;
    
    // Update velocity
    velocity_[0] += ax * dt;
    velocity_[1] += ay * dt;
    velocity_[2] += az * dt;
    
    // Update position using Verlet integration
    position_[0] += velocity_[0] * dt + 0.5 * ax * dt * dt;
    position_[1] += velocity_[1] * dt + 0.5 * ay * dt * dt;
    position_[2] += velocity_[2] * dt + 0.5 * az * dt * dt;
    
    // Apply boundary constraints
    apply_boundary_constraints();
}

void CoffeeParticle3D::update_position_3d(double dx, double dy, double dz) {
    position_[0] += dx;
    position_[1] += dy;
    position_[2] += dz;
    apply_boundary_constraints();
}

void CoffeeParticle3D::set_position_3d(double x, double y, double z) {
    position_[0] = x;
    position_[1] = y;
    position_[2] = z;
    apply_boundary_constraints();
}

double CoffeeParticle3D::get_volume() const {
    return (4.0 / 3.0) * pi_ * radius_ * radius_ * radius_;
}

double CoffeeParticle3D::distance_to(const CoffeeParticle3D& other) const {
    auto [ox, oy, oz] = other.get_position_3d();
    double dx = position_[0] - ox;
    double dy = position_[1] - oy;
    double dz = position_[2] - oz;
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

bool CoffeeParticle3D::overlaps_with(const CoffeeParticle3D& other) const {
    double distance = distance_to(other);
    double combined_radius = radius_ + other.radius_;
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

double CoffeeParticle3D::calculate_mass() const {
    return particle_density_ * get_volume();
}

void CoffeeParticle3D::apply_boundary_constraints() {
    // Ensure particle stays within reasonable bounds
    // Z should be non-negative (above bed bottom)
    position_[2] = std::max(0.0, position_[2]);
    
    // Limit maximum depth (prevent particles from falling through bed)
    position_[2] = std::min(0.02, position_[2]); // 20mm max depth
    
    // Radial constraint (stay within portafilter)
    double r = std::sqrt(position_[0] * position_[0] + position_[1] * position_[1]);
    const double max_radius = 0.029; // 58mm diameter / 2
    if (r > max_radius) {
        double scale = max_radius / r;
        position_[0] *= scale;
        position_[1] *= scale;
    }
}

void CoffeeParticle3D::update_physics_state(double dt) {
    // Apply damping to velocity (friction)
    const double damping = 0.95;
    velocity_[0] *= damping;
    velocity_[1] *= damping;
    velocity_[2] *= damping;
}

namespace particle_3d_utils {

std::vector<CoffeeParticle3D> generate_realistic_distribution(
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
    
    // Position distributions
    std::uniform_real_distribution<double> angle_dist(0.0, 2.0 * M_PI);
    std::uniform_real_distribution<double> radius_dist(0.0, 1.0);
    std::uniform_real_distribution<double> depth_dist(0.0, bed_depth);
    
    for (int i = 0; i < num_particles; i++) {
        // Generate particle size
        double particle_size = size_dist(gen);
        particle_size = std::clamp(particle_size, size_range.first, size_range.second);
        double radius = particle_size / 2.0;
        
        // Generate 3D position
        double r_fraction = std::sqrt(radius_dist(gen)); // Square root for uniform area distribution
        double r = bed_radius * r_fraction;
        double theta = angle_dist(gen);
        double z = depth_dist(gen);
        
        double x = r * std::cos(theta);
        double y = r * std::sin(theta);
        
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
            double mass = particle.calculate_mass();
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