#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "sprosim/CoffeeParticle3D.h"
#include <cmath>
#include <numeric>

using namespace sprosim;
using namespace sprosim::particle_3d_utils;

TEST_CASE("generate_random_3d_distribution creates particles", "[particle3d]") {
    const int num_particles = 100;
    const double bed_radius = 0.02;  // 2 cm
    const double bed_depth = 0.01;   // 1 cm
    const std::pair<double, double> size_range = {200e-6, 800e-6};
    
    auto particles = generate_random_3d_distribution(num_particles, bed_radius, bed_depth, size_range);
    
    // Due to rejection sampling, we might get slightly fewer particles than requested
    REQUIRE(particles.size() >= static_cast<size_t>(num_particles * 0.95));
    REQUIRE(particles.size() <= static_cast<size_t>(num_particles));
}

TEST_CASE("generate_random_3d_distribution respects bounds", "[particle3d]") {
    const int num_particles = 50;
    const double bed_radius = 0.015;
    const double bed_depth = 0.008;
    const std::pair<double, double> size_range = {300e-6, 600e-6};
    
    auto particles = generate_random_3d_distribution(num_particles, bed_radius, bed_depth, size_range);
    
    for (const auto& particle : particles) {
        auto [x, y, z] = particle.get_position_3d();
        
        // Check radial bounds
        double r = std::sqrt(x*x + y*y);
        REQUIRE(r <= bed_radius);
        
        // Check depth bounds
        REQUIRE(z >= 0.0);
        REQUIRE(z <= bed_depth);
        
        // Check size bounds
        double size = particle.get_size();
        REQUIRE(size >= size_range.first);
        REQUIRE(size <= size_range.second);
    }
}

TEST_CASE("generate_random_3d_distribution uses normal distribution for X,Y", "[particle3d]") {
    const int num_particles = 200;
    const double bed_radius = 0.02;
    const double bed_depth = 0.01;
    
    auto particles = generate_random_3d_distribution(num_particles, bed_radius, bed_depth);
    
    // Collect X and Y coordinates
    std::vector<double> x_coords, y_coords;
    for (const auto& particle : particles) {
        auto [x, y, z] = particle.get_position_3d();
        x_coords.push_back(x);
        y_coords.push_back(y);
    }
    
    // Check that mean is close to zero (normal distribution centered at 0)
    double x_mean = std::accumulate(x_coords.begin(), x_coords.end(), 0.0) / num_particles;
    double y_mean = std::accumulate(y_coords.begin(), y_coords.end(), 0.0) / num_particles;
    
    REQUIRE_THAT(x_mean, Catch::Matchers::WithinAbs(0.0, 0.003));
    REQUIRE_THAT(y_mean, Catch::Matchers::WithinAbs(0.0, 0.003));
    
    // Check that most particles are within 2 standard deviations (expected ~95%)
    double expected_std = bed_radius * 0.3;  // From implementation
    int within_2std_x = 0, within_2std_y = 0;
    
    for (double x : x_coords) {
        if (std::abs(x) <= 2 * expected_std) within_2std_x++;
    }
    for (double y : y_coords) {
        if (std::abs(y) <= 2 * expected_std) within_2std_y++;
    }
    
    // At least 90% should be within 2 standard deviations
    REQUIRE(within_2std_x >= static_cast<int>(0.90 * num_particles));
    REQUIRE(within_2std_y >= static_cast<int>(0.90 * num_particles));
}