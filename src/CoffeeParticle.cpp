#include "sprosim/CoffeeParticle.h"
#include <algorithm>
// #include <cmath>

namespace sprosim {

CoffeeParticle::CoffeeParticle(double x, double y, double radius)
    : position_{x, y}, radius_(radius), extraction_state_(0.0), concentration_(0.0) {}

std::pair<double, double> CoffeeParticle::get_position() const {
    return {position_[0], position_[1]};
}

double CoffeeParticle::get_size() const {
    return radius_ * 2.0;
}

double CoffeeParticle::get_extraction_state() const {
    return extraction_state_;
}

double CoffeeParticle::get_concentration() const {
    return concentration_;
}

void CoffeeParticle::update_extraction(double delta_conc, double dt) {
    concentration_ += delta_conc;
    extraction_state_ = std::min(1.0, extraction_state_ + delta_conc / max_extractable_);
}

void CoffeeParticle::apply_force(double fx, double fy) {
    // Simple force application - will be expanded later
    // TODO: Develop more complex force physics
    const double dt = 0.001;  // 1ms timestep
    const double mass = 1e-6; // 1mg particle mass

    // F = ma -> a = F/m
    double ax = fx / mass;
    double ay = fy / mass;

    // Update position using simple Euler integration
    position_[0] += 0.5 * ax * dt * dt;
    position_[1] += 0.5 * ay * dt * dt;
}

void CoffeeParticle::update_position(double dx, double dy) {
    position_[0] += dx;
    position_[1] += dy;
}

} // namespace sprosim