#pragma once
#include <sprosim/PhysicsSolver.h>
#include <utility>
#include <algorithm>

namespace sprosim {
namespace testing {

class MockCoffeeParticle : public ICoffeeParticle {
public:
    MockCoffeeParticle(double x, double y, double size) 
        : x_(x), y_(y), size_(size), extraction_state_(0.0), concentration_(0.0) {}
    
    std::pair<double, double> get_position() const override { return {x_, y_}; }
    double get_size() const override { return size_; }
    double get_extraction_state() const override { return extraction_state_; }
    double get_concentration() const override { return concentration_; }
    
    void update_extraction(double delta_conc, double dt) override {
        concentration_ += delta_conc;
        extraction_state_ = std::min(1.0, extraction_state_ + delta_conc / max_extractable_);
    }
    
    void apply_force(double fx, double fy) override {
        last_force_x_ = fx;
        last_force_y_ = fy;
    }

private:
    double x_, y_, size_;
    double extraction_state_;
    double concentration_;
    double last_force_x_ = 0.0, last_force_y_ = 0.0;
    const double max_extractable_ = 0.2;
};

}} // namespace sprosim::testing