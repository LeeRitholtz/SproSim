#pragma once
#include <sprosim/PhysicsSolver.h>
#include <vector>
#include <utility>

namespace sprosim {
namespace testing {

class MockWaterFlow : public IWaterFlow {
public:
    MockWaterFlow(size_t nx, size_t ny) 
        : nx_(nx), ny_(ny), 
          velocities_(nx * ny * 2),
          concentrations_(nx * ny) {}
    
    std::pair<size_t, size_t> get_grid_dimensions() const override { return {nx_, ny_}; }
    double get_cell_size() const override { return 0.001; }
    
    void set_velocity(size_t i, size_t j, double vx, double vy) override {
        if (i >= nx_ || j >= ny_) return;
        size_t idx = 2 * (j * nx_ + i);
        velocities_[idx] = vx;
        velocities_[idx + 1] = vy;
    }
    
    std::pair<double, double> get_velocity(size_t i, size_t j) const override {
        if (i >= nx_ || j >= ny_) return {0.0, 0.0};
        size_t idx = 2 * (j * nx_ + i);
        return {velocities_[idx], velocities_[idx + 1]};
    }
    
    void add_concentration(size_t i, size_t j, double delta) override {
        if (i >= nx_ || j >= ny_) return;
        concentrations_[j * nx_ + i] += delta;
    }

private:
    size_t nx_, ny_;
    std::vector<double> velocities_;
    std::vector<double> concentrations_;
};

}} // namespace sprosim::testing