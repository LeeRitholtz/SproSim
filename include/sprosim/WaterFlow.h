#pragma once
#include "sprosim/interfaces/IFlow.h"
#include <vector>

namespace sprosim {

class WaterFlow : public IWaterFlow {
public:
    WaterFlow(size_t width, size_t height);
    
    // IFlow interface
    std::pair<size_t, size_t> get_grid_dimensions() const override;
    double get_cell_size() const override;
    
    // IWaterFlow interface
    void set_velocity(size_t i, size_t j, double vx, double vy) override;
    std::pair<double, double> get_velocity(size_t i, size_t j) const override;
    void add_concentration(size_t i, size_t j, double delta) override;
    
    // Additional methods
    void update_pressure(double dt, double input_pressure);
    double pressure_at(size_t x, size_t y) const;

private:
    size_t width_, height_;
    std::vector<double> pressure_;
    std::vector<double> velocity_x_;
    std::vector<double> velocity_y_;
    std::vector<double> concentration_;
    const double cell_size_ = 0.001; // 1mm cells
};

} // namespace sprosim