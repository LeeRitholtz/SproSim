#include "sprosim/WaterFlow.h"

namespace sprosim {

WaterFlow::WaterFlow(size_t width, size_t height):
    width_(width),
    height_(height),
    pressure_(width * height, 0.0),
    velocity_x_(width * height, 0.0),
    velocity_y_(width * height, 0.0),
    concentration_(width * height, 0.0)
{}

std::pair<size_t, size_t> WaterFlow::get_grid_dimensions() const {
    return {width_, height_};
}

double WaterFlow::get_cell_size() const {
    return cell_size_;
}

void WaterFlow::set_velocity(size_t i, size_t j, double vx, double vy) {
    if (i >= width_ || j >= height_) return;
    size_t idx = j * width_ + i;
    velocity_x_[idx] = vx;
    velocity_y_[idx] = vy;
}

std::pair<double, double> WaterFlow::get_velocity(size_t i, size_t j) const {
    if (i >= width_ || j >= height_) return {0.0, 0.0};
    size_t idx = j * width_ + i;
    return {velocity_x_[idx], velocity_y_[idx]};
}

void WaterFlow::add_concentration(size_t i, size_t j, double delta) {
    if (i >= width_ || j >= height_) return;
    concentration_[j * width_ + i] += delta;
}

void WaterFlow::update_pressure(double dt, double input_pressure) {
    for (size_t y = 0; y < height_; ++y) {
        for (size_t x = 0; x < width_; ++x) {
            pressure_[y * width_ + x] = input_pressure * (1.0 - static_cast<double>(y) / height_);
        }
    }
}

double WaterFlow::pressure_at(size_t x, size_t y) const {
    if (x >= width_ || y >= height_) return 0.0;
    return pressure_[y * width_ + x];
}

} // namespace sprosim