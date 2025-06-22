#include "sprosim/Parameters.h"
#include "sprosim/interfaces/ICoffeeBed.h"
#include "sprosim/interfaces/IFlow.h"
#include "sprosim/interfaces/IParticle.h"
#include "sprosim/models/flow/DarcyFlowModel.h"
#include "sprosim/models/permeability/ConstantPermeabilityModel.h"
#include "sprosim/models/permeability/KozenyCarmanPermeabilityModel.h"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <memory>

using namespace sprosim;

// Mock IWaterFlow implementation for testing
class MockWaterFlow : public IWaterFlow {
  public:
    MockWaterFlow(size_t width, size_t height, double cell_size)
        : width_(width), height_(height), cell_size_(cell_size),
          velocities_(width * height, {0.0, 0.0}) {}

    std::pair<size_t, size_t> get_grid_dimensions() const override {
        return {width_, height_};
    }

    double get_cell_size() const override {
        return cell_size_;
    }

    void set_velocity(size_t i, size_t j, double vx, double vy) override {
        if (i >= width_ || j >= height_)
            return;
        velocities_[j * width_ + i] = {vx, vy};
    }

    std::pair<double, double> get_velocity(size_t i, size_t j) const override {
        if (i >= width_ || j >= height_)
            return {0.0, 0.0};
        return velocities_[j * width_ + i];
    }

    void add_concentration(size_t, size_t, double) override {
        // no-op for test
    }

    void update_pressure(double, double) override {
        // no-op for test
    }

  private:
    size_t width_, height_;
    double cell_size_;
    std::vector<std::pair<double, double>> velocities_;
};

// Minimal CoffeeParticle mock implementing ICoffeeParticle for bed height and porosity tests if
// needed
class MockCoffeeParticle : public ICoffeeParticle {
  public:
    MockCoffeeParticle(double x, double y) : pos_{x, y} {}
    std::pair<double, double> get_position() const {
        return pos_;
    }
    double get_size() const {
        return 0.001;
    }
    double get_extraction_state() const {
        return 0.0;
    }
    double get_concentration() const {
        return 0.0;
    }
    void update_extraction(double, double) {}
    void apply_force(double, double) {}

  private:
    std::pair<double, double> pos_;
};

class MockCoffeeBed : public ICoffeeBed {
  public:
    MockCoffeeBed() : porosity_(0.4), bed_height_(0.058) {}

    double get_porosity() const override {
        return porosity_;
    }

    double get_bed_height() const override {
        return bed_height_;
    }

    void set_porosity(double p) {
        porosity_ = p;
    }

    void set_bed_height(double h) {
        bed_height_ = h;
    }

  private:
    double porosity_;
    double bed_height_;
};

TEST_CASE("DarcyFlowModel updates velocity field using permeability model", "[DarcyFlowModel]") {
    // Arrange
    size_t width = 5;
    size_t height = 5;
    double cell_size = 0.001; // 1mm cell

    auto water_flow = std::make_shared<MockWaterFlow>(width, height, cell_size);
    auto coffee_bed = std::make_shared<MockCoffeeBed>();

    Parameters params;
    params.fluid_viscosity = 1e-3; // water viscosity Pa.s
    params.inlet_pressure = 2e5;   // 2 bar
    params.outlet_pressure = 1e5;  // 1 bar

    // Use constant permeability for simplicity, 1e-12 m2
    auto perm_model = std::make_shared<ConstantPermeabilityModel>(1e-12);
    auto flow_model = std::make_shared<DarcyFlowModel>(perm_model);

    // Act
    flow_model->update_velocity(water_flow, coffee_bed, params);

    // Calculate expected velocity manually
    double bed_height = 0.058; // basket height ~0.058 m from CoffeeBed default
    double dp_dy = (params.outlet_pressure - params.inlet_pressure) / bed_height;
    double expected_velocity =
        -(perm_model->calculate_permeability(coffee_bed->get_porosity()) / params.fluid_viscosity) *
        dp_dy;

    // Assert
    for (size_t i = 0; i < width; ++i) {
        for (size_t j = 0; j < height; ++j) {
            auto [vx, vy] = water_flow->get_velocity(i, j);
            INFO("Grid cell (" << i << "," << j << ") velocity: vx=" << vx << ", vy=" << vy);
            REQUIRE_THAT(vx, Catch::Matchers::WithinAbs(0.0, 1e-10));
            REQUIRE_THAT(vy, Catch::Matchers::WithinAbs(expected_velocity, 1e-5));
        }
    }
}

TEST_CASE("DarcyFlowModel velocity responds to porosity changes via permeability model",
          "[DarcyFlowModel]") {
    // Arrange
    size_t width = 3;
    size_t height = 3;
    double cell_size = 0.001;

    auto water_flow = std::make_shared<MockWaterFlow>(width, height, cell_size);
    auto coffee_bed = std::make_shared<MockCoffeeBed>();

    Parameters params;
    params.fluid_viscosity = 1e-3;
    params.inlet_pressure = 200000;
    params.outlet_pressure = 100000;

    auto kc_model = std::make_shared<KozenyCarmanPermeabilityModel>(1e-12);
    auto flow_model = std::make_shared<DarcyFlowModel>(kc_model);

    // Test low porosity
    coffee_bed->set_porosity(0.3);
    flow_model->update_velocity(water_flow, coffee_bed, params);
    double dp_dy = (params.outlet_pressure - params.inlet_pressure) / coffee_bed->get_bed_height();
    double expected_velocity_low =
        -(kc_model->calculate_permeability(coffee_bed->get_porosity()) / params.fluid_viscosity) *
        dp_dy;

    for (size_t i = 0; i < width; ++i) {
        for (size_t j = 0; j < height; ++j) {
            auto [vx, vy] = water_flow->get_velocity(i, j);
            INFO("Grid cell (" << i << "," << j << ") velocity: vx=" << vx << ", vy=" << vy);
            REQUIRE_THAT(vx, Catch::Matchers::WithinAbs(0.0, 1e-10));
            REQUIRE_THAT(vy, Catch::Matchers::WithinAbs(expected_velocity_low, 1e-5));
        }
    }

    // Test higher porosity
    coffee_bed->set_porosity(0.5);
    flow_model->update_velocity(water_flow, coffee_bed, params);
    double expected_velocity_high =
        -(kc_model->calculate_permeability(coffee_bed->get_porosity()) / params.fluid_viscosity) *
        dp_dy;

    for (size_t i = 0; i < width; ++i) {
        for (size_t j = 0; j < height; ++j) {
            auto [vx, vy] = water_flow->get_velocity(i, j);
            INFO("Grid cell (" << i << "," << j << ") velocity: vx=" << vx << ", vy=" << vy);
            REQUIRE_THAT(vx, Catch::Matchers::WithinAbs(0.0, 1e-10));
            REQUIRE_THAT(vy, Catch::Matchers::WithinAbs(expected_velocity_high, 1e-5));
        }
    }
}
