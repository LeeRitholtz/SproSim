#include "sprosim/PhysicsSolver.h"
#include "sprosim/interfaces/IFlowModel.h"
#include "sprosim/interfaces/IPermeabilityModel.h"
#include "sprosim/models/flow/DarcyFlowModel.h"
#include "sprosim/models/permeability/ConstantPermeabilityModel.h"
#include <cmath>

namespace sprosim {

PhysicsSolver::PhysicsSolver(std::shared_ptr<CoffeeBed> bed, std::shared_ptr<IWaterFlow> flow,
                             Parameters params, std::shared_ptr<IPermeabilityModel> perm_model,
                             std::shared_ptr<IFlowModel> flow_model)
    : coffee_bed_(bed), water_flow_(flow), params_(params), permeability_model_(perm_model),
      flow_model_(flow_model) {}

// Backwards-compatible constructor with default models
PhysicsSolver::PhysicsSolver(std::shared_ptr<CoffeeBed> bed, std::shared_ptr<IWaterFlow> flow,
                             Parameters params)
    : coffee_bed_(bed), water_flow_(flow), params_(params),
      permeability_model_(std::make_shared<ConstantPermeabilityModel>(params.permeability)),
      flow_model_(std::make_shared<DarcyFlowModel>(
          std::make_shared<ConstantPermeabilityModel>(params.permeability))) {}

void PhysicsSolver::simulate_step(double dt) {
    update_flow_field();
    handle_particle_interaction();
    update_extraction(dt);
}

void PhysicsSolver::update_flow_field() {
    flow_model_->update_velocity(water_flow_, coffee_bed_, params_);
}

double PhysicsSolver::calculate_temperature_factor() const {
    return std::exp(params_.temperature_factor * (params_.temperature - 373.15));
}

void PhysicsSolver::handle_particle_interaction() {
    for (const auto& particle : coffee_bed_->get_particles()) {
        auto [px, py] = particle->get_position();
        const auto [nx, ny] = water_flow_->get_grid_dimensions();
        const double dx = water_flow_->get_cell_size();

        size_t i = static_cast<size_t>(px / dx);
        size_t j = static_cast<size_t>(py / dx);

        if (i >= nx || j >= ny)
            continue;

        auto [vx, vy] = water_flow_->get_velocity(i, j);
        double particle_size = particle->get_size();
        double extraction_state = particle->get_extraction_state();

        modify_local_flow(i, j, particle_size, extraction_state);
        update_particle_state(particle, vx, vy);
    }

    coffee_bed_->update_compaction(params_.inlet_pressure);
}

void PhysicsSolver::update_extraction(double dt) {
    double temp_factor = calculate_temperature_factor();

    for (const auto& particle : coffee_bed_->get_particles()) {
        auto [px, py] = particle->get_position();
        const auto [nx, ny] = water_flow_->get_grid_dimensions();
        const double dx = water_flow_->get_cell_size();

        size_t i = static_cast<size_t>(px / dx);
        size_t j = static_cast<size_t>(py / dx);

        if (i >= nx || j >= ny)
            continue;

        auto [vx, vy] = water_flow_->get_velocity(i, j);
        double flow_magnitude = std::sqrt(vx * vx + vy * vy);

        double current_concentration = particle->get_concentration();
        double extraction_state = particle->get_extraction_state();

        double effective_rate = params_.extraction_rate * temp_factor;
        double flow_enhancement = std::sqrt(1.0 + flow_magnitude);
        effective_rate *= flow_enhancement;

        double concentration_gradient =
            params_.saturation_concentration * (1.0 - extraction_state) - current_concentration;

        double delta_concentration = effective_rate * concentration_gradient * dt;

        particle->update_extraction(delta_concentration, dt);
        water_flow_->add_concentration(i, j, delta_concentration);
    }
}

void PhysicsSolver::modify_local_flow(size_t i, size_t j, double particle_size,
                                      double extraction_state) {
    auto [vx, vy] = water_flow_->get_velocity(i, j);
    double resistance_factor = params_.flow_resistance * particle_size * (1.0 - extraction_state);
    vy *= (1.0 - resistance_factor);
    water_flow_->set_velocity(i, j, vx, vy);
}

void PhysicsSolver::update_particle_state(const std::shared_ptr<ICoffeeParticle>& particle,
                                          double vx, double vy) {
    double flow_magnitude = std::sqrt(vx * vx + vy * vy);
    particle->apply_force(0, params_.particle_drag * flow_magnitude);
}

} // namespace sprosim
