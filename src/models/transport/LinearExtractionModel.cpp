#include "sprosim/models/transport/LinearExtractionModel.h"
#include <cmath>

namespace sprosim {

void LinearExtractionModel::update_extraction(std::shared_ptr<IWaterFlow> water_flow,
                                              std::shared_ptr<CoffeeBed> bed,
                                              const Parameters& params, double dt) {
    double temp_factor = calculate_temperature_factor(params);

    for (const auto& particle : bed->get_particles()) {
        auto [px, py] = particle->get_position();
        const auto [nx, ny] = water_flow->get_grid_dimensions();
        const double dx = water_flow->get_cell_size();

        size_t i = static_cast<size_t>(px / dx);
        size_t j = static_cast<size_t>(py / dx);

        if (i >= nx || j >= ny)
            continue;

        auto [vx, vy] = water_flow->get_velocity(i, j);
        double flow_magnitude = std::sqrt(vx * vx + vy * vy);

        double current_concentration = particle->get_concentration();
        double extraction_state = particle->get_extraction_state();

        double effective_rate = params.extraction_rate * temp_factor;
        double flow_enhancement = std::sqrt(1.0 + flow_magnitude);
        effective_rate *= flow_enhancement;

        // First-order kinetics: Rate proportional to remaining extractable material
        double remaining_extractable = (1.0 - extraction_state);

        // Simple extraction rate with temperature and flow enhancement
        double extraction_rate = effective_rate * remaining_extractable;

        // Convert to concentration units (particle max_extractable = 0.2)
        double particle_max_extractable = 0.2;
        double delta_concentration = extraction_rate * particle_max_extractable * dt;

        particle->update_extraction(delta_concentration, dt);
        water_flow->add_concentration(i, j, delta_concentration);
    }
}

double LinearExtractionModel::calculate_temperature_factor(const Parameters& params) const {
    return std::exp(params.temperature_factor * (params.temperature - 373.15));
}

} // namespace sprosim
