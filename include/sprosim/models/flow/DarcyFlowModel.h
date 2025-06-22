#pragma once
#include "sprosim/Parameters.h"
#include "sprosim/interfaces/ICoffeeBed.h"
#include "sprosim/interfaces/IFlow.h"
#include "sprosim/interfaces/IFlowModel.h"
#include "sprosim/interfaces/IPermeabilityModel.h"
#include <memory>

namespace sprosim {

/**
 * @brief Darcy flow model implementation
 *
 * Calculates water velocity through the coffee bed using Darcy's law,
 * utilizing an injected permeability model to obtain effective permeability.
 */
class DarcyFlowModel : public IFlowModel {
  public:
    explicit DarcyFlowModel(std::shared_ptr<IPermeabilityModel> permeability_model)
        : permeability_model_(permeability_model) {}

    void update_velocity(std::shared_ptr<IWaterFlow> water_flow, std::shared_ptr<ICoffeeBed> bed,
                         const Parameters& params) override {

        auto [nx, ny] = water_flow->get_grid_dimensions();
        double bed_height = bed->get_bed_height();
        if (bed_height <= 0.0)
            return;

        double delta_p = params.outlet_pressure - params.inlet_pressure;
        double dp_dy = delta_p / bed_height;

        double porosity = bed->get_porosity();
        double effective_k = permeability_model_->calculate_permeability(porosity);

        double velocity = -(effective_k / params.fluid_viscosity) * dp_dy;

        for (size_t i = 0; i < nx; ++i) {
            for (size_t j = 0; j < ny; ++j) {
                water_flow->set_velocity(i, j, 0.0, velocity);
            }
        }
    }

  private:
    std::shared_ptr<IPermeabilityModel> permeability_model_;
};

} // namespace sprosim
