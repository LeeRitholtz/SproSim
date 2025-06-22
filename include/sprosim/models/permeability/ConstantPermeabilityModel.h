#pragma once
#include "sprosim/interfaces/IPermeabilityModel.h"

namespace sprosim {

/**
 * @brief Constant permeability model implementation
 *
 * Returns a fixed permeability value regardless of porosity.
 */
class ConstantPermeabilityModel : public IPermeabilityModel {
  public:
    explicit ConstantPermeabilityModel(double permeability) : permeability_(permeability) {}

    double calculate_permeability(double /*porosity*/) const override {
        return permeability_;
    }

  private:
    double permeability_;
};

} // namespace sprosim
