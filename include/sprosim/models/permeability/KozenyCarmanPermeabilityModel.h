#pragma once
#include "sprosim/interfaces/IPermeabilityModel.h"
#include <cmath>

namespace sprosim {

/**
 * @brief Kozeny-Carman permeability model implementation
 *
 * Provides the effective permeability calculation based on
 * the Kozeny-Carman empirical relation as a function of porosity.
 */
class KozenyCarmanPermeabilityModel : public IPermeabilityModel {
  public:
    /**
     * @param base_permeability Base permeability constant [m^2]
     */
    explicit KozenyCarmanPermeabilityModel(double base_permeability)
        : base_permeability_(base_permeability) {}

    /**
     * Calculate effective permeability using Kozeny-Carman formula
     * @param porosity Bed porosity (void fraction), range 0.0 to 1.0
     * @return Effective permeability [m^2]
     */
    double calculate_permeability(double porosity) const override {
        return base_permeability_ * std::pow(porosity, 3) / std::pow(1.0 - porosity, 2);
    }

  private:
    double base_permeability_;
};

} // namespace sprosim
