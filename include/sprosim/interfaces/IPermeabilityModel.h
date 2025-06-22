#pragma once

namespace sprosim {

/**
 * @brief Interface for permeability models in the coffee brewing simulation
 *
 * Defines the contract for calculating effective permeability
 * based on factors like porosity or other physical parameters.
 */
class IPermeabilityModel {
  public:
    virtual ~IPermeabilityModel() = default;

    /**
     * @brief Calculate the effective permeability given the porosity of the bed
     * @param porosity Void fraction of the coffee bed (0.0 - 1.0)
     * @return Effective permeability value (in units consistent with simulation)
     */
    virtual double calculate_permeability(double porosity) const = 0;
};

} // namespace sprosim
