#pragma once
#include "sprosim/CoffeeBed.h"
#include "sprosim/Parameters.h"
#include "sprosim/interfaces/IFlow.h"
#include "sprosim/interfaces/ITransportModel.h"
#include <cmath>
#include <memory>

namespace sprosim {

/**
 * @brief Linear extraction model implementing first-order extraction kinetics
 *
 * This model simulates coffee extraction using physically-based first-order kinetics
 * where the extraction rate is proportional to the remaining extractable material.
 *
 * Physics Implementation:
 * - **First-order kinetics**: extraction_rate = effective_rate * (1 - extraction_state)
 * - **Temperature dependence**: Arrhenius factor exp(temperature_factor * (T - 373.15))
 * - **Flow enhancement**: sqrt(1 + flow_magnitude) multiplier for mass transfer
 * - **Particle maximum**: 20% of particle mass is extractable (realistic for coffee)
 *
 * The model processes each particle individually, calculating local flow conditions
 * and updating both particle extraction state and water concentration grid.
 *
 * @note This model was corrected to use proper first-order kinetics instead of
 *       the previously incorrect concentration gradient approach.
 */
class LinearExtractionModel : public ITransportModel {
  public:
    /**
     * @brief Constructor
     */
    LinearExtractionModel() = default;

    /**
     * @brief Destructor
     */
    ~LinearExtractionModel() override = default;

    /**
     * @brief Update extraction physics for one time step
     *
     * For each particle in the bed:
     * 1. Maps particle position to flow grid coordinates
     * 2. Calculates local flow velocity and magnitude
     * 3. Applies temperature factor using Arrhenius equation
     * 4. Applies flow enhancement factor: sqrt(1 + flow_magnitude)
     * 5. Calculates extraction rate using first-order kinetics
     * 6. Updates particle concentration and extraction state
     * 7. Adds extracted material to water flow concentration grid
     *
     * @param water_flow Pointer to the water flow grid interface
     * @param bed Pointer to the coffee bed interface
     * @param params Simulation parameters (extraction_rate, temperature, temperature_factor)
     * @param dt Time step size [s]
     */
    void update_extraction(std::shared_ptr<IWaterFlow> water_flow, std::shared_ptr<CoffeeBed> bed,
                           const Parameters& params, double dt) override;

  private:
    /**
     * @brief Calculate temperature factor using Arrhenius equation
     *
     * Computes temperature enhancement factor using:
     * factor = exp(temperature_factor * (T - 373.15))
     *
     * Where 373.15 K is the reference temperature (100°C) and temperature_factor
     * controls sensitivity (typical values: 0.005-0.01 K^-1).
     *
     * @param params Simulation parameters containing temperature [K] and temperature_factor [K^-1]
     * @return Temperature enhancement factor (dimensionless, typically 0.5-2.0)
     */
    double calculate_temperature_factor(const Parameters& params) const;
};

} // namespace sprosim
