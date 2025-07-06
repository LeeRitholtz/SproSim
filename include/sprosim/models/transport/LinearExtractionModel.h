#pragma once
#include "sprosim/CoffeeBed.h"
#include "sprosim/Parameters.h"
#include "sprosim/interfaces/IFlow.h"
#include "sprosim/interfaces/ITransportModel.h"
#include <cmath>
#include <memory>

namespace sprosim {

/**
 * @brief Linear extraction model with flow enhancement and temperature dependence
 *
 * Implements first-order extraction kinetics with:
 * - Arrhenius temperature dependence
 * - Flow-enhanced mass transfer (sqrt(1 + flow_magnitude))
 * - Linear concentration gradient driving force
 * - Saturation concentration limiting
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
     * @param water_flow Pointer to the water flow grid interface
     * @param bed Pointer to the coffee bed interface
     * @param params Simulation parameters containing extraction rates and limits
     * @param dt Time step size [s]
     */
    void update_extraction(std::shared_ptr<IWaterFlow> water_flow, std::shared_ptr<CoffeeBed> bed,
                           const Parameters& params, double dt) override;

  private:
    /**
     * @brief Calculate temperature factor using Arrhenius equation
     * @param params Simulation parameters containing temperature and temperature_factor
     * @return Temperature enhancement factor (dimensionless)
     */
    double calculate_temperature_factor(const Parameters& params) const;
};

} // namespace sprosim
