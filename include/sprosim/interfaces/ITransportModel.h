#pragma once
#include "sprosim/CoffeeBed.h"
#include "sprosim/Parameters.h"
#include "sprosim/interfaces/IFlow.h"
#include <memory>

namespace sprosim {

/**
 * @brief Interface for solute transport models in the coffee brewing simulation
 *
 * Defines the contract for updating extraction physics, including:
 * - Mass transfer from coffee particles to water flow
 * - Cup yield calculation (water + dissolved solids)
 * - Temperature-dependent extraction kinetics
 */
class ITransportModel {
  public:
    virtual ~ITransportModel() = default;

    /**
     * @brief Update extraction physics for one time step
     * @param water_flow Pointer to the water flow grid interface
     * @param bed Pointer to the coffee bed interface
     * @param params Simulation parameters containing extraction rates and limits
     * @param dt Time step size [s]
     */
    virtual void update_extraction(std::shared_ptr<IWaterFlow> water_flow,
                                   std::shared_ptr<CoffeeBed> bed, const Parameters& params,
                                   double dt) = 0;
};

} // namespace sprosim
