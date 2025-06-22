#pragma once
#include "sprosim/Parameters.h"
#include <memory>

namespace sprosim {

class ICoffeeBed;
class IWaterFlow;

/**
 * @brief Interface for fluid flow models in the coffee brewing simulation
 *
 * Defines the contract for updating the flow velocity field
 * based on the current state of the coffee bed, water flow grid,
 * and simulation parameters.
 */
class IFlowModel {
  public:
    virtual ~IFlowModel() = default;

    /**
     * @brief Update the velocity field for water flow through the coffee bed
     * @param water_flow Pointer to the water flow grid interface
     * @param bed Pointer to the coffee bed interface
     * @param params Simulation parameters containing physical constants and settings
     */
    virtual void update_velocity(std::shared_ptr<IWaterFlow> water_flow,
                                 std::shared_ptr<ICoffeeBed> bed, const Parameters& params) = 0;
};

} // namespace sprosim
