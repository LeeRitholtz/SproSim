#pragma once
#include "sprosim/CoffeeBed.h"
#include "sprosim/Parameters.h"
#include "sprosim/interfaces/IFlow.h"
#include "sprosim/interfaces/IFlowModel.h"
#include "sprosim/interfaces/IParticle.h"
#include "sprosim/interfaces/ITransportModel.h"

#include <memory>

namespace sprosim {

/**
 * @brief Core physics simulation engine for coffee brewing
 *
 * The PhysicsSolver integrates fluid flow, particle dynamics, and extraction
 * kinetics to simulate the coffee brewing process. It couples water flow
 * through the coffee bed with particle-level extraction and handles the
 * complex interactions between flow resistance, extraction rates, and
 * temperature effects.
 *
 * @example
 * ```cpp
 * PhysicsSolver::Parameters params{
 *     .permeability = 1e-12,     // m²
 *     .fluid_viscosity = 1e-3,   // Pa·s (water)
 *     .extraction_rate = 0.1,    // 1/s
 *     .temperature = 363.15      // K (90°C)
 * };
 * PhysicsSolver solver(bed, flow, params);
 * solver.simulate_step(0.01);  // 10ms timestep
 * ```
 */
class PhysicsSolver {
  public:
    PhysicsSolver(std::shared_ptr<CoffeeBed> bed, std::shared_ptr<IWaterFlow> flow,
                  Parameters params, std::shared_ptr<IFlowModel> flow_model,
                  std::shared_ptr<ITransportModel> transport_model);

    // Backwards-compatible constructor with default models
    PhysicsSolver(std::shared_ptr<CoffeeBed> bed, std::shared_ptr<IWaterFlow> flow,
                  Parameters params);

    virtual void simulate_step(double dt);

  protected:
    virtual void update_flow_field();
    virtual void handle_particle_interaction();
    virtual void update_extraction(double dt);

  private:
    void modify_local_flow(size_t i, size_t j, double particle_size, double extraction_state);
    void update_particle_state(const std::shared_ptr<ICoffeeParticle>& particle, double vx,
                               double vy);

    std::shared_ptr<CoffeeBed> coffee_bed_;
    std::shared_ptr<IWaterFlow> water_flow_;
    Parameters params_;

    std::shared_ptr<IFlowModel> flow_model_;
    std::shared_ptr<ITransportModel> transport_model_;
};

} // namespace sprosim
