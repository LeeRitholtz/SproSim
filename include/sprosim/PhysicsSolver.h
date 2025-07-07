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
 * Uses pluggable transport models (ITransportModel) for customizable extraction
 * physics, allowing different extraction kinetics implementations.
 *
 * @example
 * ```cpp
 * PhysicsSolver::Parameters params{
 *     .permeability = 1e-12,     // m²
 *     .fluid_viscosity = 1e-3,   // Pa·s (water)
 *     .extraction_rate = 0.1,    // 1/s
 *     .temperature = 363.15      // K (90°C)
 * };
 * // Using backwards-compatible constructor with default models
 * PhysicsSolver solver(bed, flow, params);
 * solver.simulate_step(0.01);  // 10ms timestep
 * ```
 */
class PhysicsSolver {
  public:
    /**
     * @brief Constructor with explicit flow and transport models
     * @param bed Coffee bed containing particles to simulate
     * @param flow Water flow field for fluid dynamics
     * @param params Physics parameters (permeability, temperature, etc.)
     * @param flow_model Flow physics model (e.g., DarcyFlowModel)
     * @param transport_model Extraction physics model (e.g., LinearExtractionModel)
     */
    PhysicsSolver(std::shared_ptr<CoffeeBed> bed, std::shared_ptr<IWaterFlow> flow,
                  Parameters params, std::shared_ptr<IFlowModel> flow_model,
                  std::shared_ptr<ITransportModel> transport_model);

    /**
     * @brief Backwards-compatible constructor with default models
     * @param bed Coffee bed containing particles to simulate
     * @param flow Water flow field for fluid dynamics
     * @param params Physics parameters (permeability, temperature, etc.)
     *
     * Uses default DarcyFlowModel and LinearExtractionModel automatically.
     */
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
