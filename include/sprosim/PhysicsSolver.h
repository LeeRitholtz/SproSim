#pragma once
#include <memory>
#include "sprosim/interfaces/IFlow.h"
#include "sprosim/interfaces/IParticle.h"
#include "sprosim/CoffeeBed.h"

namespace sprosim {

class PhysicsSolver {
public:
    struct Parameters {
        double permeability;           // bed permeability [m²]
        double fluid_viscosity;        // water viscosity [Pa·s]
        double extraction_rate;        // first-order extraction rate [1/s]
        double temperature;            // uniform temperature [K]
        double inlet_pressure;         // pressure at top of bed [Pa]
        double outlet_pressure;        // pressure at bottom of bed [Pa]
        double particle_drag;          // particle drag coefficient
        double flow_resistance;        // particle-induced flow resistance
        double saturation_concentration; // max concentration [kg/m³]
        double temperature_factor;     // Arrhenius temperature factor [K⁻¹]
    };

    PhysicsSolver(std::shared_ptr<CoffeeBed> bed,
                 std::shared_ptr<IWaterFlow> flow,
                 Parameters params);

    virtual void simulate_step(double dt);

protected:
    virtual void update_flow_field();
    virtual void handle_particle_interaction();
    virtual void update_extraction(double dt);

private:
    void modify_local_flow(size_t i, size_t j, double particle_size, double extraction_state);
    void update_particle_state(const std::shared_ptr<ICoffeeParticle>& particle, double vx, double vy);
    double calculate_effective_permeability() const;
    double calculate_temperature_factor() const;

    std::shared_ptr<CoffeeBed> coffee_bed_;
    std::shared_ptr<IWaterFlow> water_flow_;
    Parameters params_;
};

} // namespace sprosim