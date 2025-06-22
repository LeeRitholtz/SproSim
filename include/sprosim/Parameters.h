#pragma once

namespace sprosim {

/**
 * @brief Physical and chemical parameters for brewing simulation
 */
struct Parameters {
    double permeability;             // bed permeability [m²]
    double fluid_viscosity;          // water viscosity [Pa·s]
    double extraction_rate;          // first-order extraction rate [1/s]
    double temperature;              // uniform temperature [K]
    double inlet_pressure;           // pressure at top of bed [Pa]
    double outlet_pressure;          // pressure at bottom of bed [Pa]
    double particle_drag;            // particle drag coefficient
    double flow_resistance;          // particle-induced flow resistance
    double saturation_concentration; // max concentration [kg/m³]
    double temperature_factor;       // Arrhenius temperature factor [K⁻¹]
};

} // namespace sprosim
