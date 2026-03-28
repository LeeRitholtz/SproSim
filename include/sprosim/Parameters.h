#pragma once

namespace sprosim {

/**
 * @brief Physical and chemical parameters for brewing simulation
 *
 * All units use human-friendly conventions:
 * - Temperature in degrees Celsius
 * - Pressure in bar
 * - Permeability in m²
 * - Viscosity in Pa·s
 * - Concentration in kg/m³
 */
struct Parameters {
    double permeability;             // bed permeability [m²]
    double fluid_viscosity;          // water viscosity [Pa·s]
    double extraction_rate;          // first-order extraction rate [1/s]
    double temperature;              // uniform temperature [°C]
    double inlet_pressure;           // pressure at top of bed [bar]
    double outlet_pressure;          // pressure at bottom of bed [bar]
    double particle_drag;            // particle drag coefficient
    double flow_resistance;          // particle-induced flow resistance
    double saturation_concentration; // max concentration [kg/m³]
    double temperature_factor;       // Arrhenius temperature factor [°C⁻¹]
};

} // namespace sprosim