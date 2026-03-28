/// @file custom_transport_model.cpp
/// @brief Demonstrates implementing a custom ITransportModel and using it with PhysicsSolver.

#include <cmath>
#include <iostream>
#include <memory>

#include "sprosim/CoffeeBed.h"
#include "sprosim/CoffeeParticle.h"
#include "sprosim/PhysicsSolver.h"
#include "sprosim/WaterFlow.h"
#include "sprosim/interfaces/ITransportModel.h"
#include "sprosim/models/flow/DarcyFlowModel.h"
#include "sprosim/models/permeability/KozenyCarmanPermeabilityModel.h"

using namespace sprosim;

/// A custom transport model that uses exponential decay extraction kinetics.
/// Extraction slows down exponentially as the particle approaches full extraction,
/// with a configurable rate constant and temperature sensitivity.
class ExponentialExtractionModel : public ITransportModel {
  public:
    ExponentialExtractionModel(double rate_constant, double reference_temp_celsius)
        : rate_constant_(rate_constant), reference_temp_(reference_temp_celsius) {}

    void update_extraction(std::shared_ptr<IWaterFlow> water_flow,
                           std::shared_ptr<CoffeeBed> bed, const Parameters& params,
                           double dt) override {
        // Temperature factor: extraction speeds up ~5% per degree above reference
        double temp_factor = 1.0 + 0.05 * (params.temperature - reference_temp_);

        for (auto& particle : bed->get_particles()) {
            double state = particle->get_extraction_state();
            double remaining = 1.0 - state;

            // Exponential decay: extraction rate proportional to remaining extractable material
            double delta_conc = rate_constant_ * temp_factor * remaining * dt;
            delta_conc = std::min(delta_conc, remaining);

            particle->update_extraction(delta_conc, dt);
        }
    }

  private:
    double rate_constant_;
    double reference_temp_;
};

int main() {
    // Set up physics parameters
    Parameters params{};
    params.permeability = 1e-12;
    params.fluid_viscosity = 0.001;
    params.extraction_rate = 0.01;
    params.temperature = 95.0;                // °C
    params.inlet_pressure = 10.01325;         // 9 bar + atmospheric [bar]
    params.outlet_pressure = 1.01325;         // atmospheric [bar]
    params.particle_drag = 0.1;
    params.flow_resistance = 0.5;
    params.saturation_concentration = 0.30;
    params.temperature_factor = 0.01;

    // Create simulation components
    auto bed = std::make_shared<CoffeeBed>(18.0, 58.0);  // 18g dose, 58mm basket

    // Add particles on the grid (58x30 cells at 1mm each = 0.058m x 0.030m)
    for (int i = 0; i < 350; ++i) {
        double x = 0.001 * static_cast<double>(i % 58);         // 0.0 to 0.057m
        double y = 0.001 * static_cast<double>(i / 58);         // 0.0 to ~0.006m
        double radius = 0.0002 + 0.0001 * (static_cast<double>(i % 10) / 10.0);
        auto particle = std::make_shared<CoffeeParticle>(x, y, radius);
        bed->add_particle(particle);
    }

    auto flow = std::make_shared<WaterFlow>(58, 30);

    // Create pluggable models
    auto permeability_model = std::make_shared<KozenyCarmanPermeabilityModel>(1e-12);
    auto flow_model = std::make_shared<DarcyFlowModel>(permeability_model);
    auto transport_model = std::make_shared<ExponentialExtractionModel>(0.002, 95.0);

    // Create solver with custom transport model
    PhysicsSolver solver(bed, flow, params, flow_model, transport_model);

    // Run 30 seconds of simulation
    double dt = 0.01;
    int steps = static_cast<int>(30.0 / dt);

    for (int i = 0; i < steps; ++i) {
        solver.simulate_step(dt);
    }

    // Report results
    double total_extraction = 0.0;
    auto& particles = bed->get_particles();
    for (const auto& p : particles) {
        total_extraction += p->get_extraction_state();
    }
    double avg_extraction = total_extraction / static_cast<double>(particles.size());

    std::cout << "Custom transport model results:" << std::endl;
    std::cout << "  Average extraction: " << (avg_extraction * 100.0) << "%" << std::endl;
    std::cout << "  Particles: " << particles.size() << std::endl;
    std::cout << "  Steps: " << steps << std::endl;

    return 0;
}