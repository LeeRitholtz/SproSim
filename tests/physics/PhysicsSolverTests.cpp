#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "sprosim/PhysicsSolver.h"
#include "sprosim/CoffeeBed.h"
#include "../mocks/MockCoffeeParticle.h"
#include "../mocks/MockWaterFlow.h"
//#include <iostream>
//#include <cmath>

using namespace sprosim;
using namespace sprosim::testing;

TEST_CASE("PhysicsSolver basic flow field test", "[physics]") {
    auto bed = std::make_shared<CoffeeBed>(18.0, 58.0);
    auto particle = std::make_shared<MockCoffeeParticle>(0.005, 0.02, 0.0005);
    bed->add_particle(particle);

    // Verify particle position
    auto [px, py] = particle->get_position();
    printf("\nParticle check:\n");
    printf("  Direct position: %.6f, %.6f\n", px, py);
    printf("  Bed height before solver: %.6f\n", bed->get_bed_height());

    auto flow = std::make_shared<MockWaterFlow>(10, 20);

    PhysicsSolver::Parameters params{
        .permeability = 1e-12,
        .fluid_viscosity = 1e-3,
        .extraction_rate = 0.1,
        .temperature = 368.15,
        .inlet_pressure = 9e5,
        .outlet_pressure = 1e5,
        .particle_drag = 0.1,
        .flow_resistance = 0.2,
        .saturation_concentration = 0.2,
        .temperature_factor = 0.01
    };

    PhysicsSolver solver(bed, flow, params);

    // Check bed height right before test
    printf("  Bed height right before test: %.6f\n", bed->get_bed_height());

    SECTION("Flow field follows Darcy's law") {
        solver.simulate_step(0.001);
        // ... rest of test
    }
}

TEST_CASE("PhysicsSolver extraction test", "[physics]") {
    auto bed = std::make_shared<CoffeeBed>(18.0, 58.0);
    auto particle = std::make_shared<MockCoffeeParticle>(0.005, 0.01, 0.0005);
    bed->add_particle(particle);

    auto flow = std::make_shared<MockWaterFlow>(10, 20);

    PhysicsSolver::Parameters params{
        .permeability = 1e-12,
        .fluid_viscosity = 1e-3,
        .extraction_rate = 0.1,
        .temperature = 368.15,
        .inlet_pressure = 9e5,
        .outlet_pressure = 1e5,
        .particle_drag = 0.1,
        .flow_resistance = 0.2,
        .saturation_concentration = 0.2,
        .temperature_factor = 0.01
    };

    PhysicsSolver solver(bed, flow, params);

    SECTION("Extraction increases with time") {
        double initial_state = particle->get_extraction_state();

        for (int i = 0; i < 1000; i++) {
            solver.simulate_step(0.001);
        }

        REQUIRE(particle->get_extraction_state() > initial_state);
        REQUIRE(particle->get_extraction_state() < 1.0);
    }
}
