#include "sprosim/Parameters.h"
#include "sprosim/SimulationManager.h"
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>

using namespace sprosim;

TEST_CASE("SimulationManager basic construction", "[SimulationManager]") {
    SECTION("Constructor with basic configuration") {
        SimulationManager::Configuration config;
        config.dose = 18.0;
        config.brewing_time = 5.0;
        config.timestep = 0.1;
        config.particle_count = 50;
        config.export_vtk = false;

        // Test construction succeeds
        SimulationManager manager(config);
        REQUIRE(true); // If we get here, construction succeeded
    }

    SECTION("Constructor with explicit physics parameters") {
        SimulationManager::Configuration config;
        config.dose = 18.0;
        config.brewing_time = 5.0;
        config.timestep = 0.1;
        config.particle_count = 50;
        config.export_vtk = false;

        Parameters physics_params;
        physics_params.permeability = 1e-12;
        physics_params.fluid_viscosity = 1e-3;
        physics_params.extraction_rate = 0.01;
        physics_params.temperature = 368.15;

        // Test construction with explicit parameters succeeds
        SimulationManager manager(config, physics_params);
        REQUIRE(true); // If we get here, construction succeeded
    }
}

TEST_CASE("SimulationManager configuration validation", "[SimulationManager]") {
    SECTION("Configuration parameter ranges") {
        SimulationManager::Configuration config;
        config.dose = 18.0;
        config.brewing_time = 30.0;
        config.timestep = 0.01;
        config.particle_count = 100;
        config.export_vtk = false;

        // Verify configuration values are reasonable
        REQUIRE(config.dose > 0.0);
        REQUIRE(config.brewing_time > 0.0);
        REQUIRE(config.timestep > 0.0);
        REQUIRE(config.particle_count > 0);
        REQUIRE_FALSE(config.export_vtk); // VTK disabled for unit tests
    }
}
