#include "../mocks/MockCoffeeParticle.h"
#include "../mocks/MockWaterFlow.h"
#include "sprosim/CoffeeBed.h"
#include "sprosim/Parameters.h"
#include "sprosim/models/transport/LinearExtractionModel.h"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <memory>

using namespace sprosim;
using namespace sprosim::testing;

TEST_CASE("LinearExtractionModel temperature factor calculation", "[LinearExtractionModel]") {
    LinearExtractionModel model;

    SECTION("Temperature factor validation through extraction") {
        Parameters params;
        params.temperature = 373.15;
        params.temperature_factor = 0.1;
        params.extraction_rate = 1.0;
        params.saturation_concentration = 10.0;

        // Create test setup
        auto water_flow = std::make_shared<MockWaterFlow>(2, 2);
        auto bed = std::make_shared<CoffeeBed>(18.0, 0.058);

        // Add a particle to test
        auto particle = std::make_shared<MockCoffeeParticle>(0.0005, 0.0005, 0.0005);
        bed->add_particle(particle);

        // Set up flow conditions
        water_flow->set_velocity(0, 0, 0.1, 0.0);

        model.update_extraction(water_flow, bed, params, 0.1);

        // Verify extraction occurred
        const auto& particles = bed->get_particles();
        REQUIRE(particles.size() == 1);
        // Temperature factor affects extraction rate, so we can verify it's working
        REQUIRE(particles[0]->get_extraction_state() > 0.0);
    }
}

TEST_CASE("LinearExtractionModel extraction physics", "[LinearExtractionModel]") {
    LinearExtractionModel model;
    Parameters params;
    params.temperature = 373.15; // 100°C
    params.temperature_factor = 0.1;
    params.extraction_rate = 0.5;
    params.saturation_concentration = 10.0;

    auto water_flow = std::make_shared<MockWaterFlow>(3, 3);
    auto bed = std::make_shared<CoffeeBed>(18.0, 0.058);

    SECTION("Basic extraction with no flow") {
        // Add particle at center
        auto particle = std::make_shared<MockCoffeeParticle>(0.0015, 0.0015, 0.0005);
        bed->add_particle(particle);

        // No flow velocity
        water_flow->set_velocity(1, 1, 0.0, 0.0);

        double dt = 0.01;
        model.update_extraction(water_flow, bed, params, dt);

        const auto& particles = bed->get_particles();
        REQUIRE(particles.size() == 1);

        // Verify extraction occurred
        REQUIRE(particles[0]->get_extraction_state() > 0.0);
        REQUIRE(particles[0]->get_concentration() > 0.0);
    }

    SECTION("Flow enhancement increases extraction") {
        // Add two identical particles
        auto particle1 = std::make_shared<MockCoffeeParticle>(0.0005, 0.0005, 0.0005);
        auto particle2 = std::make_shared<MockCoffeeParticle>(0.0015, 0.0015, 0.0005);
        bed->add_particle(particle1);
        bed->add_particle(particle2);

        // Set different flow velocities
        water_flow->set_velocity(0, 0, 0.0, 0.0); // No flow
        water_flow->set_velocity(1, 1, 0.5, 0.0); // High flow

        double dt = 0.001;
        model.update_extraction(water_flow, bed, params, dt);

        const auto& particles = bed->get_particles();
        REQUIRE(particles.size() == 2);

        // Particle in high flow should have higher extraction
        double low_flow_extraction = particles[0]->get_extraction_state();
        double high_flow_extraction = particles[1]->get_extraction_state();

        REQUIRE(high_flow_extraction > low_flow_extraction);
    }

    SECTION("Concentration gradient drives extraction") {
        auto particle = std::make_shared<MockCoffeeParticle>(0.0005, 0.0005, 0.0005);
        bed->add_particle(particle);
        water_flow->set_velocity(0, 0, 0.1, 0.0);

        const auto& particles = bed->get_particles();
        double initial_concentration = particles[0]->get_concentration();
        double initial_extraction_state = particles[0]->get_extraction_state();

        // Run extraction
        double dt = 0.01;
        model.update_extraction(water_flow, bed, params, dt);

        // Verify changes
        REQUIRE(particles[0]->get_concentration() > initial_concentration);
        REQUIRE(particles[0]->get_extraction_state() > initial_extraction_state);
    }
}

TEST_CASE("LinearExtractionModel boundary conditions", "[LinearExtractionModel]") {
    LinearExtractionModel model;
    Parameters params;
    params.temperature = 373.15;
    params.temperature_factor = 0.1;
    params.extraction_rate = 1.0;
    params.saturation_concentration = 10.0;

    auto water_flow = std::make_shared<MockWaterFlow>(2, 2);
    auto bed = std::make_shared<CoffeeBed>(18.0, 0.058);

    SECTION("Particles outside grid are ignored") {
        // Add particle outside grid boundaries
        auto particle = std::make_shared<MockCoffeeParticle>(0.003, 0.003, 0.0005);
        bed->add_particle(particle);

        water_flow->set_velocity(0, 0, 0.1, 0.0);

        double dt = 0.01;
        model.update_extraction(water_flow, bed, params, dt);

        const auto& particles = bed->get_particles();
        REQUIRE(particles.size() == 1);

        // Particle should remain unchanged
        REQUIRE(particles[0]->get_extraction_state() == 0.0);
        REQUIRE(particles[0]->get_concentration() == 0.0);
    }

    SECTION("Zero time step produces no change") {
        auto particle = std::make_shared<MockCoffeeParticle>(0.0005, 0.0005, 0.0005);
        bed->add_particle(particle);
        water_flow->set_velocity(0, 0, 0.1, 0.0);

        double dt = 0.0;
        model.update_extraction(water_flow, bed, params, dt);

        const auto& particles = bed->get_particles();
        REQUIRE(particles[0]->get_extraction_state() == 0.0);
        REQUIRE(particles[0]->get_concentration() == 0.0);
    }
}

TEST_CASE("LinearExtractionModel parameter sensitivity", "[LinearExtractionModel]") {
    LinearExtractionModel model;

    SECTION("Higher extraction rate increases extraction") {
        auto water_flow = std::make_shared<MockWaterFlow>(2, 2);
        auto bed1 = std::make_shared<CoffeeBed>(18.0, 0.058);
        auto bed2 = std::make_shared<CoffeeBed>(18.0, 0.058);

        auto particle1 = std::make_shared<MockCoffeeParticle>(0.0005, 0.0005, 0.0005);
        auto particle2 = std::make_shared<MockCoffeeParticle>(0.0005, 0.0005, 0.0005);
        bed1->add_particle(particle1);
        bed2->add_particle(particle2);

        water_flow->set_velocity(0, 0, 0.1, 0.0);

        // Low extraction rate
        Parameters params1;
        params1.temperature = 373.15;
        params1.temperature_factor = 0.1;
        params1.extraction_rate = 1.0;
        params1.saturation_concentration = 10.0;

        // High extraction rate
        Parameters params2 = params1;
        params2.extraction_rate = 5.0;

        double dt = 0.01;
        model.update_extraction(water_flow, bed1, params1, dt);
        model.update_extraction(water_flow, bed2, params2, dt);

        double low_rate_extraction = bed1->get_particles()[0]->get_extraction_state();
        double high_rate_extraction = bed2->get_particles()[0]->get_extraction_state();

        REQUIRE(high_rate_extraction > low_rate_extraction);
    }

    SECTION("Higher temperature increases extraction") {
        auto water_flow = std::make_shared<MockWaterFlow>(2, 2);
        auto bed1 = std::make_shared<CoffeeBed>(18.0, 0.058);
        auto bed2 = std::make_shared<CoffeeBed>(18.0, 0.058);

        auto particle1 = std::make_shared<MockCoffeeParticle>(0.0005, 0.0005, 0.0005);
        auto particle2 = std::make_shared<MockCoffeeParticle>(0.0005, 0.0005, 0.0005);
        bed1->add_particle(particle1);
        bed2->add_particle(particle2);

        water_flow->set_velocity(0, 0, 0.1, 0.0);

        // Lower temperature
        Parameters params1;
        params1.temperature = 363.15; // 90°C
        params1.temperature_factor = 0.1;
        params1.extraction_rate = 1.0;
        params1.saturation_concentration = 10.0;

        // Higher temperature
        Parameters params2 = params1;
        params2.temperature = 383.15; // 110°C

        double dt = 0.01;
        model.update_extraction(water_flow, bed1, params1, dt);
        model.update_extraction(water_flow, bed2, params2, dt);

        double low_temp_extraction = bed1->get_particles()[0]->get_extraction_state();
        double high_temp_extraction = bed2->get_particles()[0]->get_extraction_state();

        REQUIRE(high_temp_extraction > low_temp_extraction);
    }
}
