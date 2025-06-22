#include "sprosim/models/permeability/KozenyCarmanPermeabilityModel.h"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>

using namespace sprosim;

TEST_CASE("KozenyCarmanPermeabilityModel calculates effective permeability",
          "[KozenyCarmanPermeabilityModel]") {
    const double base_permeability = 1e-12;
    KozenyCarmanPermeabilityModel model(base_permeability);

    SECTION("Valid porosity values") {
        double porosities[] = {0.1, 0.3, 0.5, 0.7, 0.9};
        for (auto porosity : porosities) {
            double k_eff = model.calculate_permeability(porosity);
            REQUIRE(k_eff > 0.0);

            // Verify Kozeny-Carman formula: k = k₀ * φ³ / (1-φ)²
            double expected =
                base_permeability * std::pow(porosity, 3) / std::pow(1.0 - porosity, 2);
            REQUIRE_THAT(k_eff, Catch::Matchers::WithinAbs(expected, 1e-15));
        }
    }

    SECTION("Edge cases") {
        // Porosity near zero should yield permeability near zero
        double near_zero = 1e-6;
        double k_eff = model.calculate_permeability(near_zero);
        REQUIRE_THAT(k_eff, Catch::Matchers::WithinAbs(0.0, 1e-15));

        // Porosity near one should be larger than base permeability
        double near_one = 0.999999;
        k_eff = model.calculate_permeability(near_one);
        REQUIRE(k_eff > base_permeability);

        // Porosity zero returns zero
        k_eff = model.calculate_permeability(0.0);
        REQUIRE_THAT(k_eff, Catch::Matchers::WithinAbs(0.0, 1e-15));
    }
}
