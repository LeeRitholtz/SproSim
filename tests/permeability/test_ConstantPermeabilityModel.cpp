#include "sprosim/models/permeability/ConstantPermeabilityModel.h"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using namespace sprosim;

TEST_CASE("ConstantPermeabilityModel returns fixed permeability", "[ConstantPermeabilityModel]") {
    const double fixed_permeability = 2e-12;
    ConstantPermeabilityModel model(fixed_permeability);

    SECTION("Calculate permeability ignores porosity") {
        double porosities[] = {0.0, 0.1, 0.5, 0.9, 1.0, 100.0};
        for (auto porosity : porosities) {
            double k = model.calculate_permeability(porosity);
            REQUIRE_THAT(k, Catch::Matchers::WithinAbs(fixed_permeability, 1e-15));
        }
    }
}
