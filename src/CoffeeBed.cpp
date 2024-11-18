#include "sprosim/CoffeeBed.h"
#include "sprosim/interfaces/IParticle.h"
#include <algorithm>

namespace sprosim {

CoffeeBed::CoffeeBed(double initial_mass_g, double bed_diameter_mm):
    initial_mass_g_(initial_mass_g),
    bed_diameter_(bed_diameter_mm / 1000.0),
    compaction_(1.0)
{}

void CoffeeBed::add_particle(std::shared_ptr<ICoffeeParticle> particle) {
    particles_.push_back(particle);
}

const std::vector<std::shared_ptr<ICoffeeParticle>>& CoffeeBed::particles() const {
    return particles_;
}

double CoffeeBed::get_compaction() const {
    return compaction_;
}

double CoffeeBed::get_diameter() const {
    return bed_diameter_;
}

double CoffeeBed::get_total_dissolved_solids() const {
    double tds = 0.0;
    for (const auto& particle : particles_) {
        tds += particle->get_extraction_state() * (initial_mass_g_ / particles_.size());
    }
    return tds;
}

double CoffeeBed::get_extraction_yield() const {
    return get_total_dissolved_solids() / initial_mass_g_ * 100.0;
}

void CoffeeBed::update_compaction(double pressure_pa) {
    const double k_compact = 1e-7;
    compaction_ = 1.0 + k_compact * pressure_pa;
}

double CoffeeBed::get_bed_height() const {
    if (particles_.empty()) return 0.0;

    auto [_, max_y] = particles_[0]->get_position();
    for (const auto& particle : particles_) {
        auto [__, y] = particle->get_position();
        max_y = std::max(max_y, y);
    }
    return max_y;
}

double CoffeeBed::get_porosity() const {
    return 0.4 / compaction_;
}

} // namespace sprosim
