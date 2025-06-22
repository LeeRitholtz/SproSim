#include "sprosim/CoffeeBed.h"
#include "sprosim/CoffeeParticle3D.h"
#include "sprosim/interfaces/IParticle.h"
#include <algorithm>
#include <memory>

namespace sprosim {

CoffeeBed::CoffeeBed(double initial_mass, double bed_diameter)
    : initial_mass_g_(initial_mass * 1000.0), // Convert kg to g for internal storage
      bed_diameter_(bed_diameter),            // Already in meters
      compaction_(1.0) {}                     // 1.0 = Not compacted at all

void CoffeeBed::add_particle(std::shared_ptr<ICoffeeParticle> particle) {
    particles_.push_back(particle);
}

const std::vector<std::shared_ptr<ICoffeeParticle>>& CoffeeBed::get_particles() const {
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
        /* Assumption: Each particle contributes equally to the total mass
         * TODO: Update particle creation to have different sizes. This would
         *      be more realistic since there are always fines
         */
        tds += particle->get_extraction_state() * (initial_mass_g_ / particles_.size());
    }
    return tds;
}

double CoffeeBed::get_extraction_yield() const {
    return get_total_dissolved_solids() / initial_mass_g_ * 100.0;
}

void CoffeeBed::update_compaction(double pressure_pa) {
    /* Provides a scaling factor for pressure_pa, which is in Pascals
     *  TODO: This should just be in bars to begin with. No one would input
     *      a number using Pascals. Pressure in bars is more familiar to
     *      coffee people.
     */
    const double k_compact = 1e-7;
    compaction_ = 1.0 + k_compact * pressure_pa;
}

double CoffeeBed::get_bed_height() const {
    if (particles_.empty())
        return 0.0;

    double max_height = 0.0;

    for (const auto& particle : particles_) {
        // Try to cast to 3D particle first
        auto particle_3d = std::dynamic_pointer_cast<CoffeeParticle3D>(particle);
        if (particle_3d) {
            // For 3D particles, use Z coordinate (depth) as height
            auto [x, y, z] = particle_3d->get_position_3d();
            double particle_top = z + particle_3d->get_size() / 2.0; // Include particle radius
            max_height = std::max(max_height, particle_top);
        } else {
            // Fallback to 2D particle Y coordinate
            auto [x, y] = particle->get_position();
            max_height = std::max(max_height, y);
        }
    }

    return max_height;
}

double CoffeeBed::get_porosity() const {
    return 0.4 / compaction_;
}

} // namespace sprosim
