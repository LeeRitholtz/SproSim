#pragma once
#include <vector>
#include <memory>
#include "sprosim/interfaces/IParticle.h"

namespace sprosim {

class CoffeeBed {
public:
    CoffeeBed(double initial_mass_g, double bed_diameter_mm);
    
    void add_particle(std::shared_ptr<ICoffeeParticle> particle);
    const std::vector<std::shared_ptr<ICoffeeParticle>>& particles() const;
    
    double get_compaction() const;
    double get_diameter() const;
    double get_total_dissolved_solids() const;
    double get_extraction_yield() const;
    double get_bed_height() const;
    double get_porosity() const;
    
    void update_compaction(double pressure_pa);

private:
    std::vector<std::shared_ptr<ICoffeeParticle>> particles_;
    double initial_mass_g_;
    double bed_diameter_;
    double compaction_;
};

} // namespace sprosim