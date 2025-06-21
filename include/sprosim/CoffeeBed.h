#pragma once
#include "sprosim/interfaces/IParticle.h"
#include <memory>
#include <vector>

namespace sprosim {

/**
 * @brief Container for coffee particles that represents the physical coffee bed
 *
 * The CoffeeBed manages a collection of coffee particles and tracks their
 * collective properties during brewing simulation. It handles bed compaction
 * under pressure, calculates extraction metrics, and provides physical
 * properties like porosity and bed height.
 *
 * @example
 * ```cpp
 * CoffeeBed bed(18.0, 0.058);  // 18g coffee, 58mm basket
 * bed.add_particle(std::make_shared<CoffeeParticle3D>(...));
 * double yield = bed.get_extraction_yield();
 * ```
 */
class CoffeeBed {
public:
  CoffeeBed(double initial_mass, double bed_diameter);

  void add_particle(std::shared_ptr<ICoffeeParticle> particle);
  const std::vector<std::shared_ptr<ICoffeeParticle>> &get_particles() const;

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
