#pragma once

namespace sprosim {

/**
 * @brief Abstract interface for coffee bed implementations
 *
 * Defines the essential properties and behaviors that flow models
 * need to access from a coffee bed during brewing simulation.
 * This interface enables dependency inversion and better testability
 * for flow model implementations.
 */
class ICoffeeBed {
  public:
    virtual ~ICoffeeBed() = default;

    /**
     * @brief Get the current height of the coffee bed
     * @return Bed height in meters
     */
    virtual double get_bed_height() const = 0;

    /**
     * @brief Get the porosity of the coffee bed
     * @return Porosity as a fraction (0.0 to 1.0)
     */
    virtual double get_porosity() const = 0;
};

} // namespace sprosim
