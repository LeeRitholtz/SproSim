// Copyright (c) 2024  SproSim Labs (United States).
// All rights reserved.
//
// This file is part of SproSim (TODO: Add URL).
//
// SPDX-License-Identifier: GPL-3.0-or-later
//
// Author(s)     : Lee Ritholtz <lee@sprosim.com>

#pragma once
#include <utility>

namespace sprosim {

/**
 * @brief Interface for particles in the simulation
 *
 * Base interface defining the minimum functionality required for
 * any particle type in the simulation. Provides methods
 * for accessing physical properties like position and size.
 */
class IParticle {
  public:
    /** @brief Virtual destructor ensuring proper cleanup of derived classes */
    virtual ~IParticle() = default;

    /**
     * @brief Get particle position in 2D space
     * @return Position coordinates (x,y) in meters
     * TODO: This will need to eventually be made 3D
     */
    virtual std::pair<double, double> get_position() const = 0;

    /**
     * @brief Get particle diameter
     * @return Particle diameter in meters
     */
    virtual double get_size() const = 0;
};

/**
 * @brief Specialized interface for coffee particles
 *
 * Extends the base particle interface with *coffee-specific* properties
 * and behaviors like extraction state and concentration tracking.
 * Handles interactions with water flow and extraction physics.
 */
class ICoffeeParticle : public IParticle {
  public:
    /** @brief Virtual destructor ensuring proper cleanup */
    virtual ~ICoffeeParticle() = default;

    /**
     * @brief Get current extraction state of the particle
     * @return Value between 0.0 (unextracted) and 1.0 (fully extracted)
     */
    virtual double get_extraction_state() const = 0;

    /**
     * @brief Get concentration of dissolved solids in surrounding water
     * @return Local concentration of dissolved coffee compounds [kg/m³]
     *
     * Returns the concentration of dissolved coffee compounds in the water
     * immediately surrounding this particle. Used to calculate concentration
     * gradients for extraction physics.
     */
    virtual double get_concentration() const = 0;

    /**
     * @brief Update particle's extraction progress and local concentration
     * @param delta_conc Amount of coffee solids extracted in this step [kg/m³]
     * @param Time step (dt) [s]
     *
     * Updates both the local water concentration and the particle's internal
     * extraction state. The extraction state tracks cumulative extraction
     * progress from 0.0 (fresh) to 1.0 (fully extracted).
     */
    virtual void update_extraction(double delta_conc, double dt) = 0;

    /**
     * @brief Apply force to the particle
     * @param Force in x direction [N]
     * @param Force in y direction [N]
     */
    virtual void apply_force(double fx, double fy) = 0;
};

} // namespace sprosim
