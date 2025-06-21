// Copyright (c) 2024  SproSim Labs (United States).
// All rights reserved.
//
// This file is part of SproSim (TODO: Add URL).
//
// SPDX-License-Identifier: GPL-3.0-or-later
//
// Author(s)     : Lee Ritholtz <lee@sprosim.com>

#pragma once
#include <cstddef>
#include <utility>

namespace sprosim {

/**
 * @brief Base interface for flow field implementations
 *
 * Defines the common functionality for any flow field in the simulation.
 * Uses a fixed Cartesian grid for field quantities. The grid represents
 * a physical domain where each cell corresponds to a finite volume in
 * the espresso bed.
 */
class IFlow {
  public:
    /** @brief Virtual destructor ensuring proper cleanup of derived classes */
    virtual ~IFlow() = default;

    /**
     * @brief Get the dimensions of the flow field grid
     * @return Pair containing (width, height) in grid cells
     *
     * Grid dimensions define the computational domain resolution.
     * Width corresponds to basket diameter, height to bed depth.
     */
    virtual std::pair<size_t, size_t> get_grid_dimensions() const = 0;

    /**
     * @brief Get the physical size of each grid cell
     * @return Cell size in meters
     *
     * Cell size determines spatial discretization. Smaller cells increase
     * accuracy but also computational cost. Typically on the order of
     * particle diameter (~1mm) for accurate flow resolution.
     */
    virtual double get_cell_size() const = 0;
};

/**
 * @brief Specialized interface for water flow fields
 *
 * Extends the base flow interface with methods specific to water flow
 * simulation. Handles velocity field updates and concentration tracking
 * for dissolved coffee compounds. Uses a staggered grid arrangement
 * where velocities are defined at cell faces and scalars at cell centers.
 */
class IWaterFlow : public IFlow {
  public:
    /** @brief Virtual destructor ensuring proper cleanup */
    virtual ~IWaterFlow() = default;

    /**
     * @brief Set velocity components at a grid point
     * @param Grid x-coordinate
     * @param Grid y-coordinate
     * @param X-component of velocity [m/s]
     * @param Y-component of velocity [m/s]
     *
     * Updates velocity components at the specified grid location.
     * Velocities represent volume-averaged flow through cell faces,
     * consistent with the finite volume discretization.
     */
    virtual void set_velocity(size_t i, size_t j, double vx, double vy) = 0;

    /**
     * @brief Get velocity components at a grid point
     * @param Grid x-coordinate
     * @param Grid y-coordinate
     * @return Pair containing (vx, vy) velocity components [m/s]
     *
     * Retrieves local flow velocities. These determine both extraction
     * rates and particle drag forces in the simulation.
     */
    virtual std::pair<double, double> get_velocity(size_t i, size_t j) const = 0;

    /**
     * @brief Add to the concentration of dissolved solids at a grid point
     * @param Grid x-coordinate
     * @param Grid y-coordinate
     * @param Change in concentration [kg/m³]
     *
     * Accumulates dissolved coffee compounds in the flow field.
     * Concentration is cell-centered and represents the local
     * density of extracted material. Used to model saturation
     * effects and calculate extraction gradients.
     */
    virtual void add_concentration(size_t i, size_t j, double delta) = 0;
};

} // namespace sprosim
