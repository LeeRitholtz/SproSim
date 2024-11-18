#pragma once
#include "sprosim/interfaces/IParticle.h"
#include <array>

namespace sprosim {

class CoffeeParticle : public ICoffeeParticle {
public:
    CoffeeParticle(double x, double y, double radius);
    
    // IParticle interface
    std::pair<double, double> get_position() const override;
    double get_size() const override;
    
    // ICoffeeParticle interface
    double get_extraction_state() const override;
    double get_concentration() const override;
    void update_extraction(double delta_conc, double dt) override;
    void apply_force(double fx, double fy) override;
    
    // Additional methods
    void update_position(double dx, double dy);

private:
    std::array<double, 2> position_;
    double radius_;
    double extraction_state_;
    double concentration_;
    const double max_extractable_ = 0.2;  // 20% extractable content
};

} // namespace sprosim