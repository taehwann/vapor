#pragma once
#include "fluid/advection/IAdvection.hpp"
#include "fluid/MacGridState2D.hpp"

enum class MacField2D { Density, VelocityX, VelocityY };

// Source and tracing velocity must remain immutable during the operation.
// Output must not overlap either of them. Coordinates are in world units.
struct MacAdvectionOperation2D {
    const MacGridDomain2D& domain;
    std::span<const float> velocityX, velocityY;
    std::span<const float> source;
    std::span<float> output;
    float dt;
    MacField2D field;
};

class SemiLagrangian2D final : public IAdvection<MacAdvectionOperation2D> {
public:
    [[nodiscard]] std::string_view name() const noexcept override { return "2D semi-Lagrangian"; }
    void advect(const MacAdvectionOperation2D& operation) const override;
};
