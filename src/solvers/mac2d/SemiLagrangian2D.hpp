#pragma once
#include "solvers/mac2d/MacGridState2D.hpp"

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

class SemiLagrangian2D final {
public:
    void advect(const MacAdvectionOperation2D& operation) const;
};
