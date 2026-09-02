#pragma once

#include "math/Vec3.hpp"

// Simulation dimension, independent of field storage or rendering format.
enum class Dimension {
    D2 = 2,
    D3 = 3
};

struct DomainBounds {
    Vec3 min;
    Vec3 max;
};

// Geometry shared by application-facing solvers. No grid resolution, field
// storage, or choice of primal/dual variables is imposed here. A simplicial
// domain can own its oriented topology and primal/dual geometry separately.
class IDomain {
public:
    virtual ~IDomain() = default;
    [[nodiscard]] virtual Dimension dimension() const noexcept = 0;
    [[nodiscard]] virtual DomainBounds bounds() const noexcept = 0;
};
