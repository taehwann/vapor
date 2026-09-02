#pragma once

#include "fluid/MacGridState.hpp"

#include <span>

// Borrowed staggered velocity fields, usable by sampling and advection alike.
struct VelocityFieldView {
    std::span<const float> x, y, z;
};

inline VelocityFieldView velocityView(const MacGridState& state) {
    return {state.velocityX(), state.velocityY(), state.velocityZ()};
}
