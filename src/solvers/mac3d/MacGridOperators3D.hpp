#pragma once

#include "solvers/mac3d/MacGridState3D.hpp"

namespace MacGridOperators3D {
// No transfer through the six domain walls. This does not impose no-slip tangents.
void enforceClosedBoundaries(MacGridState3D& state);
void computeDivergence(const MacGridState3D& state, std::span<float> output);
[[nodiscard]] float maxDivergence(const MacGridState3D& state);
void applyPressureGradient(MacGridState3D& state, std::span<const float> pressure, float dt);
}
