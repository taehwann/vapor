#pragma once
#include "solvers/mac2d/MacGridState2D.hpp"

namespace MacGridOperators2D {
void enforceClosedBoundaries(MacGridState2D& state);
void computeDivergence(const MacGridState2D& state, std::span<float> output);
[[nodiscard]] float maxDivergence(const MacGridState2D& state);
void applyPressureGradient(MacGridState2D& state, std::span<const double> pressure, float dt);
}
