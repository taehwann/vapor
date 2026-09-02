#pragma once

#include "fluid/MacGridState.hpp"

namespace MacGridOperators {
// No transfer through the six domain walls. This does not impose no-slip tangents.
void enforceClosedBoundaries(MacGridState& state);
void computeDivergence(const MacGridState& state, std::span<float> output);
[[nodiscard]] float maxDivergence(const MacGridState& state);
void applyPressureGradient(MacGridState& state, std::span<const float> pressure, float dt);
}
