#include "fluid/MacGridOperators.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace {
float divergenceAt(const MacGridState& s, int x, int y, int z) {
    const auto u = s.velocityX(), v = s.velocityY(), w = s.velocityZ();
    return (u[s.idX(x + 1, y, z)] - u[s.idX(x, y, z)] +
            v[s.idY(x, y + 1, z)] - v[s.idY(x, y, z)] +
            w[s.idZ(x, y, z + 1)] - w[s.idZ(x, y, z)]) / s.cellSize();
}
}

void MacGridOperators::enforceClosedBoundaries(MacGridState& s) {
    const int n = s.resolution();
    auto u = s.velocityX(), v = s.velocityY(), w = s.velocityZ();
    for (int z = 0; z < n; ++z) for (int y = 0; y < n; ++y) {
        u[s.idX(0, y, z)] = u[s.idX(n, y, z)] = 0.f;
    }
    for (int z = 0; z < n; ++z) for (int x = 0; x < n; ++x) {
        v[s.idY(x, 0, z)] = v[s.idY(x, n, z)] = 0.f;
    }
    for (int y = 0; y < n; ++y) for (int x = 0; x < n; ++x) {
        w[s.idZ(x, y, 0)] = w[s.idZ(x, y, n)] = 0.f;
    }
}

void MacGridOperators::computeDivergence(const MacGridState& s, std::span<float> output) {
    if (output.size() != s.cellCount()) {
        throw std::invalid_argument("Divergence output must match the cell count");
    }
    const int n = s.resolution();
    for (int z = 0; z < n; ++z) for (int y = 0; y < n; ++y) for (int x = 0; x < n; ++x) {
        output[s.idC(x, y, z)] = divergenceAt(s, x, y, z);
    }
}

float MacGridOperators::maxDivergence(const MacGridState& s) {
    float result = 0.f;
    const int n = s.resolution();
    for (int z = 0; z < n; ++z) for (int y = 0; y < n; ++y) for (int x = 0; x < n; ++x) {
        const float value = divergenceAt(s, x, y, z);
        if (!std::isfinite(value)) throw std::invalid_argument("Non-finite MAC velocity");
        result = std::max(result, std::abs(value));
    }
    return result;
}

void MacGridOperators::applyPressureGradient(MacGridState& s, std::span<const float> p, float dt) {
    if (p.size() != s.cellCount() || !std::isfinite(dt) || dt <= 0.f) {
        throw std::invalid_argument("Invalid pressure field or time step");
    }
    const int n = s.resolution();
    const float hInv = 1.f / s.cellSize();
    auto u = s.velocityX(), v = s.velocityY(), w = s.velocityZ();
    for (int z = 0; z < n; ++z) for (int y = 0; y < n; ++y) for (int x = 1; x < n; ++x) {
        u[s.idX(x, y, z)] -= dt * (p[s.idC(x, y, z)] - p[s.idC(x - 1, y, z)]) * hInv;
    }
    for (int z = 0; z < n; ++z) for (int y = 1; y < n; ++y) for (int x = 0; x < n; ++x) {
        v[s.idY(x, y, z)] -= dt * (p[s.idC(x, y, z)] - p[s.idC(x, y - 1, z)]) * hInv;
    }
    for (int z = 1; z < n; ++z) for (int y = 0; y < n; ++y) for (int x = 0; x < n; ++x) {
        w[s.idZ(x, y, z)] -= dt * (p[s.idC(x, y, z)] - p[s.idC(x, y, z - 1)]) * hInv;
    }
}
