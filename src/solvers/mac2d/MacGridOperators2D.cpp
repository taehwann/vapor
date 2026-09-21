#include "solvers/mac2d/MacGridOperators2D.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace {
float divergenceAt(const MacGridState2D& s, int x, int y) {
    const float value = (s.velocityX()[s.idX(x + 1, y)] - s.velocityX()[s.idX(x, y)] +
                         s.velocityY()[s.idY(x, y + 1)] - s.velocityY()[s.idY(x, y)]) / s.cellSize();
    if (!std::isfinite(value)) throw std::invalid_argument("Non-finite 2D MAC divergence");
    return value;
}
}
void MacGridOperators2D::enforceClosedBoundaries(MacGridState2D& s) {
    const int n = s.resolution();
    for (int i = 0; i < n; ++i) {
        s.velocityX()[s.idX(0, i)] = s.velocityX()[s.idX(n, i)] = 0.f;
        s.velocityY()[s.idY(i, 0)] = s.velocityY()[s.idY(i, n)] = 0.f;
    }
}
void MacGridOperators2D::computeDivergence(const MacGridState2D& s, std::span<float> output) {
    if (output.size() != s.cellCount()) throw std::invalid_argument("2D divergence output size mismatch");
    for (int y = 0; y < s.resolution(); ++y) for (int x = 0; x < s.resolution(); ++x)
        output[s.idC(x, y)] = divergenceAt(s, x, y);
}
float MacGridOperators2D::maxDivergence(const MacGridState2D& s) {
    float result = 0.f;
    for (int y = 0; y < s.resolution(); ++y) for (int x = 0; x < s.resolution(); ++x)
        result = std::max(result, std::abs(divergenceAt(s, x, y)));
    return result;
}
void MacGridOperators2D::applyPressureGradient(MacGridState2D& s, std::span<const double> p, float dt) {
    if (p.size() != s.cellCount() || !std::isfinite(dt) || dt <= 0.f)
        throw std::invalid_argument("Invalid 2D pressure field or time step");
    const int n = s.resolution();
    const double scale = double(dt) / s.cellSize();
    for (int y = 0; y < n; ++y) for (int x = 1; x < n; ++x)
        s.velocityX()[s.idX(x, y)] = float(s.velocityX()[s.idX(x, y)] - scale * (p[s.idC(x, y)] - p[s.idC(x - 1, y)]));
    for (int y = 1; y < n; ++y) for (int x = 0; x < n; ++x)
        s.velocityY()[s.idY(x, y)] = float(s.velocityY()[s.idY(x, y)] - scale * (p[s.idC(x, y)] - p[s.idC(x, y - 1)]));
}
