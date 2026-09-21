#include "solvers/mac2d/SemiLagrangian2D.hpp"
#include <algorithm>
#include <cmath>
#include <functional>
#include <stdexcept>

namespace {
bool overlaps(std::span<const float> a, std::span<const float> b) {
    const std::less<const float*> less;
    return !a.empty() && !b.empty() && less(a.data(), b.data() + b.size()) && less(b.data(), a.data() + a.size());
}
// Bilinear sampling in index coordinates; clamp to the nearest boundary sample.
float sample(std::span<const float> f, int w, int h, float x, float y) {
    if (!std::isfinite(x) || !std::isfinite(y)) throw std::invalid_argument("Non-finite 2D backtrace");
    x = std::clamp(x, 0.f, float(w - 1));
    y = std::clamp(y, 0.f, float(h - 1));
    const int i = int(std::floor(x)), j = int(std::floor(y));
    const int ip = std::min(i + 1, w - 1), jp = std::min(j + 1, h - 1);
    const float tx = x - i, ty = y - j;
    return std::lerp(std::lerp(f[i + w * j], f[ip + w * j], tx),
                     std::lerp(f[i + w * jp], f[ip + w * jp], tx), ty);
}
}

void SemiLagrangian2D::advect(const MacAdvectionOperation2D& op) const {
    const int n = op.domain.resolution();
    if (op.field != MacField2D::Density && op.field != MacField2D::VelocityX && op.field != MacField2D::VelocityY)
        throw std::invalid_argument("Unknown 2D MAC field");
    const int w = n + (op.field == MacField2D::VelocityX);
    const int h = n + (op.field == MacField2D::VelocityY);
    const auto faceCount = std::size_t(n) * (n + 1);
    if (!std::isfinite(op.dt) || op.dt < 0.f || op.velocityX.size() != faceCount ||
        op.velocityY.size() != faceCount || op.source.size() != std::size_t(w) * h ||
        op.output.size() != op.source.size() || overlaps(op.source, op.output) ||
        overlaps(op.velocityX, op.output) || overlaps(op.velocityY, op.output))
        throw std::invalid_argument("Invalid 2D advection timestep, extent, or overlapping output");
    for (auto field : {op.source, op.velocityX, op.velocityY})
        for (float value : field) if (!std::isfinite(value)) throw std::invalid_argument("Non-finite 2D advection input");

    const float ox = op.field == MacField2D::VelocityX ? 0.f : .5f;
    const float oy = op.field == MacField2D::VelocityY ? 0.f : .5f;
    const float dtOverH = op.dt / op.domain.cellSize();
    for (int y = 0; y < h; ++y) for (int x = 0; x < w; ++x) {
        const float px = x + ox, py = y + oy;
        const float u = sample(op.velocityX, n + 1, n, px, py - .5f);
        const float v = sample(op.velocityY, n, n + 1, px - .5f, py);
        op.output[x + w * y] = sample(op.source, w, h,
            px - dtOverH * u - ox, py - dtOverH * v - oy);
    }
}
