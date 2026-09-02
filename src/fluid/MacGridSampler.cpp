#include "fluid/MacGridSampler.hpp"

#include <algorithm>

float MacGridSampler::sampleCell(std::span<const float> f, Vec3 p) const {
    float gx = std::clamp(p.x / h() - 0.5f, 0.f, float(state_.resolution() - 1));
    float gy = std::clamp(p.y / h() - 0.5f, 0.f, float(state_.resolution() - 1));
    float gz = std::clamp(p.z / h() - 0.5f, 0.f, float(state_.resolution() - 1));
    const int x0 = int(gx), y0 = int(gy), z0 = int(gz);
    const int x1 = std::min(x0 + 1, state_.resolution() - 1), y1 = std::min(y0 + 1, state_.resolution() - 1), z1 = std::min(z0 + 1, state_.resolution() - 1);
    const float tx = gx - x0, ty = gy - y0, tz = gz - z0;
    return f[idC(x0, y0, z0)] * (1 - tx) * (1 - ty) * (1 - tz) +
        f[idC(x1, y0, z0)] * tx * (1 - ty) * (1 - tz) +
        f[idC(x0, y1, z0)] * (1 - tx) * ty * (1 - tz) +
        f[idC(x1, y1, z0)] * tx * ty * (1 - tz) +
        f[idC(x0, y0, z1)] * (1 - tx) * (1 - ty) * tz +
        f[idC(x1, y0, z1)] * tx * (1 - ty) * tz +
        f[idC(x0, y1, z1)] * (1 - tx) * ty * tz +
        f[idC(x1, y1, z1)] * tx * ty * tz;
}

float MacGridSampler::sampleVxField(std::span<const float> f, Vec3 p) const {
    const float gx = std::clamp(p.x / h(), 0.f, float(state_.resolution()));
    const float gy = std::clamp(p.y / h() - 0.5f, 0.f, float(state_.resolution() - 1));
    const float gz = std::clamp(p.z / h() - 0.5f, 0.f, float(state_.resolution() - 1));
    const int x0 = std::min(int(gx), state_.resolution() - 1), y0 = int(gy), z0 = int(gz);
    const int x1 = x0 + 1, y1 = std::min(y0 + 1, state_.resolution() - 1), z1 = std::min(z0 + 1, state_.resolution() - 1);
    const float tx = gx - x0, ty = gy - y0, tz = gz - z0;
    return f[idX(x0, y0, z0)] * (1 - tx) * (1 - ty) * (1 - tz) +
        f[idX(x1, y0, z0)] * tx * (1 - ty) * (1 - tz) +
        f[idX(x0, y1, z0)] * (1 - tx) * ty * (1 - tz) +
        f[idX(x1, y1, z0)] * tx * ty * (1 - tz) +
        f[idX(x0, y0, z1)] * (1 - tx) * (1 - ty) * tz +
        f[idX(x1, y0, z1)] * tx * (1 - ty) * tz +
        f[idX(x0, y1, z1)] * (1 - tx) * ty * tz +
        f[idX(x1, y1, z1)] * tx * ty * tz;
}

float MacGridSampler::sampleVyField(std::span<const float> f, Vec3 p) const {
    const float gx = std::clamp(p.x / h() - 0.5f, 0.f, float(state_.resolution() - 1));
    const float gy = std::clamp(p.y / h(), 0.f, float(state_.resolution()));
    const float gz = std::clamp(p.z / h() - 0.5f, 0.f, float(state_.resolution() - 1));
    const int x0 = int(gx), y0 = std::min(int(gy), state_.resolution() - 1), z0 = int(gz);
    const int x1 = std::min(x0 + 1, state_.resolution() - 1), y1 = y0 + 1, z1 = std::min(z0 + 1, state_.resolution() - 1);
    const float tx = gx - x0, ty = gy - y0, tz = gz - z0;
    return f[idY(x0, y0, z0)] * (1 - tx) * (1 - ty) * (1 - tz) +
        f[idY(x1, y0, z0)] * tx * (1 - ty) * (1 - tz) +
        f[idY(x0, y1, z0)] * (1 - tx) * ty * (1 - tz) +
        f[idY(x1, y1, z0)] * tx * ty * (1 - tz) +
        f[idY(x0, y0, z1)] * (1 - tx) * (1 - ty) * tz +
        f[idY(x1, y0, z1)] * tx * (1 - ty) * tz +
        f[idY(x0, y1, z1)] * (1 - tx) * ty * tz +
        f[idY(x1, y1, z1)] * tx * ty * tz;
}

float MacGridSampler::sampleVzField(std::span<const float> f, Vec3 p) const {
    const float gx = std::clamp(p.x / h() - 0.5f, 0.f, float(state_.resolution() - 1));
    const float gy = std::clamp(p.y / h() - 0.5f, 0.f, float(state_.resolution() - 1));
    const float gz = std::clamp(p.z / h(), 0.f, float(state_.resolution()));
    const int x0 = int(gx), y0 = int(gy), z0 = std::min(int(gz), state_.resolution() - 1);
    const int x1 = std::min(x0 + 1, state_.resolution() - 1), y1 = std::min(y0 + 1, state_.resolution() - 1), z1 = z0 + 1;
    const float tx = gx - x0, ty = gy - y0, tz = gz - z0;
    return f[idZ(x0, y0, z0)] * (1 - tx) * (1 - ty) * (1 - tz) +
        f[idZ(x1, y0, z0)] * tx * (1 - ty) * (1 - tz) +
        f[idZ(x0, y1, z0)] * (1 - tx) * ty * (1 - tz) +
        f[idZ(x1, y1, z0)] * tx * ty * (1 - tz) +
        f[idZ(x0, y0, z1)] * (1 - tx) * (1 - ty) * tz +
        f[idZ(x1, y0, z1)] * tx * (1 - ty) * tz +
        f[idZ(x0, y1, z1)] * (1 - tx) * ty * tz +
        f[idZ(x1, y1, z1)] * tx * ty * tz;
}

Vec3 MacGridSampler::velocity(Vec3 p, VelocityFieldView field) const {
    return {sampleVxField(field.x, p), sampleVyField(field.y, p), sampleVzField(field.z, p)};
}
