#pragma once

#include "fluid/VelocityFieldView.hpp"
#include "math/Vec3.hpp"

// Read-only, clamped trilinear sampling of cubic MAC fields. Coordinates are
// world-space. Field sizes must match the grid and positions must be finite.
class MacGridSampler {
public:
    explicit MacGridSampler(const MacGridState& grid) : state_(grid) {}
    MacGridSampler(const MacGridState&&) = delete;
    float sampleCell(std::span<const float> f, Vec3 p) const;
    float sampleVxField(std::span<const float> f, Vec3 p) const;
    float sampleVyField(std::span<const float> f, Vec3 p) const;
    float sampleVzField(std::span<const float> f, Vec3 p) const;
    Vec3 velocity(Vec3 p, VelocityFieldView field) const;

private:
    const MacGridState& state_;
    inline int idC(int x, int y, int z) const { return static_cast<int>(state_.idC(x, y, z)); }
    inline int idX(int x, int y, int z) const { return static_cast<int>(state_.idX(x, y, z)); }
    inline int idY(int x, int y, int z) const { return static_cast<int>(state_.idY(x, y, z)); }
    inline int idZ(int x, int y, int z) const { return static_cast<int>(state_.idZ(x, y, z)); }
    inline float h() const { return state_.cellSize(); }
};
