#pragma once

#include "solvers/mac3d/VelocityFieldView.hpp"
#include "common/Vec3.hpp"

// Read-only, clamped trilinear sampling of cubic MAC fields. Coordinates are
// world-space. Field sizes must match the grid and positions must be finite.
class MacGridSampler3D {
public:
    explicit MacGridSampler3D(const MacGridState3D& grid) : state_(grid) {}
    MacGridSampler3D(const MacGridState3D&&) = delete;
    float sampleCell(std::span<const float> f, Vec3 p) const;
    float sampleVxField(std::span<const float> f, Vec3 p) const;
    float sampleVyField(std::span<const float> f, Vec3 p) const;
    float sampleVzField(std::span<const float> f, Vec3 p) const;
    Vec3 velocity(Vec3 p, VelocityFieldView field) const;

private:
    const MacGridState3D& state_;
    inline int idC(int x, int y, int z) const { return static_cast<int>(state_.idC(x, y, z)); }
    inline int idX(int x, int y, int z) const { return static_cast<int>(state_.idX(x, y, z)); }
    inline int idY(int x, int y, int z) const { return static_cast<int>(state_.idY(x, y, z)); }
    inline int idZ(int x, int y, int z) const { return static_cast<int>(state_.idZ(x, y, z)); }
    inline float h() const { return state_.cellSize(); }
};
