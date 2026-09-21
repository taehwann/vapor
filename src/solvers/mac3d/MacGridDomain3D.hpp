#pragma once

#include "common/DomainInfo.hpp"

#include <cstddef>

// Geometry and indexing of the current uniform cubic MAC discretization.
// Simulation fields live in MacGridState3D, not in the domain.
class MacGridDomain3D final {
public:
    explicit MacGridDomain3D(int resolution, float boxSize = 4.5f);
    [[nodiscard]] Dimension dimension() const noexcept { return Dimension::D3; }
    [[nodiscard]] DomainBounds bounds() const noexcept {
        return {{0.f, 0.f, 0.f}, {boxSize_, boxSize_, boxSize_}};
    }
    [[nodiscard]] int resolution() const noexcept { return resolution_; }
    [[nodiscard]] float boxSize() const noexcept { return boxSize_; }
    [[nodiscard]] float cellSize() const noexcept { return boxSize_ / resolution_; }
    [[nodiscard]] std::size_t cellCount() const noexcept {
        return std::size_t(resolution_) * resolution_ * resolution_;
    }
    [[nodiscard]] std::size_t idC(int x, int y, int z) const noexcept {
        return x + std::size_t(resolution_) * (y + std::size_t(resolution_) * z);
    }
    [[nodiscard]] std::size_t idX(int x, int y, int z) const noexcept {
        return x + std::size_t(resolution_ + 1) * (y + std::size_t(resolution_) * z);
    }
    [[nodiscard]] std::size_t idY(int x, int y, int z) const noexcept {
        return x + std::size_t(resolution_) * (y + std::size_t(resolution_ + 1) * z);
    }
    [[nodiscard]] std::size_t idZ(int x, int y, int z) const noexcept { return idC(x, y, z); }
    void setBoxSize(float size);

private:
    int resolution_;
    float boxSize_;
};
