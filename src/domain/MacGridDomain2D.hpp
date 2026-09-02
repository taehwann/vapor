#pragma once

#include "domain/IDomain.hpp"
#include <cstddef>

// Uniform square grid in the XY plane. Y points upward; no Z field is stored.
class MacGridDomain2D final : public IDomain {
public:
    explicit MacGridDomain2D(int resolution, float boxSize = 4.f);
    [[nodiscard]] Dimension dimension() const noexcept override { return Dimension::D2; }
    [[nodiscard]] DomainBounds bounds() const noexcept override {
        return {{0.f, 0.f, 0.f}, {boxSize_, boxSize_, 0.f}};
    }
    [[nodiscard]] int resolution() const noexcept { return resolution_; }
    [[nodiscard]] float boxSize() const noexcept { return boxSize_; }
    [[nodiscard]] float cellSize() const noexcept { return boxSize_ / resolution_; }
    [[nodiscard]] std::size_t idC(int x, int y) const noexcept { return x + std::size_t(resolution_) * y; }
    [[nodiscard]] std::size_t idX(int x, int y) const noexcept { return x + std::size_t(resolution_ + 1) * y; }
    [[nodiscard]] std::size_t idY(int x, int y) const noexcept { return idC(x, y); }
    void setBoxSize(float size);
private:
    int resolution_;
    float boxSize_;
};
