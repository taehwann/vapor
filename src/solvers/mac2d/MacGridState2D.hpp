#pragma once

#include "solvers/mac2d/MacGridDomain2D.hpp"
#include <span>
#include <vector>

// Physical fields: N*N density, (N+1)*N U faces, N*(N+1) V faces.
class MacGridState2D {
public:
    explicit MacGridState2D(int resolution, float boxSize = 4.f);
    explicit MacGridState2D(MacGridDomain2D domain);
    [[nodiscard]] const MacGridDomain2D& domain() const noexcept { return domain_; }
    [[nodiscard]] int resolution() const noexcept { return domain_.resolution(); }
    [[nodiscard]] float boxSize() const noexcept { return domain_.boxSize(); }
    [[nodiscard]] float cellSize() const noexcept { return domain_.cellSize(); }
    [[nodiscard]] std::size_t cellCount() const noexcept { return density_.size(); }
    [[nodiscard]] std::size_t idC(int x, int y) const noexcept { return domain_.idC(x, y); }
    [[nodiscard]] std::size_t idX(int x, int y) const noexcept { return domain_.idX(x, y); }
    [[nodiscard]] std::size_t idY(int x, int y) const noexcept { return domain_.idY(x, y); }
    [[nodiscard]] std::span<float> density() noexcept { return density_; }
    [[nodiscard]] std::span<const float> density() const noexcept { return density_; }
    [[nodiscard]] std::span<float> velocityX() noexcept { return u_; }
    [[nodiscard]] std::span<const float> velocityX() const noexcept { return u_; }
    [[nodiscard]] std::span<float> velocityY() noexcept { return v_; }
    [[nodiscard]] std::span<const float> velocityY() const noexcept { return v_; }
    void reset() noexcept;
    void setBoxSize(float size);
private:
    MacGridDomain2D domain_;
    std::vector<float> density_, u_, v_;
};
