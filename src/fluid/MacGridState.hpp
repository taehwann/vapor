#pragma once

#include "domain/MacGridDomain.hpp"

#include <cstddef>
#include <span>
#include <vector>

// Physical state only: no OpenGL, UI, solver policy, or algorithm scratch memory.
// This first implementation deliberately represents a cubic, uniformly spaced grid.
class MacGridState {
public:
    explicit MacGridState(int resolution, float boxSize = 4.5f);
    explicit MacGridState(MacGridDomain domain);

    [[nodiscard]] const MacGridDomain& domain() const noexcept { return domain_; }
    [[nodiscard]] int resolution() const noexcept { return domain_.resolution(); }
    [[nodiscard]] float boxSize() const noexcept { return domain_.boxSize(); }
    [[nodiscard]] float cellSize() const noexcept { return domain_.cellSize(); }
    [[nodiscard]] std::size_t cellCount() const noexcept { return density_.size(); }

    [[nodiscard]] std::size_t idC(int x, int y, int z) const noexcept {
        return domain_.idC(x, y, z);
    }
    [[nodiscard]] std::size_t idX(int x, int y, int z) const noexcept {
        return domain_.idX(x, y, z);
    }
    [[nodiscard]] std::size_t idY(int x, int y, int z) const noexcept {
        return domain_.idY(x, y, z);
    }
    [[nodiscard]] std::size_t idZ(int x, int y, int z) const noexcept {
        return idC(x, y, z);
    }

    [[nodiscard]] std::span<float> density() noexcept { return density_; }
    [[nodiscard]] std::span<const float> density() const noexcept { return density_; }
    [[nodiscard]] std::span<float> velocityX() noexcept { return vx_; }
    [[nodiscard]] std::span<const float> velocityX() const noexcept { return vx_; }
    [[nodiscard]] std::span<float> velocityY() noexcept { return vy_; }
    [[nodiscard]] std::span<const float> velocityY() const noexcept { return vy_; }
    [[nodiscard]] std::span<float> velocityZ() noexcept { return vz_; }
    [[nodiscard]] std::span<const float> velocityZ() const noexcept { return vz_; }

    void reset() noexcept;
    // A domain-size change starts a fresh state, as in the original application.
    void setBoxSize(float boxSize);

private:
    MacGridDomain domain_;
    std::vector<float> density_, vx_, vy_, vz_;
};
