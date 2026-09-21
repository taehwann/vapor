#include "solvers/mac3d/MacGridState3D.hpp"

#include <algorithm>
#include <utility>

MacGridState3D::MacGridState3D(int resolution, float boxSize)
    : MacGridState3D(MacGridDomain3D(resolution, boxSize)) {}

MacGridState3D::MacGridState3D(MacGridDomain3D domain)
    : domain_(std::move(domain)) {
    const std::size_t n = static_cast<std::size_t>(domain_.resolution());
    density_.assign(n * n * n, 0.f);
    vx_.assign((n + 1) * n * n, 0.f);
    vy_.assign(n * (n + 1) * n, 0.f);
    vz_.assign(n * n * (n + 1), 0.f);
}

void MacGridState3D::reset() noexcept {
    std::fill(density_.begin(), density_.end(), 0.f);
    std::fill(vx_.begin(), vx_.end(), 0.f);
    std::fill(vy_.begin(), vy_.end(), 0.f);
    std::fill(vz_.begin(), vz_.end(), 0.f);
}

void MacGridState3D::setBoxSize(float boxSize) {
    const float previous = domain_.boxSize();
    domain_.setBoxSize(boxSize);
    if (boxSize != previous) reset();
}
