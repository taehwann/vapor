#include "fluid/MacGridState.hpp"

#include <algorithm>
#include <utility>

MacGridState::MacGridState(int resolution, float boxSize)
    : MacGridState(MacGridDomain(resolution, boxSize)) {}

MacGridState::MacGridState(MacGridDomain domain)
    : domain_(std::move(domain)) {
    const std::size_t n = static_cast<std::size_t>(domain_.resolution());
    density_.assign(n * n * n, 0.f);
    vx_.assign((n + 1) * n * n, 0.f);
    vy_.assign(n * (n + 1) * n, 0.f);
    vz_.assign(n * n * (n + 1), 0.f);
}

void MacGridState::reset() noexcept {
    std::fill(density_.begin(), density_.end(), 0.f);
    std::fill(vx_.begin(), vx_.end(), 0.f);
    std::fill(vy_.begin(), vy_.end(), 0.f);
    std::fill(vz_.begin(), vz_.end(), 0.f);
}

void MacGridState::setBoxSize(float boxSize) {
    const float previous = domain_.boxSize();
    domain_.setBoxSize(boxSize);
    if (boxSize != previous) reset();
}
