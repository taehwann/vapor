#include "domain/MacGridDomain.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>

MacGridDomain::MacGridDomain(int resolution, float boxSize)
    : resolution_(resolution), boxSize_(boxSize) {
    if (resolution <= 0) throw std::invalid_argument("MAC grid resolution must be positive");
    const std::size_t n = static_cast<std::size_t>(resolution);
    const auto limit = static_cast<std::size_t>(std::numeric_limits<int>::max());
    if (n > limit / n || n * n > limit / (n + 1)) {
        throw std::length_error("MAC grid exceeds the supported index range");
    }
    setBoxSize(boxSize);
}

void MacGridDomain::setBoxSize(float size) {
    if (!std::isfinite(size) || size <= 0.f || size / resolution_ <= 0.f) {
        throw std::invalid_argument("MAC grid requires a finite positive box size and spacing");
    }
    boxSize_ = size;
}
