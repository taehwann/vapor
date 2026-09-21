#include "solvers/mac2d/MacGridDomain2D.hpp"
#include <cmath>
#include <limits>
#include <stdexcept>

MacGridDomain2D::MacGridDomain2D(int resolution, float boxSize)
    : resolution_(resolution), boxSize_(boxSize) {
    if (resolution <= 0) throw std::invalid_argument("2D MAC resolution must be positive");
    const auto n = static_cast<std::size_t>(resolution);
    if (n > static_cast<std::size_t>(std::numeric_limits<int>::max()) / (n + 1))
        throw std::length_error("2D MAC grid exceeds the supported index range");
    setBoxSize(boxSize);
}

void MacGridDomain2D::setBoxSize(float size) {
    if (!std::isfinite(size) || size <= 0.f || size / resolution_ <= 0.f ||
        !std::isfinite(1.f / (size / resolution_)))
        throw std::invalid_argument("2D MAC requires finite positive size and spacing");
    boxSize_ = size;
}
