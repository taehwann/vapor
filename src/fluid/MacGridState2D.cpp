#include "fluid/MacGridState2D.hpp"
#include <algorithm>
#include <utility>

MacGridState2D::MacGridState2D(int resolution, float boxSize)
    : MacGridState2D(MacGridDomain2D(resolution, boxSize)) {}

MacGridState2D::MacGridState2D(MacGridDomain2D domain) : domain_(std::move(domain)) {
    const auto n = static_cast<std::size_t>(domain_.resolution());
    density_.resize(n * n);
    u_.resize((n + 1) * n);
    v_.resize(n * (n + 1));
}
void MacGridState2D::reset() noexcept {
    for (auto field : {density(), velocityX(), velocityY()}) std::fill(field.begin(), field.end(), 0.f);
}
void MacGridState2D::setBoxSize(float size) {
    const float previous = boxSize();
    domain_.setBoxSize(size);
    if (size != previous) reset();
}
