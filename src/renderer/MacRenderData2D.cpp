#include "renderer/MacRenderData2D.hpp"
#include "solvers/mac2d/MacGridState2D.hpp"

ImageRenderData imageData(const MacGridState2D& state) noexcept {
    return {state.density(), state.resolution(), state.resolution(), state.boxSize(), state.boxSize()};
}
