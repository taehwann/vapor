#include "renderer/MacRenderData.hpp"
#include "solvers/mac3d/MacGridState3D.hpp"

RenderData volumeData(const MacGridState3D& state) noexcept {
    return {state.density(), state.resolution(), state.resolution(), state.resolution(), state.boxSize()};
}
