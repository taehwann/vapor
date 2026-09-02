#include "presentation/MacGridPresentationAdapter2D.hpp"
#include "fluid/MacGridState2D.hpp"

ImageRenderData MacGridPresentationAdapter2D::renderData() const noexcept {
    return {state_.density(), state_.resolution(), state_.resolution(), state_.boxSize(), state_.boxSize()};
}
