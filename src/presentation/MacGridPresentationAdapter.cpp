#include "presentation/MacGridPresentationAdapter.hpp"
#include "fluid/MacGridState.hpp"

MacGridPresentationAdapter::MacGridPresentationAdapter(const MacGridState& state) noexcept
    : state_(state) {}

RenderData MacGridPresentationAdapter::renderData() const noexcept {
    const int resolution = state_.resolution();
    return {state_.density(), resolution, resolution, resolution, state_.boxSize()};
}
