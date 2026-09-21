#pragma once

#include "renderer/ImageRenderData.hpp"

#include <vector>

class SimplicialFluidState2D;

// Resamples vertex dye from the square triangular mesh into an orthogonal
// world-space raster so physical directions are displayed without shear.
class SimplicialRasterizer2D final {
public:
    explicit SimplicialRasterizer2D(const SimplicialFluidState2D& state);
    SimplicialRasterizer2D(const SimplicialFluidState2D&&) = delete;
    [[nodiscard]] ImageRenderData renderData() const noexcept;
private:
    const SimplicialFluidState2D& state_;
    int height_ = 0;
    int width_ = 0;
    mutable std::vector<float> raster_;
};
