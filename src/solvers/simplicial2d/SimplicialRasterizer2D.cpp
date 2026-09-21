#include "solvers/simplicial2d/SimplicialRasterizer2D.hpp"

#include "solvers/simplicial2d/SimplicialFluidState2D.hpp"

#include <algorithm>
#include <cmath>

SimplicialRasterizer2D::SimplicialRasterizer2D(
    const SimplicialFluidState2D& state)
    : state_(state), height_(state.domain().isTeapot() ? 512 : state.resolution() + 1),
      width_(height_),
      raster_(static_cast<std::size_t>(width_) * height_) {}

ImageRenderData SimplicialRasterizer2D::renderData() const noexcept {
    const auto& mesh = state_.domain();
    const double box = state_.boxSize();
    const double worldWidth = box;
    const double worldHeight = box;
    const auto density = state_.density();
    std::fill(raster_.begin(), raster_.end(), 0.f);

    for (int row = 0; row < height_; ++row) {
        const double y = (row + 0.5) * worldHeight / height_;
        for (int column = 0; column < width_; ++column) {
            const double x = (column + 0.5) * worldWidth / width_;
            const SimplicialPoint2D point{x, y};
            const int triangleIndex = mesh.locateTriangle(point);
            if (triangleIndex < 0) continue;
            const auto barycentric = mesh.barycentric(triangleIndex, point);
            const auto& triangle = mesh.triangles()[triangleIndex];
            raster_[static_cast<std::size_t>(column) + static_cast<std::size_t>(width_) * row] =
                static_cast<float>(
                    barycentric[0] * density[triangle.vertices[0]] +
                    barycentric[1] * density[triangle.vertices[1]] +
                    barycentric[2] * density[triangle.vertices[2]]);
        }
    }
    return {raster_, width_, height_, static_cast<float>(worldWidth),
            static_cast<float>(worldHeight), true};
}
