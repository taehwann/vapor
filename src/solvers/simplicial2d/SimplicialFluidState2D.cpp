#include "solvers/simplicial2d/SimplicialFluidState2D.hpp"

#include <algorithm>
#include <utility>

SimplicialFluidState2D::SimplicialFluidState2D(int resolution, float boxSize)
    : SimplicialFluidState2D(SimplicialMesh2D(resolution, boxSize)) {}

SimplicialFluidState2D::SimplicialFluidState2D(SimplicialMesh2D mesh)
    : mesh_(std::move(mesh)), density_(mesh_.vertexCount()), flux_(mesh_.edgeCount()),
      vorticity_(mesh_.vertexCount()), streamFunction_(mesh_.vertexCount()) {}

void SimplicialFluidState2D::reset() noexcept {
    std::fill(density_.begin(), density_.end(), 0.f);
    std::fill(flux_.begin(), flux_.end(), 0.0);
    std::fill(vorticity_.begin(), vorticity_.end(), 0.0);
    std::fill(streamFunction_.begin(), streamFunction_.end(), 0.0);
}

void SimplicialFluidState2D::setBoxSize(float size) {
    const float previous = boxSize();
    mesh_.setBoxSize(size);
    if (previous != size) reset();
}
