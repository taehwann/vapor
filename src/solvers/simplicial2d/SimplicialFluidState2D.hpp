#pragma once

#include "solvers/simplicial2d/SimplicialMesh2D.hpp"

#include <span>
#include <vector>

// Integrated DEC fields: flux is a primal 1-form and vorticity is the dual
// 2-form associated with each primal vertex. Dye is a passive vertex scalar.
class SimplicialFluidState2D {
public:
    explicit SimplicialFluidState2D(int resolution = 32, float boxSize = 4.f);
    explicit SimplicialFluidState2D(SimplicialMesh2D mesh);

    [[nodiscard]] const SimplicialMesh2D& domain() const noexcept { return mesh_; }
    [[nodiscard]] int resolution() const noexcept { return mesh_.resolution(); }
    [[nodiscard]] float boxSize() const noexcept { return mesh_.boxSize(); }
    [[nodiscard]] std::span<float> density() noexcept { return density_; }
    [[nodiscard]] std::span<const float> density() const noexcept { return density_; }
    [[nodiscard]] std::span<double> edgeFlux() noexcept { return flux_; }
    [[nodiscard]] std::span<const double> edgeFlux() const noexcept { return flux_; }
    [[nodiscard]] std::span<double> vorticity() noexcept { return vorticity_; }
    [[nodiscard]] std::span<const double> vorticity() const noexcept { return vorticity_; }
    [[nodiscard]] std::span<double> streamFunction() noexcept { return streamFunction_; }
    [[nodiscard]] std::span<const double> streamFunction() const noexcept { return streamFunction_; }

    void reset() noexcept;
    void setBoxSize(float size);

private:
    SimplicialMesh2D mesh_;
    std::vector<float> density_;
    std::vector<double> flux_, vorticity_, streamFunction_;
};
