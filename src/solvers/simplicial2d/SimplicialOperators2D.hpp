#pragma once

#include "solvers/simplicial2d/SimplicialMesh2D.hpp"
#include "common/LinearSolveOptions.hpp"

#include <span>
#include <vector>

class SimplicialOperators2D {
public:
    static void fluxFromStreamFunction(const SimplicialMesh2D& mesh,
                                       std::span<const double> streamFunction,
                                       std::span<double> flux);
    static void vorticityFromFlux(const SimplicialMesh2D& mesh,
                                  std::span<const double> flux,
                                  std::span<double> vorticity);
    [[nodiscard]] static std::vector<SimplicialPoint2D> reconstructTriangleVelocities(
        const SimplicialMesh2D& mesh, std::span<const double> flux);
    [[nodiscard]] static SimplicialPoint2D sampleVelocity(
        const SimplicialMesh2D& mesh, std::span<const SimplicialPoint2D> triangleVelocity,
        SimplicialPoint2D point) noexcept;
    [[nodiscard]] static float sampleVertexScalar(
        const SimplicialMesh2D& mesh, std::span<const float> values,
        SimplicialPoint2D point) noexcept;
    [[nodiscard]] static double maxDivergence(
        const SimplicialMesh2D& mesh, std::span<const double> flux);
    [[nodiscard]] static LinearSolveResult recoverFlux(
        const SimplicialMesh2D& mesh, std::span<const double> vorticity,
        std::span<double> streamFunction, std::span<double> flux,
        const LinearSolveOptions& options);
};
