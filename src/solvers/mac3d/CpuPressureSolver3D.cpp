#include "solvers/mac3d/CpuPressureSolver3D.hpp"
#include "solvers/mac3d/MacGridOperators3D.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>

void CpuPressureSolver3D::prepareSystem(int n) {
    if (system_.resolution == n) return;
    system_.resize(n);
    divergence_.assign(system_.x.size(), 0.f);
    // A = -h^2 Laplacian with homogeneous Neumann walls. Boundary rows have
    // fewer neighbors. Constants are its null space; no arbitrary pin is needed
    // by SOR, and project() chooses a zero-mean representative after solving.
    for (int z = 0; z < n; ++z) for (int y = 0; y < n; ++y) for (int x = 0; x < n; ++x) {
        auto& a = system_.matrix[system_.index(x, y, z)];
        a.diagonal = float((x > 0) + (x + 1 < n) + (y > 0) + (y + 1 < n) + (z > 0) + (z + 1 < n));
        a.right = x + 1 < n ? -1.f : 0.f;
        a.up = y + 1 < n ? -1.f : 0.f;
        a.front = z + 1 < n ? -1.f : 0.f;
    }
}

ProjectionResult CpuPressureSolver3D::project(MacGridState3D& state, const ProjectionOptions& options) {
    sor_.setOmega(options.relaxation);
    if (!std::isfinite(options.dt) || options.dt <= 0.f) {
        throw std::invalid_argument("Projection time step must be finite and positive");
    }
    const float scale = state.cellSize() * state.cellSize() / options.dt;
    if (!std::isfinite(scale) || scale <= 0.f) {
        throw std::invalid_argument("Projection scale is outside the supported float range");
    }
    prepareSystem(state.resolution());
    MacGridOperators3D::enforceClosedBoundaries(state);
    MacGridOperators3D::computeDivergence(state, divergence_);
    ProjectionResult result;
    for (std::size_t i = 0; i < divergence_.size(); ++i) {
        if (!std::isfinite(divergence_[i])) throw std::invalid_argument("Non-finite divergence");
        result.divergenceBefore = std::max(result.divergenceBefore, std::abs(divergence_[i]));
        system_.b[i] = -divergence_[i] * scale;
    }
    // Closed walls imply a compatible (zero-sum) RHS. Remove floating-point drift.
    result.removedRhsMean = std::accumulate(system_.b.begin(), system_.b.end(), 0.0) / system_.b.size();
    for (auto& value : system_.b) value -= static_cast<float>(result.removedRhsMean);
    std::fill(system_.x.begin(), system_.x.end(), 0.f);
    result.linearSolve = sor_.solve(system_, options.linearSolve);
    result.residualAvailable = true;

    const double pressureMean = std::accumulate(system_.x.begin(), system_.x.end(), 0.0) / system_.x.size();
    for (auto& value : system_.x) value -= static_cast<float>(pressureMean);
    result.linearSolve.finalResidual = system_.residualInfinityNorm();
    const double target = std::max(options.linearSolve.absoluteTolerance,
                                  options.linearSolve.relativeTolerance * result.linearSolve.initialResidual);
    result.linearSolve.converged = result.linearSolve.finalResidual <= target;
    MacGridOperators3D::applyPressureGradient(state, system_.x, options.dt);
    MacGridOperators3D::enforceClosedBoundaries(state);
    result.divergenceAfter = MacGridOperators3D::maxDivergence(state);
    return result;
}

void CpuPressureSolver3D::reset() noexcept {
    std::fill(system_.x.begin(), system_.x.end(), 0.f);
    std::fill(system_.b.begin(), system_.b.end(), 0.f);
    std::fill(divergence_.begin(), divergence_.end(), 0.f);
}
