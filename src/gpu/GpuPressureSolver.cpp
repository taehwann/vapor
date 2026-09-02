#include "gpu/GpuPressureSolver.hpp"
#include "fluid/MacGridOperators.hpp"

#include <cmath>
#include <stdexcept>

ProjectionResult GpuPressureSolver::project(MacGridState& state, const ProjectionOptions& options) {
    backend_.validateGrid(state);
    const auto& solve = options.linearSolve;
    const float scale = state.cellSize() * state.cellSize() / options.dt;
    if (!std::isfinite(options.dt) || options.dt <= 0.f ||
        !std::isfinite(scale) || scale <= 0.f ||
        !std::isfinite(options.relaxation) || options.relaxation <= 0.f || options.relaxation >= 2.f ||
        solve.maxIterations < 0 || solve.residualCheckInterval <= 0 ||
        !std::isfinite(solve.absoluteTolerance) || solve.absoluteTolerance < 0.0 ||
        !std::isfinite(solve.relativeTolerance) || solve.relativeTolerance < 0.0) {
        throw std::invalid_argument("Invalid GPU pressure options");
    }
    MacGridOperators::enforceClosedBoundaries(state);
    ProjectionResult result;
    result.divergenceBefore = MacGridOperators::maxDivergence(state);
    backend_.project(options.dt, solve.maxIterations, options.relaxation,
        1.f / state.cellSize(), state.cellSize(),
        state.velocityX().data(), state.velocityY().data(), state.velocityZ().data(),
        state.velocityX().data(), state.velocityY().data(), state.velocityZ().data());
    MacGridOperators::enforceClosedBoundaries(state);
    result.divergenceAfter = MacGridOperators::maxDivergence(state);
    // The existing GPU algorithm uses a fixed iteration budget. Do not report
    // a fabricated residual or convergence decision through the common API.
    result.linearSolve.iterations = solve.maxIterations;
    return result;
}
