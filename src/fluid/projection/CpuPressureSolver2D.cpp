#include "fluid/projection/CpuPressureSolver2D.hpp"
#include "fluid/MacGridOperators2D.hpp"
#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>

ProjectionResult CpuPressureSolver2D::project(MacGridState2D& s, const ProjectionOptions& options) {
    const auto& control = options.linearSolve;
    const float scale = s.cellSize() * s.cellSize() / options.dt;
    if (!std::isfinite(options.dt) || options.dt <= 0.f || !std::isfinite(scale) || scale <= 0.f ||
        !std::isfinite(options.relaxation) || options.relaxation <= 0.f || options.relaxation >= 2.f ||
        control.maxIterations < 0 || control.residualCheckInterval <= 0 ||
        !std::isfinite(control.absoluteTolerance) || control.absoluteTolerance < 0.0 ||
        !std::isfinite(control.relativeTolerance) || control.relativeTolerance < 0.0)
        throw std::invalid_argument("Invalid 2D SOR projection options");
    for (auto field : {s.velocityX(), s.velocityY()})
        for (float value : field) if (!std::isfinite(value)) throw std::invalid_argument("Non-finite 2D velocity");

    pressure_.assign(s.cellCount(), 0.f);
    rhs_.resize(s.cellCount());
    MacGridOperators2D::enforceClosedBoundaries(s);
    MacGridOperators2D::computeDivergence(s, rhs_);
    ProjectionResult result;
    for (float& value : rhs_) {
        result.divergenceBefore = std::max(result.divergenceBefore, std::abs(value));
        value *= -scale; // A = -h^2 Laplacian; A p = -h^2 div(u*) / dt.
        if (!std::isfinite(value)) throw std::invalid_argument("2D pressure RHS overflow");
    }
    // Closed-wall flux telescopes to zero. Remove roundoff from the null space.
    result.removedRhsMean = std::accumulate(rhs_.begin(), rhs_.end(), 0.0) / rhs_.size();
    for (float& value : rhs_) value -= float(result.removedRhsMean);
    const int n = s.resolution();
    auto neighbors = [&](int x, int y) {
        double sum = 0.0;
        if (x > 0) sum += pressure_[s.idC(x - 1, y)];
        if (x + 1 < n) sum += pressure_[s.idC(x + 1, y)];
        if (y > 0) sum += pressure_[s.idC(x, y - 1)];
        if (y + 1 < n) sum += pressure_[s.idC(x, y + 1)];
        return sum;
    };
    auto diagonal = [n](int x, int y) { return (x > 0) + (x + 1 < n) + (y > 0) + (y + 1 < n); };
    auto residual = [&] {
        double value = 0.0;
        for (int y = 0; y < n; ++y) for (int x = 0; x < n; ++x) {
            const auto i = s.idC(x, y);
            value = std::max(value, std::abs(double(rhs_[i]) + neighbors(x, y) - double(diagonal(x, y)) * pressure_[i]));
        }
        return value;
    };
    auto& solve = result.linearSolve;
    solve.initialResidual = solve.finalResidual = residual();
    const double target = std::max(control.absoluteTolerance, control.relativeTolerance * solve.initialResidual);
    solve.converged = solve.finalResidual <= target;
    for (int iter = 0; iter < control.maxIterations && (control.fixedIterations || !solve.converged); ++iter) {
        for (int parity = 0; parity < 2; ++parity) {
            for (int y = 0; y < n; ++y) for (int x = (y + parity) & 1; x < n; x += 2) {
                const int d = diagonal(x, y);
                if (d == 0) continue; // The single-cell closed domain has no pressure gradient.
                const auto i = s.idC(x, y);
                pressure_[i] += options.relaxation * ((rhs_[i] + neighbors(x, y)) / d - pressure_[i]);
            }
        }
        solve.iterations = iter + 1;
        if (solve.iterations % control.residualCheckInterval == 0 || solve.iterations == control.maxIterations) {
            solve.finalResidual = residual();
            solve.converged = solve.finalResidual <= target;
        }
    }
    const double mean = std::accumulate(pressure_.begin(), pressure_.end(), 0.0) / pressure_.size();
    for (double& value : pressure_) {
        value -= mean;
        if (!std::isfinite(value)) throw std::runtime_error("2D SOR produced non-finite pressure");
    }
    solve.finalResidual = residual();
    solve.converged = solve.finalResidual <= target;
    result.residualAvailable = true;
    MacGridOperators2D::applyPressureGradient(s, pressure_, options.dt);
    MacGridOperators2D::enforceClosedBoundaries(s);
    result.divergenceAfter = MacGridOperators2D::maxDivergence(s);
    return result;
}

void CpuPressureSolver2D::reset() noexcept {
    std::fill(pressure_.begin(), pressure_.end(), 0.f);
    std::fill(rhs_.begin(), rhs_.end(), 0.f);
}
