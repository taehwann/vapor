#include "numerics/SOR.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

SOR::SOR(float omega, float relaxationDiagonalFloor)
    : omega_(omega), relaxationDiagonalFloor_(relaxationDiagonalFloor) {
    setOmega(omega);
    if (!std::isfinite(relaxationDiagonalFloor) || relaxationDiagonalFloor < 0.f) {
        throw std::invalid_argument("SOR diagonal floor must be finite and nonnegative");
    }
}

void SOR::setOmega(float omega) {
    if (!std::isfinite(omega) || omega <= 0.f || omega >= 2.f) {
        throw std::invalid_argument("SOR relaxation factor must be between 0 and 2");
    }
    omega_ = omega;
}

LinearSolveResult SOR::solve(StencilLinearSystem& s, const LinearSolveOptions& options) const {
    if (options.maxIterations < 0 || options.residualCheckInterval <= 0 ||
        !std::isfinite(options.absoluteTolerance) || options.absoluteTolerance < 0.0 ||
        !std::isfinite(options.relativeTolerance) || options.relativeTolerance < 0.0) {
        throw std::invalid_argument("Invalid linear solve options");
    }
    LinearSolveResult result;
    result.initialResidual = result.finalResidual = s.residualInfinityNorm();
    const double target = std::max(options.absoluteTolerance,
                                  options.relativeTolerance * result.initialResidual);
    result.converged = result.finalResidual <= target;
    if (result.converged && !options.fixedIterations) return result;

    const int n = s.resolution;
    for (int iter = 0; iter < options.maxIterations; ++iter) {
        for (int parity = 0; parity < 2; ++parity) {
            for (int k = 0; k < n; ++k) for (int j = 0; j < n; ++j) {
                for (int i = (j + k + parity) & 1; i < n; i += 2) {
                    const auto idx = s.index(i, j, k);
                    const auto& a = s.matrix[idx];
                    if (a.diagonal == 0.f) continue;
                    float offDiagonal = 0.f;
                    if (i > 0) offDiagonal += s.matrix[idx - 1].right * s.x[idx - 1];
                    if (i + 1 < n) offDiagonal += a.right * s.x[idx + 1];
                    if (j > 0) offDiagonal += s.matrix[idx - n].up * s.x[idx - n];
                    if (j + 1 < n) offDiagonal += a.up * s.x[idx + n];
                    if (k > 0) offDiagonal += s.matrix[idx - n * n].front * s.x[idx - n * n];
                    if (k + 1 < n) offDiagonal += a.front * s.x[idx + n * n];
                    const float diagonal = std::max(a.diagonal, relaxationDiagonalFloor_);
                    s.x[idx] += omega_ * (s.b[idx] - offDiagonal - a.diagonal * s.x[idx]) / diagonal;
                }
            }
        }
        result.iterations = iter + 1;
        if (result.iterations % options.residualCheckInterval == 0 ||
            result.iterations == options.maxIterations) {
            result.finalResidual = s.residualInfinityNorm();
            result.converged = result.finalResidual <= target;
            if (result.converged && !options.fixedIterations) break;
        }
    }
    return result;
}
