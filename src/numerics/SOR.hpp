#pragma once

#include "numerics/ILinearSolver.hpp"

// Intended for positive-definite or compatible positive-semidefinite stencils
// such as the pressure Laplacian. Validation does not prove matrix definiteness.
class SOR final : public ILinearSolver {
public:
    // A nonzero diagonal floor damps rows with smaller diagonals. A floor of 6
    // reproduces Vapor's old lagged-Neumann-ghost relaxation on its unscaled grid
    // Laplacian; the default 0 uses ordinary red-black SOR.
    explicit SOR(float omega = 1.95f, float relaxationDiagonalFloor = 0.f);
    void setOmega(float omega);
    LinearSolveResult solve(StencilLinearSystem& system,
                            const LinearSolveOptions& options) const override;

private:
    float omega_;
    float relaxationDiagonalFloor_;
};
