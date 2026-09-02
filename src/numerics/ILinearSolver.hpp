#pragma once

#include "numerics/StencilLinearSystem.hpp"
#include "numerics/LinearSolveOptions.hpp"

class ILinearSolver {
public:
    virtual ~ILinearSolver() = default;
    virtual LinearSolveResult solve(StencilLinearSystem& system,
                                   const LinearSolveOptions& options) const = 0;
};
