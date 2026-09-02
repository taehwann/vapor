#pragma once

#include "numerics/LinearSolveOptions.hpp"

struct ProjectionOptions {
    float dt = 1.f / 60.f;
    float relaxation = 1.95f;
    LinearSolveOptions linearSolve;
};

struct ProjectionResult {
    LinearSolveResult linearSolve;
    bool residualAvailable = false;
    float divergenceBefore = 0.f;
    float divergenceAfter = 0.f;
    double removedRhsMean = 0.0;
};

// Full incompressibility projection: divergence, pressure solve, and velocity
// correction. State is a discretization-specific type: MAC fields today, or
// primal/dual cochains in a future mesh solver. No cast to a common fake grid.
template<class State>
class IPressureSolver {
public:
    virtual ~IPressureSolver() = default;
    virtual ProjectionResult project(State& state, const ProjectionOptions& options) = 0;
    virtual void reset() noexcept = 0;
};
