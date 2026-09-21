#pragma once

#include "common/LinearSolveOptions.hpp"

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

