#pragma once

// Solver controls/diagnostics shared across discretizations; no matrix layout.
struct LinearSolveOptions {
    int maxIterations = 200;
    int residualCheckInterval = 8;
    double absoluteTolerance = 1e-6;
    double relativeTolerance = 1e-5;
    bool fixedIterations = false;
};

struct LinearSolveResult {
    int iterations = 0;
    double initialResidual = 0.0;
    double finalResidual = 0.0;
    bool converged = false;
};
