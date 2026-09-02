#pragma once

#include "domain/IDomain.hpp"

struct FluidDiagnostics {
    float kineticEnergy = 0.f;
    float maxVelocity = 0.f;
    float maxDivergence = 0.f;
};

// Controls shared by application-facing fluid solvers. A concrete solver may
// ignore optional controls that do not apply to its algorithm.
struct FluidSolverControls {
    float sourceStrength = 0.f;
    float emitterRadius = 0.f;
    float buoyancy = 0.f;
    float smokeDecay = 0.f;
    int pressureIterations = 0;
    float pressureRelaxation = 1.f;
    float stirStrength = 0.f;
    float advectionCfl = 1.f;
};

struct PressureSolveDiagnostics {
    bool residualAvailable = false;
    bool converged = false;
    int iterations = 0;
    double finalResidual = 0.0;
};

// The app advances and inspects a simulation without knowing its discretization
// or compute backend. Concrete state remains private to registered factories and
// their presentation adapters.
class IFluidSolver {
public:
    virtual ~IFluidSolver() = default;
    virtual void advance(float dt) = 0;
    virtual void resetState() = 0;
    [[nodiscard]] virtual const IDomain& domain() const noexcept = 0;
    [[nodiscard]] virtual FluidDiagnostics diagnostics() const = 0;
    [[nodiscard]] virtual FluidSolverControls controls() const noexcept = 0;
    virtual void setControls(const FluidSolverControls& controls) = 0;
    [[nodiscard]] virtual PressureSolveDiagnostics pressureDiagnostics() const noexcept = 0;
    [[nodiscard]] virtual float boxSize() const noexcept = 0;
    virtual void setBoxSize(float size) = 0;
};
