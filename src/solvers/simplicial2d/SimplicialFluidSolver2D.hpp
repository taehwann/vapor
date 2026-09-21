#pragma once

#include "common/FluidDiagnostics.hpp"
#include "solvers/simplicial2d/SimplicialFluidState2D.hpp"
#include "common/LinearSolveOptions.hpp"

#include <cstdint>
#include <span>
#include <vector>

struct SimplicialFluidParameters2D {
    float sourceStrength = 4.f;
    float emitterRadius = .08f;
    float emitterCenterS = .5f;
    float emitterCenterT = .2f;
    float buoyancy = 3.f;
    float smokeDecay = .08f;
    int cgIterations = 600;
    double cgAbsoluteTolerance = 1e-10;
    double cgRelativeTolerance = 1e-9;
    float advectionCfl = .5f;
    float stirStrength = 0.f; // Emitter perturbation velocity magnitude (world units/sec); 0 disables.
    std::uint32_t randomSeed = 0x9E3779B9u; // Applied at reset.
};

// Two-dimensional reduction of Elcott et al.: U is a primal edge flux,
// Omega=d0^T*star1*U is dual-cell circulation, and U=d0*Phi is recovered by CG.
class SimplicialFluidSolver2D final {
public:
    explicit SimplicialFluidSolver2D(int resolution = 64, float boxSize = 4.f);
    explicit SimplicialFluidSolver2D(SimplicialFluidState2D state);

    void advance(float dt);
    void resetState();
    [[nodiscard]] const SimplicialMesh2D& domain() const noexcept { return state_.domain(); }
    [[nodiscard]] FluidDiagnostics diagnostics() const;
    [[nodiscard]] float boxSize() const noexcept { return state_.boxSize(); }
    void setBoxSize(float size);

    [[nodiscard]] SimplicialFluidState2D& state() noexcept { return state_; }
    [[nodiscard]] const SimplicialFluidState2D& state() const noexcept { return state_; }
    [[nodiscard]] SimplicialFluidParameters2D& parameters() noexcept { return parameters_; }
    [[nodiscard]] const SimplicialFluidParameters2D& parameters() const noexcept { return parameters_; }
    [[nodiscard]] const LinearSolveResult& lastRecovery() const noexcept { return lastRecovery_; }
    [[nodiscard]] int lastSubstepCount() const noexcept { return lastSubstepCount_; }

    // Numerical hooks for manufactured DEC and circulation tests.
    void synchronizeFromStreamFunction();
    void advectVorticity(float dt);
    void recoverFluxFromVorticity();

private:
    [[nodiscard]] SimplicialPoint2D traceBack(
        SimplicialPoint2D point, float dt,
        std::span<const SimplicialPoint2D> triangleVelocity) const noexcept;
    float randomUnit() noexcept;
    void emitDensity();
    void stirEmitter(float dt);
    void addBuoyancy(float dt);
    void advectDensity(float dt);
    void advanceSubstep(float dt);

    SimplicialFluidState2D state_;
    SimplicialFluidParameters2D parameters_;
    std::vector<SimplicialPoint2D> triangleVelocity_, backtracked_, backtrackedVelocity_;
    std::vector<double> forceFlux_, forceVorticity_;
    std::vector<float> densitySource_, densityNext_;
    LinearSolveResult lastRecovery_;
    std::uint32_t randomState_ = 0x9E3779B9u;
    int lastSubstepCount_ = 0;
};
