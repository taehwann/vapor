#pragma once
#include "common/FluidDiagnostics.hpp"
#include "solvers/mac2d/SemiLagrangian2D.hpp"
#include "solvers/mac2d/CpuPressureSolver2D.hpp"

struct MacGridParameters2D {
    bool reflection = false;
    float sourceStrength = 4.f;
    float emitterKickMultiplier = 1.f;
    float emitterRadius = .08f; // Fraction of box size.
    float emitterCenterX = .5f, emitterCenterY = .2f;
    float buoyancy = 3.f;
    float smokeDecay = .08f;
    float sorOmega = 1.9f;
    int projectIterations = 600;
};

// Simple or midpoint-reflection stepping on a closed 2D MAC grid.
class MacGridFluidSolver2D final {
public:
    MacGridFluidSolver2D(); // Ready-to-run simple semi-Lagrangian plume.
    explicit MacGridFluidSolver2D(int resolution, float boxSize = 4.f);
    explicit MacGridFluidSolver2D(MacGridState2D state);
    MacGridFluidSolver2D(const MacGridFluidSolver2D&) = delete;
    MacGridFluidSolver2D& operator=(const MacGridFluidSolver2D&) = delete;
    MacGridFluidSolver2D(MacGridFluidSolver2D&&) noexcept = default;
    MacGridFluidSolver2D& operator=(MacGridFluidSolver2D&&) noexcept = default;
    void advance(float dt);
    void resetState();
    [[nodiscard]] const MacGridDomain2D& domain() const noexcept { return state_.domain(); }
    [[nodiscard]] FluidDiagnostics diagnostics() const;
    [[nodiscard]] const MacGridState2D& state() const noexcept { return state_; }
    [[nodiscard]] MacGridParameters2D& parameters() noexcept { return parameters_; }
    [[nodiscard]] const MacGridParameters2D& parameters() const noexcept { return parameters_; }
    [[nodiscard]] const ProjectionResult& lastProjection() const noexcept { return lastProjection_; }
    [[nodiscard]] const ProjectionResult& midpointProjection() const noexcept { return midpointProjection_; }
    [[nodiscard]] float boxSize() const noexcept { return state_.boxSize(); }
    void setBoxSize(float size);
private:
    void advanceReflection(float dt);
    void advectVelocity(std::span<const float> u, std::span<const float> v, float dt);
    void advectDensity(float dt);
    void emit();
    void addBuoyancy(float dt);
    ProjectionResult project(float dt);
    MacGridState2D state_;
    MacGridParameters2D parameters_;
    CpuPressureSolver2D pressure_;
    SemiLagrangian2D advection_;
    std::vector<float> uNext_, vNext_, densityNext_, uHat_, vHat_;
    ProjectionResult lastProjection_, midpointProjection_;
};
