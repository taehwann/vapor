#pragma once

#include "fluid/IFluidSolver.hpp"
#include "fluid/MacGridState.hpp"
#include "fluid/advection/MacCormack.hpp"
#include "fluid/advection/SemiLagrangian.hpp"
#include "fluid/projection/CpuPressureSolver.hpp"

#include <span>
#include <vector>

// Algorithm-specific controls stay outside the generic IFluidSolver interface.
struct MacGridParameters {
    float emitterRadius = 0.035f;
    float emitterCenterX = 0.5f, emitterCenterY = 0.5f, emitterCenterZ = 0.1f;
    float emitterKickBase = 0.8f;
    float emitterKickScale = 0.8f;
    float stirStrength = 0.5f;
    float sorOmega = 1.95f;
    float buoyancySplit = 0.5f;
    float cflMc = 3.0f;
    bool macCormackSmoke = true, macCormackVel = true;

    // Dissipation tuned to fade out before reaching the top wall
    float smokeDecay = 0.06f;

    bool reflection = true;
    float buoyancy = 2.5f;
    float sourceStrength = 1.f;
    int projectIterations = 64;
};

// Coordinates MAC stages using injected pressure/advection implementations.
// Default CPU implementations are owned. Overrides are borrowed and must
// outlive their use by this solver; no OpenGL lifecycle belongs here.
class MacGridFluidSolver final : public IFluidSolver {
public:

    explicit MacGridFluidSolver(int res);
    explicit MacGridFluidSolver(MacGridState state);

    MacGridFluidSolver(const MacGridFluidSolver&) = delete;
    MacGridFluidSolver& operator=(const MacGridFluidSolver&) = delete;

    MacGridParameters& parameters() noexcept { return parameters_; }
    const MacGridParameters& parameters() const noexcept { return parameters_; }
    const MacGridState& state() const noexcept { return state_; }
    const IDomain& domain() const noexcept override { return state_.domain(); }

    void setPressureSolver(IPressureSolver<MacGridState>& solver) noexcept;
    void useDefaultPressureSolver() noexcept;
    void setVelocityAdvection(const IAdvection<MacAdvectionOperation>& method) noexcept;
    void setDensityAdvection(const IAdvection<MacAdvectionOperation>& method) noexcept;
    void setVelocityAdvection(const IAdvection<MacAdvectionOperation>&&) = delete;
    void setDensityAdvection(const IAdvection<MacAdvectionOperation>&&) = delete;
    void useDefaultAdvection() noexcept;

    void advance(float dt) override;

    void step(float dt, float buoyancy, float sourceStrength, int projectIterations);
    void stepNoReflect(float dt, float buoyancy, float sourceStrength, int projectIterations);
    void resetState() override;
    void setBoxSize(float size) override;

    [[nodiscard]] int resolution() const;
    [[nodiscard]] float boxSize() const noexcept override;
    [[nodiscard]] const ProjectionResult& lastProjection() const;
    [[nodiscard]] FluidDiagnostics diagnostics() const override;
    [[nodiscard]] FluidSolverControls controls() const noexcept override;
    void setControls(const FluidSolverControls& controls) override;
    [[nodiscard]] PressureSolveDiagnostics pressureDiagnostics() const noexcept override;
    [[nodiscard]] float kineticEnergy() const;
    [[nodiscard]] float maxVelocity() const;
    [[nodiscard]] float computeDivergenceNorm() const;

private:
    void allocateBuffers();
    int idC(int x, int y, int z) const;
    int idX(int x, int y, int z) const;
    int idY(int x, int y, int z) const;
    int idZ(int x, int y, int z) const;
    float h() const;

    void advectScalar(std::span<float> field, float dt);

    void advectVx(float dt, std::span<const float> vxSrc,
        std::span<const float> vxVel, std::span<const float> vyVel, std::span<const float> vzVel,
        std::span<float> out);

    void advectVy(float dt, std::span<const float> vySrc,
        std::span<const float> vxVel, std::span<const float> vyVel, std::span<const float> vzVel,
        std::span<float> out);

    void advectVz(float dt, std::span<const float> vzSrc,
        std::span<const float> vxVel, std::span<const float> vyVel, std::span<const float> vzVel,
        std::span<float> out);

    void enforceVelocityBoundaries();

    void project(float dt, int iterations);

    void emit(float strength, float dt);

    void applyBuoyancy(float dtStep, float buoyancyVal);

    // Smoke density dissipation does not change velocity.
    void applyDissipation(float dt);

    [[nodiscard]] const IAdvection<MacAdvectionOperation>& velocityAdvection() const;

    [[nodiscard]] const IAdvection<MacAdvectionOperation>& densityAdvection() const;

    MacGridParameters parameters_;
    MacGridState state_;
    // Advection/reflection workspace is separate from physical state.
    AdvectionWorkspace advectionWorkspace_;
    std::vector<float> vx0, vy0, vz0;
    std::vector<float> vxTilde, vyTilde, vzTilde;
    std::vector<float> vxHat, vyHat, vzHat;
    std::vector<float> vxDiv0, vyDiv0, vzDiv0;
    CpuPressureSolver defaultPressure_{6.f};
    IPressureSolver<MacGridState>* pressureSolver_ = &defaultPressure_;
    ProjectionResult lastProjection_;
    SemiLagrangian semiLagrangian;
    MacCormack macCormack;
    const IAdvection<MacAdvectionOperation>* velocityAdvection_ = nullptr;
    const IAdvection<MacAdvectionOperation>* densityAdvection_ = nullptr;
};
