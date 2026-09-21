#pragma once

#include "common/FluidDiagnostics.hpp"
#include "solvers/mac3d/MacGridState3D.hpp"
#include "solvers/mac3d/CpuAdvection.hpp"
#include <memory>
#include "solvers/mac3d/CpuPressureSolver3D.hpp"

#include <span>
#include <cstdint>
#include <vector>

// Parameters of the 3D MAC timestep.
struct MacGridParameters3D {
    MacGridParameters3D() = default;
    explicit MacGridParameters3D(MacAlgorithm algorithm);
    float emitterRadius = 0.035f;
    float emitterCenterX = 0.5f, emitterCenterY = 0.5f, emitterCenterZ = 0.1f;
    float emitterKickMultiplier = 1.f;
    float emitterKickBase = 0.8f;
    float emitterKickScale = 0.8f;
    float stirStrength = 0.5f;
    std::uint32_t randomSeed = 0x9E3779B9u; // Applied at construction and reset.
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

// Owns physical state, timestep workspace, and selected CPU/GPU stage implementations.
class MacGridFluidSolver3D final {
public:

    explicit MacGridFluidSolver3D(MacAlgorithm algorithm); // Ready-to-run CPU simulation.
    explicit MacGridFluidSolver3D(int res);
    explicit MacGridFluidSolver3D(MacGridState3D state);

    MacGridFluidSolver3D(const MacGridFluidSolver3D&) = delete;
    MacGridFluidSolver3D& operator=(const MacGridFluidSolver3D&) = delete;

    MacGridParameters3D& parameters() noexcept { return parameters_; }
    const MacGridParameters3D& parameters() const noexcept { return parameters_; }
    const MacGridState3D& state() const noexcept { return state_; }
    const MacGridDomain3D& domain() const noexcept { return state_.domain(); }

    // Ownership transfers to the solver. Changing a backend preserves physical fields.
    void setPressureSolver(std::unique_ptr<Projection3D> solver);
    void setAdvection(std::unique_ptr<Advection3D> method);
    void useDefaultPressureSolver();
    void useDefaultAdvection();

    void advance(float dt);

    void step(float dt, float buoyancy, float sourceStrength, int projectIterations);
    void stepNoReflect(float dt, float buoyancy, float sourceStrength, int projectIterations);
    void resetState();
    void setBoxSize(float size);

    [[nodiscard]] int resolution() const;
    [[nodiscard]] float boxSize() const noexcept;
    [[nodiscard]] const ProjectionResult& lastProjection() const;
    [[nodiscard]] FluidDiagnostics diagnostics() const;
    [[nodiscard]] float kineticEnergy() const;
    [[nodiscard]] float maxVelocity() const;
    [[nodiscard]] float computeDivergenceNorm() const;

private:
    float randomUnit();
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

    MacGridParameters3D parameters_;
    MacGridState3D state_;
    std::uint32_t randomState_ = parameters_.randomSeed;
    // Advection/reflection workspace is separate from physical state.
    AdvectionWorkspace advectionWorkspace_;
    std::vector<float> vx0, vy0, vz0;
    std::vector<float> vxTilde, vyTilde, vzTilde;
    std::vector<float> vxHat, vyHat, vzHat;
    std::vector<float> vxDiv0, vyDiv0, vzDiv0;
    CpuPressureSolver3D cpuPressure_{6.f};
    CpuAdvection cpuAdvection_;
    // Optional experimental stage overrides; normal CPU construction uses concrete members.
    std::unique_ptr<Projection3D> pressureSolver_;
    std::unique_ptr<Advection3D> advection_;
    ProjectionResult lastProjection_;
};
