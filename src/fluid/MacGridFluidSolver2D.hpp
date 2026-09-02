#pragma once
#include "fluid/IFluidSolver.hpp"
#include "fluid/advection/SemiLagrangian2D.hpp"
#include "fluid/projection/CpuPressureSolver2D.hpp"

struct MacGridParameters2D {
    float sourceStrength = 4.f;
    float emitterRadius = .08f; // Fraction of box size.
    float emitterCenterX = .5f, emitterCenterY = .12f;
    float buoyancy = 3.f;
    float smokeDecay = .08f;
    float sorOmega = 1.9f;
    int projectIterations = 600;
};

// One first-order step: emit, advect velocity, apply buoyancy, project, advect density.
class MacGridFluidSolver2D final : public IFluidSolver {
public:
    explicit MacGridFluidSolver2D(int resolution = 64, float boxSize = 4.f);
    explicit MacGridFluidSolver2D(MacGridState2D state);
    MacGridFluidSolver2D(const MacGridFluidSolver2D&) = delete;
    MacGridFluidSolver2D& operator=(const MacGridFluidSolver2D&) = delete;
    void advance(float dt) override;
    void resetState() override;
    [[nodiscard]] const IDomain& domain() const noexcept override { return state_.domain(); }
    [[nodiscard]] FluidDiagnostics diagnostics() const override;
    [[nodiscard]] const MacGridState2D& state() const noexcept { return state_; }
    [[nodiscard]] MacGridParameters2D& parameters() noexcept { return parameters_; }
    [[nodiscard]] const MacGridParameters2D& parameters() const noexcept { return parameters_; }
    [[nodiscard]] const ProjectionResult& lastProjection() const noexcept { return lastProjection_; }
    [[nodiscard]] FluidSolverControls controls() const noexcept override;
    void setControls(const FluidSolverControls& controls) override;
    [[nodiscard]] PressureSolveDiagnostics pressureDiagnostics() const noexcept override;
    [[nodiscard]] float boxSize() const noexcept override { return state_.boxSize(); }
    void setBoxSize(float size) override;
    // Overrides are borrowed and must outlive their use by the solver.
    void setPressureSolver(IPressureSolver<MacGridState2D>& solver) noexcept;
    void useDefaultPressureSolver() noexcept { setPressureSolver(defaultPressure_); }
    void setAdvection(const IAdvection<MacAdvectionOperation2D>& method) noexcept { advection_ = &method; }
    void setAdvection(const IAdvection<MacAdvectionOperation2D>&&) = delete;
    void useDefaultAdvection() noexcept { advection_ = &defaultAdvection_; }
private:
    MacGridState2D state_;
    MacGridParameters2D parameters_;
    CpuPressureSolver2D defaultPressure_;
    SemiLagrangian2D defaultAdvection_;
    IPressureSolver<MacGridState2D>* pressure_ = &defaultPressure_;
    const IAdvection<MacAdvectionOperation2D>* advection_ = &defaultAdvection_;
    std::vector<float> uNext_, vNext_, densityNext_;
    ProjectionResult lastProjection_;
};
