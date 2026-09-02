#pragma once

#include "fluid/projection/IPressureSolver.hpp"
#include "numerics/SOR.hpp"

#include <span>
#include <vector>

class MacGridState;

// Closed-box MAC pressure projection. The default SOR is owned; an optionally
// injected numerical solver is borrowed and must outlive this object.
class CpuPressureSolver final : public IPressureSolver<MacGridState> {
public:
    // Relaxation is selected per call in ProjectionOptions. An injected linear
    // solver instead retains its own solver-specific configuration.
    explicit CpuPressureSolver(float diagonalFloor = 0.f)
        : defaultSolver_(1.95f, diagonalFloor), linearSolver_(&defaultSolver_) {}
    explicit CpuPressureSolver(const ILinearSolver& solver) : linearSolver_(&solver) {}
    CpuPressureSolver(const ILinearSolver&&) = delete;
    CpuPressureSolver(const CpuPressureSolver&) = delete;
    CpuPressureSolver& operator=(const CpuPressureSolver&) = delete;
    ProjectionResult project(MacGridState& state, const ProjectionOptions& options) override;
    [[nodiscard]] std::span<const float> pressure() const noexcept { return system_.x; }
    [[nodiscard]] std::span<const float> inputDivergence() const noexcept { return divergence_; }
    void reset() noexcept override;

private:
    void prepareSystem(int resolution);
    SOR defaultSolver_;
    const ILinearSolver* linearSolver_;
    StencilLinearSystem system_;
    std::vector<float> divergence_;
};
