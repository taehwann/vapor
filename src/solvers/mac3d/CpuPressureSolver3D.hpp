#pragma once

#include "solvers/mac3d/Projection3D.hpp"
#include "solvers/mac3d/SOR.hpp"

#include <span>
#include <vector>

class MacGridState3D;

// Closed-box MAC pressure projection with an owned SOR workspace.
class CpuPressureSolver3D final : public Projection3D {
public:
    // Relaxation is selected per call in ProjectionOptions.
    explicit CpuPressureSolver3D(float diagonalFloor = 0.f)
        : sor_(1.95f, diagonalFloor) {}
    CpuPressureSolver3D(const CpuPressureSolver3D&) = delete;
    CpuPressureSolver3D& operator=(const CpuPressureSolver3D&) = delete;
    ProjectionResult project(MacGridState3D& state, const ProjectionOptions& options) override;
    [[nodiscard]] std::span<const float> pressure() const noexcept { return system_.x; }
    [[nodiscard]] std::span<const float> inputDivergence() const noexcept { return divergence_; }
    void reset() noexcept override;

    [[nodiscard]] const StencilLinearSystem& system() const noexcept { return system_; }
private:
    void prepareSystem(int resolution);
    SOR sor_;
    StencilLinearSystem system_;
    std::vector<float> divergence_;
};
