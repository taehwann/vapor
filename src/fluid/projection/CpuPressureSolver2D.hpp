#pragma once
#include "fluid/projection/IPressureSolver.hpp"
#include "fluid/MacGridState2D.hpp"
#include <vector>

// Five-point, red-black SOR with homogeneous Neumann pressure walls.
// Owns only pressure workspace, never the simulation fields.
class CpuPressureSolver2D final : public IPressureSolver<MacGridState2D> {
public:
    ProjectionResult project(MacGridState2D& state, const ProjectionOptions& options) override;
    void reset() noexcept override;
    [[nodiscard]] std::span<const double> pressure() const noexcept { return pressure_; }
private:
    // Double pressure workspace avoids a float residual floor under sustained buoyancy.
    // Transported physical fields remain float.
    std::vector<double> pressure_;
    std::vector<float> rhs_;
};
