#pragma once
#include "common/ProjectionOptions.hpp"
#include "solvers/mac2d/MacGridState2D.hpp"
#include <vector>

// Five-point, red-black SOR with homogeneous Neumann pressure walls.
// Owns only pressure workspace, never the simulation fields.
class CpuPressureSolver2D final {
public:
    ProjectionResult project(MacGridState2D& state, const ProjectionOptions& options);
    void reset() noexcept;
    [[nodiscard]] std::span<const double> pressure() const noexcept { return pressure_; }
private:
    // Double pressure workspace avoids a float residual floor under sustained buoyancy.
    // Transported physical fields remain float.
    std::vector<double> pressure_;
    std::vector<float> rhs_;
};
