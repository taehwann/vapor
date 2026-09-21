#pragma once
#include <memory>
#include <stdexcept>

#include "solvers/mac3d/Projection3D.hpp"
#include "solvers/mac3d/gpu/GpuMacBackend.hpp"

// Preserves the existing fixed-budget GPU solve; residualAvailable is false.
class GpuPressureSolver final : public Projection3D {
public:
    explicit GpuPressureSolver(std::shared_ptr<GpuMacBackend> backend) : backend_(std::move(backend)) {
        if (!backend_) throw std::invalid_argument("GPU backend cannot be null");
    }
    ProjectionResult project(MacGridState3D& state, const ProjectionOptions& options) override;
    // Device pressure is cleared at the start of every projection.
    void reset() noexcept override {}

private:
    std::shared_ptr<GpuMacBackend> backend_;
};
