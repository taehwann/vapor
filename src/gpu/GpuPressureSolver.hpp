#pragma once

#include "fluid/projection/IPressureSolver.hpp"
#include "gpu/GpuMacBackend.hpp"

// Preserves the existing fixed-budget GPU solve; residualAvailable is false.
class GpuPressureSolver final : public IPressureSolver<MacGridState> {
public:
    explicit GpuPressureSolver(GpuMacBackend& backend) : backend_(backend) {}
    ProjectionResult project(MacGridState& state, const ProjectionOptions& options) override;
    // Device pressure is cleared at the start of every projection.
    void reset() noexcept override {}

private:
    GpuMacBackend& backend_;
};
