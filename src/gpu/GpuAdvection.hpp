#pragma once

#include "fluid/advection/MacGridAdvection.hpp"
#include "gpu/GpuMacBackend.hpp"

class GpuAdvection final : public MacGridAdvection {
public:
    explicit GpuAdvection(GpuMacBackend& backend, bool correction)
        : backend_(backend), correction_(correction) {}
    std::string_view name() const noexcept override {
        return correction_ ? "GPU MacCormack" : "GPU Semi-Lagrangian";
    }
    bool appliesCorrection() const noexcept override { return correction_; }
    void advectScalar(const AdvectionContext& context, std::span<float> field,
                      AdvectionWorkspace& workspace) const override;
    void advectFace(const AdvectionContext& context, FaceAxis axis,
                    std::span<const float> source, std::span<float> output,
                    AdvectionWorkspace& workspace) const override;

private:
    GpuMacBackend& backend_;
    bool correction_;
};
