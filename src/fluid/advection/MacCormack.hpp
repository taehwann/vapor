#pragma once

#include "fluid/advection/MacGridAdvection.hpp"

// MacCormack is a correction built around a semi-Lagrangian forward/backward trace.
class MacCormack final : public MacGridAdvection {
public:
    [[nodiscard]] std::string_view name() const noexcept override {
        return "MacCormack";
    }

    [[nodiscard]] bool appliesCorrection() const noexcept override {
        return true;
    }

    void advectScalar(const AdvectionContext& context, std::span<float> field,
                      AdvectionWorkspace& workspace) const override;
    void advectFace(const AdvectionContext& context, FaceAxis axis,
                    std::span<const float> source, std::span<float> output,
                    AdvectionWorkspace& workspace) const override;
};
