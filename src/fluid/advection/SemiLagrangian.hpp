#pragma once

#include "fluid/advection/MacGridAdvection.hpp"

class SemiLagrangian final : public MacGridAdvection {
public:
    [[nodiscard]] std::string_view name() const noexcept override {
        return "Semi-Lagrangian";
    }

    [[nodiscard]] bool appliesCorrection() const noexcept override {
        return false;
    }

    void advectScalar(const AdvectionContext& context, std::span<float> field,
                      AdvectionWorkspace& workspace) const override;
    void advectFace(const AdvectionContext& context, FaceAxis axis,
                    std::span<const float> source, std::span<float> output,
                    AdvectionWorkspace& workspace) const override;
};
