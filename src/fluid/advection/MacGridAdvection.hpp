#pragma once

#include "fluid/advection/AdvectionContext.hpp"
#include "fluid/advection/IAdvection.hpp"

#include <optional>
#include <stdexcept>

struct MacAdvectionOperation {
    AdvectionContext context;
    std::span<const float> source;
    std::span<float> output;
    AdvectionWorkspace& workspace;
    // No axis means in-place scalar transport; otherwise source/output are
    // distinct face fields. These rules belong to MAC, not generic IAdvection.
    std::optional<FaceAxis> axis;
};

// Convenience API shared by CPU/GPU implementations of MAC transport.
class MacGridAdvection : public IAdvection<MacAdvectionOperation> {
public:
    void advect(const MacAdvectionOperation& operation) const final {
        if (operation.axis) {
            advectFace(operation.context, *operation.axis, operation.source,
                       operation.output, operation.workspace);
        } else {
            if (operation.source.data() != operation.output.data() ||
                operation.source.size() != operation.output.size()) {
                throw std::invalid_argument("MAC scalar transport requires in-place source/output");
            }
            advectScalar(operation.context, operation.output, operation.workspace);
        }
    }
    [[nodiscard]] virtual bool appliesCorrection() const noexcept = 0;
    virtual void advectScalar(const AdvectionContext& context, std::span<float> field,
                              AdvectionWorkspace& workspace) const = 0;
    virtual void advectFace(const AdvectionContext& context, FaceAxis axis,
                            std::span<const float> source, std::span<float> output,
                            AdvectionWorkspace& workspace) const = 0;
};
