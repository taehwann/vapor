#pragma once
#include <memory>
#include <stdexcept>

#include "solvers/mac3d/Advection3D.hpp"
#include "solvers/mac3d/gpu/GpuMacBackend.hpp"

class GpuAdvection final : public Advection3D {
public:
    explicit GpuAdvection(std::shared_ptr<GpuMacBackend> backend) : backend_(std::move(backend)) {
        if (!backend_) throw std::invalid_argument("GPU backend cannot be null");
    }
    void advectScalar(const AdvectionContext& context, std::span<float> field,
                      AdvectionWorkspace& workspace, AdvectionMethod method) const override;
    void advectFace(const AdvectionContext& context, FaceAxis axis,
                    std::span<const float> source, std::span<float> output,
                    AdvectionWorkspace& workspace, AdvectionMethod method) const override;

private:
    std::shared_ptr<GpuMacBackend> backend_;
};
