#pragma once
#include "solvers/mac3d/Advection3D.hpp"

class CpuAdvection final : public Advection3D {
public:
    void advectScalar(const AdvectionContext& context, std::span<float> field,
                      AdvectionWorkspace& workspace, AdvectionMethod method) const override;
    void advectFace(const AdvectionContext& context, FaceAxis axis, std::span<const float> source,
                    std::span<float> output, AdvectionWorkspace& workspace, AdvectionMethod method) const override;
};
