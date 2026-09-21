#pragma once
#include "solvers/mac3d/AdvectionContext.hpp"

#include "common/MacAlgorithm.hpp"

class Advection3D {
public:
    virtual ~Advection3D() = default;
    virtual void advectScalar(const AdvectionContext& context, std::span<float> field,
                              AdvectionWorkspace& workspace, AdvectionMethod method) const = 0;
    virtual void advectFace(const AdvectionContext& context, FaceAxis axis,
                            std::span<const float> source, std::span<float> output,
                            AdvectionWorkspace& workspace, AdvectionMethod method) const = 0;
};
