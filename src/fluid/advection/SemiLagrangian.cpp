#include "fluid/advection/SemiLagrangian.hpp"
#include "fluid/advection/CpuAdvectionKernels.hpp"

void SemiLagrangian::advectScalar(const AdvectionContext& context, std::span<float> field,
                                AdvectionWorkspace& workspace) const {
    CpuAdvectionKernels::advectScalar(context, field, workspace, false);
}

void SemiLagrangian::advectFace(const AdvectionContext& context, FaceAxis axis,
                              std::span<const float> source, std::span<float> output,
                              AdvectionWorkspace& workspace) const {
    CpuAdvectionKernels::advectFace(context, axis, source, output, workspace, false);
}
