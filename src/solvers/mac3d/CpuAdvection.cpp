#include "solvers/mac3d/CpuAdvection.hpp"
#include "solvers/mac3d/CpuAdvectionKernels.hpp"

void CpuAdvection::advectScalar(const AdvectionContext& context, std::span<float> field,
                               AdvectionWorkspace& workspace, AdvectionMethod method) const {
    CpuAdvectionKernels::advectScalar(context, field, workspace, method == AdvectionMethod::MacCormack);
}
void CpuAdvection::advectFace(const AdvectionContext& context, FaceAxis axis, std::span<const float> source,
                             std::span<float> output, AdvectionWorkspace& workspace, AdvectionMethod method) const {
    CpuAdvectionKernels::advectFace(context, axis, source, output, workspace, method == AdvectionMethod::MacCormack);
}
