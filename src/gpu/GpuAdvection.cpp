#include "gpu/GpuAdvection.hpp"

void GpuAdvection::advectScalar(const AdvectionContext& context, std::span<float> field,
                               AdvectionWorkspace&) const {
    validateScalarAdvection(context, field);
    backend_.validateGrid(context.grid);
    backend_.advectSmokeGPU(context.dt, field.data(), context.flow.x.data(),
        context.flow.y.data(), context.flow.z.data(), context.grid.boxSize(), context.cflLimit, correction_);
}

void GpuAdvection::advectFace(const AdvectionContext& context, FaceAxis axis,
                             std::span<const float> source, std::span<float> output,
                             AdvectionWorkspace&) const {
    validateFaceAdvection(context, axis, source, output);
    backend_.validateGrid(context.grid);
    const auto flow = context.flow;
    const int count = static_cast<int>(source.size());
    switch (axis) {
    case FaceAxis::X:
        backend_.advectVxGPU(context.dt, source.data(), count, flow.x.data(), flow.y.data(), flow.z.data(),
            output.data(), context.grid.boxSize(), context.cflLimit, correction_);
        break;
    case FaceAxis::Y:
        backend_.advectVyGPU(context.dt, source.data(), count, flow.x.data(), flow.y.data(), flow.z.data(),
            output.data(), context.grid.boxSize(), context.cflLimit, correction_);
        break;
    case FaceAxis::Z:
        backend_.advectVzGPU(context.dt, source.data(), count, flow.x.data(), flow.y.data(), flow.z.data(),
            output.data(), context.grid.boxSize(), context.cflLimit, correction_);
        break;
    }
}
