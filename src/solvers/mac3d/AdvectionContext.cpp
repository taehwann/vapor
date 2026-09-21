#include "solvers/mac3d/AdvectionContext.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <stdexcept>

namespace {
bool overlaps(std::span<const float> a, std::span<const float> b) {
    if (a.empty() || b.empty()) return false;
    const std::less<const float*> less;
    return less(a.data(), b.data() + b.size()) && less(b.data(), a.data() + a.size());
}

void validate(const AdvectionContext& context) {
    const auto& grid = context.grid;
    if (!std::isfinite(context.dt) || context.dt < 0.f ||
        !std::isfinite(context.cflLimit) || context.cflLimit <= 0.f ||
        context.flow.x.size() != grid.velocityX().size() ||
        context.flow.y.size() != grid.velocityY().size() ||
        context.flow.z.size() != grid.velocityZ().size()) {
        throw std::invalid_argument("Invalid advection time step, CFL limit, or tracer field sizes");
    }
}

void rejectTracerOverlap(std::span<const float> output, VelocityFieldView flow) {
    if (overlaps(output, flow.x) || overlaps(output, flow.y) || overlaps(output, flow.z)) {
        throw std::invalid_argument("Advection output must not alias tracer velocities");
    }
}

}

void AdvectionWorkspace::resize(const MacGridState3D& grid) {
    const auto count = std::max({grid.velocityX().size(), grid.velocityY().size(), grid.velocityZ().size()});
    forward.resize(count);
    backward.resize(count);
}

void AdvectionWorkspace::reset() {
    std::fill(forward.begin(), forward.end(), 0.f);
    std::fill(backward.begin(), backward.end(), 0.f);
}

void validateScalarAdvection(const AdvectionContext& context, std::span<float> field) {
    validate(context);
    if (field.size() != context.grid.cellCount()) {
        throw std::invalid_argument("Scalar field must match the grid cell count");
    }
    rejectTracerOverlap(field, context.flow);
}

void validateFaceAdvection(const AdvectionContext& context, FaceAxis axis,
                           std::span<const float> source, std::span<float> output) {
    validate(context);
    std::size_t count;
    switch (axis) {
        case FaceAxis::X: count = context.grid.velocityX().size(); break;
        case FaceAxis::Y: count = context.grid.velocityY().size(); break;
        case FaceAxis::Z: count = context.grid.velocityZ().size(); break;
        default: throw std::invalid_argument("Invalid staggered field axis");
    }
    if (source.size() != count || output.size() != count || overlaps(source, output)) {
        throw std::invalid_argument("Velocity fields must match their grid layout and not overlap");
    }
    rejectTracerOverlap(output, context.flow);
}
