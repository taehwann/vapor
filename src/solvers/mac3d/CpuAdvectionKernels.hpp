#pragma once

#include "solvers/mac3d/AdvectionContext.hpp"

// Shared implementation details for the two CPU strategies. No FluidSolver
// dependency: only grid/field views and reusable workspace.
namespace CpuAdvectionKernels {
void advectScalar(const AdvectionContext& context, std::span<float> field,
                  AdvectionWorkspace& workspace, bool correction);
void advectFace(const AdvectionContext& context, FaceAxis axis,
                std::span<const float> source, std::span<float> output,
                AdvectionWorkspace& workspace, bool correction);
}
