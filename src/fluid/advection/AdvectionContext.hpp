#pragma once

#include "fluid/VelocityFieldView.hpp"

#include <span>
#include <vector>

enum class FaceAxis { X, Y, Z };

// Borrowed inputs for one advection operation; no UI or backend settings.
struct AdvectionContext {
    const MacGridState& grid;
    VelocityFieldView flow;
    float dt;
    float cflLimit = 3.f;
};

// Reusable algorithm workspace, separate from physical state. Do not alias these
// buffers with inputs or outputs. One workspace per concurrent operation.
struct AdvectionWorkspace {
    std::vector<float> forward, backward;
    void resize(const MacGridState& grid);
    void reset();
};

// The same field-size and aliasing contract is enforced by both backends.
void validateScalarAdvection(const AdvectionContext& context, std::span<float> field);
void validateFaceAdvection(const AdvectionContext& context, FaceAxis axis,
                           std::span<const float> source, std::span<float> output);
