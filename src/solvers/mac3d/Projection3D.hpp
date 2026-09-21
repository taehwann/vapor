#pragma once
#include "common/ProjectionOptions.hpp"

// CPU and GPU implementations of the existing 3D MAC projection.
class MacGridState3D;
class Projection3D {
public:
    virtual ~Projection3D() = default;
    virtual ProjectionResult project(MacGridState3D& state, const ProjectionOptions& options) = 0;
    virtual void reset() noexcept = 0;
};
