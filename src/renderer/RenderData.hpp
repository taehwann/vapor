#pragma once

#include <span>

struct RenderData {
    std::span<const float> density;
    int width = 0;
    int height = 0;
    int depth = 0;
    float boxSize = 1.0f;
    // Optional borrowed OpenGL buffer with width*height*depth floats.
    // When nonzero, density can be empty; no CPU transfer is needed to render.
    unsigned int densityBuffer = 0;
};
