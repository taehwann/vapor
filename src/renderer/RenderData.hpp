#pragma once

#include <span>

struct RenderData {
    std::span<const float> density;
    int width = 0;
    int height = 0;
    int depth = 0;
    float boxSize = 1.0f;
};
