#pragma once
#include <span>

// Scalar image, bottom row first. Borrowed storage, consumed before state mutation.
struct ImageRenderData {
    std::span<const float> density;
    int width = 0, height = 0;
    float worldWidth = 1.f, worldHeight = 1.f;
    bool historicalPalette = false;
};
