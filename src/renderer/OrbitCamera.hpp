#pragma once

// Build the existing orbit camera's view-projection matrix and eye position.
void buildMVP(float* m, float& cx, float& cy, float& cz, float azimuth, float elevation, float dist, float aspect, float boxSize);
