#include "renderer/OrbitCamera.hpp"

#include <cmath>

void buildMVP(float* m, float& cx, float& cy, float& cz, float azimuth, float elevation, float dist, float aspect, float boxSize) {
    float cax = cosf(azimuth), sax = sinf(azimuth), cel = cosf(elevation), sel = sinf(elevation);
    float half = boxSize * 0.5f;
    cx = half + dist * cel * cax;
    cy = half + dist * cel * sax;
    cz = half + dist * sel;
    float upx = -sel * cax, upy = -sel * sax, upz = cel;
    float fx = half - cx, fy = half - cy, fz = half - cz, fl = sqrtf(fx * fx + fy * fy + fz * fz);
    fx /= fl; fy /= fl; fz /= fl;
    float rx = fy * upz - fz * upy, ry = fz * upx - fx * upz, rz = fx * upy - fy * upx;
    float proj[16] = {};
    float f = 1.0f / tanf(0.8f * 0.5f);
    proj[0] = f / aspect; proj[5] = f; proj[10] = -101.f / 99.f; proj[11] = -1.f; proj[14] = -202.f / 99.f;
    float view[16] = {};
    view[0] = rx; view[4] = ry; view[8] = rz; view[12] = -(rx * cx + ry * cy + rz * cz);
    view[1] = upx; view[5] = upy; view[9] = upz; view[13] = -(upx * cx + upy * cy + upz * cz);
    view[2] = -fx; view[6] = -fy; view[10] = -fz; view[14] = fx * cx + fy * cy + fz * cz;
    view[3] = 0; view[7] = 0; view[11] = 0; view[15] = 1.f;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j) {
            m[i + 4 * j] = 0;
            for (int k = 0; k < 4; ++k) m[i + 4 * j] += proj[i + 4 * k] * view[k + 4 * j];
        }
}
