#include <algorithm>
#include <cmath>
#include <cstdarg>
#include <cstdint>
#include <cstdio>
#include <vector>

#include "gl_loader.h"
#include <GLFW/glfw3.h>
#include <imgui.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_opengl3.h>

struct Vec3 { float x, y, z; };
static Vec3 operator-(Vec3 a, Vec3 b) { return { a.x - b.x, a.y - b.y, a.z - b.z }; }
static Vec3 operator+(Vec3 a, Vec3 b) { return { a.x + b.x, a.y + b.y, a.z + b.z }; }
static Vec3 operator*(Vec3 a, float s) { return { a.x * s, a.y * s, a.z * s }; }

static float frand() {
    static uint32_t state = 0x9E3779B9u;
    state = state * 1664525u + 1013904223u;
    return (state >> 8) * (1.0f / 16777216.0f);
}

struct SmokeSim3D {
    int n = 64;
    float boxSize = 3.0f;
    float emitterRadius = 0.05f;
    float emitterCenterX = 0.5f, emitterCenterY = 0.5f, emitterCenterZ = 0.15f;
    float emitterKickBase = 1.2f;
    float emitterKickScale = 1.0f;
    float stirStrength = 0.6f;
    float sorOmega = 1.95f;
    float buoyancySplit = 0.5f;
    float cflMc = 3.0f;
    bool macCormackSmoke = true, macCormackVel = true;
    bool openXm = true, openXp = true, openYm = true, openYp = true, openZm = false, openZp = true;
    std::vector<float> smoke, vx, vy, vz, pressure, divergence;

    explicit SmokeSim3D(int res) {
        n = res;
        smoke.assign(n * n * n, 0);
        vx.assign((n + 1) * n * n, 0);
        vy.assign(n * (n + 1) * n, 0);
        vz.assign(n * n * (n + 1), 0);
        pressure.assign(n * n * n, 0);
        divergence.assign(n * n * n, 0);
    }

    int idC(int x, int y, int z) const { return x + n * (y + n * z); }
    int idX(int x, int y, int z) const { return x + (n + 1) * (y + n * z); }
    int idY(int x, int y, int z) const { return x + n * (y + (n + 1) * z); }
    int idZ(int x, int y, int z) const { return x + n * (y + n * z); }
    float h() const { return boxSize / n; }
    bool inside(int x, int y, int z) const { return x >= 0 && y >= 0 && z >= 0 && x < n && y < n && z < n; }
    bool insideVx(int x, int y, int z) const { return x >= 0 && y >= 0 && z >= 0 && x <= n && y < n && z < n; }
    bool insideVy(int x, int y, int z) const { return x >= 0 && y >= 0 && z >= 0 && x < n && y <= n && z < n; }
    bool insideVz(int x, int y, int z) const { return x >= 0 && y >= 0 && z >= 0 && x < n && y < n && z <= n; }

    float sampleCell(const std::vector<float>& f, Vec3 p) const {
        const float gx = std::clamp(p.x / h() - 0.5f, 0.f, float(n - 1));
        const float gy = std::clamp(p.y / h() - 0.5f, 0.f, float(n - 1));
        const float gz = std::clamp(p.z / h() - 0.5f, 0.f, float(n - 1));
        const int x0 = int(gx), y0 = int(gy), z0 = int(gz);
        const int x1 = std::min(x0 + 1, n - 1), y1 = std::min(y0 + 1, n - 1), z1 = std::min(z0 + 1, n - 1);
        const float tx = gx - x0, ty = gy - y0, tz = gz - z0;
        return f[idC(x0, y0, z0)] * (1 - tx) * (1 - ty) * (1 - tz) +
            f[idC(x1, y0, z0)] * tx * (1 - ty) * (1 - tz) +
            f[idC(x0, y1, z0)] * (1 - tx) * ty * (1 - tz) +
            f[idC(x1, y1, z0)] * tx * ty * (1 - tz) +
            f[idC(x0, y0, z1)] * (1 - tx) * (1 - ty) * tz +
            f[idC(x1, y0, z1)] * tx * (1 - ty) * tz +
            f[idC(x0, y1, z1)] * (1 - tx) * ty * tz +
            f[idC(x1, y1, z1)] * tx * ty * tz;
    }

    float sampleVxField(const std::vector<float>& f, Vec3 p) const {
        const float gx = std::clamp(p.x / h(), 0.f, float(n));
        const float gy = std::clamp(p.y / h() - 0.5f, 0.f, float(n - 1));
        const float gz = std::clamp(p.z / h() - 0.5f, 0.f, float(n - 1));
        const int x0 = std::min(int(gx), n - 1), y0 = int(gy), z0 = int(gz);
        const int x1 = x0 + 1, y1 = std::min(y0 + 1, n - 1), z1 = std::min(z0 + 1, n - 1);
        const float tx = gx - x0, ty = gy - y0, tz = gz - z0;
        return f[idX(x0, y0, z0)] * (1 - tx) * (1 - ty) * (1 - tz) +
            f[idX(x1, y0, z0)] * tx * (1 - ty) * (1 - tz) +
            f[idX(x0, y1, z0)] * (1 - tx) * ty * (1 - tz) +
            f[idX(x1, y1, z0)] * tx * ty * (1 - tz) +
            f[idX(x0, y0, z1)] * (1 - tx) * (1 - ty) * tz +
            f[idX(x1, y0, z1)] * tx * (1 - ty) * tz +
            f[idX(x0, y1, z1)] * (1 - tx) * ty * tz +
            f[idX(x1, y1, z1)] * tx * ty * tz;
    }

    float sampleVyField(const std::vector<float>& f, Vec3 p) const {
        const float gx = std::clamp(p.x / h() - 0.5f, 0.f, float(n - 1));
        const float gy = std::clamp(p.y / h(), 0.f, float(n));
        const float gz = std::clamp(p.z / h() - 0.5f, 0.f, float(n - 1));
        const int x0 = int(gx), y0 = std::min(int(gy), n - 1), z0 = int(gz);
        const int x1 = std::min(x0 + 1, n - 1), y1 = y0 + 1, z1 = std::min(z0 + 1, n - 1);
        const float tx = gx - x0, ty = gy - y0, tz = gz - z0;
        return f[idY(x0, y0, z0)] * (1 - tx) * (1 - ty) * (1 - tz) +
            f[idY(x1, y0, z0)] * tx * (1 - ty) * (1 - tz) +
            f[idY(x0, y1, z0)] * (1 - tx) * ty * (1 - tz) +
            f[idY(x1, y1, z0)] * tx * ty * (1 - tz) +
            f[idY(x0, y0, z1)] * (1 - tx) * (1 - ty) * tz +
            f[idY(x1, y0, z1)] * tx * (1 - ty) * tz +
            f[idY(x0, y1, z1)] * (1 - tx) * ty * tz +
            f[idY(x1, y1, z1)] * tx * ty * tz;
    }

    float sampleVzField(const std::vector<float>& f, Vec3 p) const {
        const float gx = std::clamp(p.x / h() - 0.5f, 0.f, float(n - 1));
        const float gy = std::clamp(p.y / h() - 0.5f, 0.f, float(n - 1));
        const float gz = std::clamp(p.z / h(), 0.f, float(n));
        const int x0 = int(gx), y0 = int(gy), z0 = std::min(int(gz), n - 1);
        const int x1 = std::min(x0 + 1, n - 1), y1 = std::min(y0 + 1, n - 1), z1 = z0 + 1;
        const float tx = gx - x0, ty = gy - y0, tz = gz - z0;
        return f[idZ(x0, y0, z0)] * (1 - tx) * (1 - ty) * (1 - tz) +
            f[idZ(x1, y0, z0)] * tx * (1 - ty) * (1 - tz) +
            f[idZ(x0, y1, z0)] * (1 - tx) * ty * (1 - tz) +
            f[idZ(x1, y1, z0)] * tx * ty * (1 - tz) +
            f[idZ(x0, y0, z1)] * (1 - tx) * (1 - ty) * tz +
            f[idZ(x1, y0, z1)] * tx * (1 - ty) * tz +
            f[idZ(x0, y1, z1)] * (1 - tx) * ty * tz +
            f[idZ(x1, y1, z1)] * tx * ty * tz;
    }

    float sampleVx(Vec3 p) const { return sampleVxField(vx, p); }
    float sampleVy(Vec3 p) const { return sampleVyField(vy, p); }
    float sampleVz(Vec3 p) const { return sampleVzField(vz, p); }
    Vec3 velocity(Vec3 p) const { return { sampleVx(p), sampleVy(p), sampleVz(p) }; }
    Vec3 velocity(Vec3 p, const std::vector<float>& vxSrc, const std::vector<float>& vySrc, const std::vector<float>& vzSrc) const {
        return { sampleVxField(vxSrc, p), sampleVyField(vySrc, p), sampleVzField(vzSrc, p) };
    }

    float sampleCellDx(const std::vector<float>& f, Vec3 p, float dx) const {
        const float gx = std::clamp(p.x / dx - 0.5f, 0.f, float(n - 1));
        const float gy = std::clamp(p.y / dx - 0.5f, 0.f, float(n - 1));
        const float gz = std::clamp(p.z / dx - 0.5f, 0.f, float(n - 1));
        const int x0 = int(gx), y0 = int(gy), z0 = int(gz);
        const int x1 = std::min(x0 + 1, n - 1), y1 = std::min(y0 + 1, n - 1), z1 = std::min(z0 + 1, n - 1);
        const float tx = gx - x0, ty = gy - y0, tz = gz - z0;
        return f[idC(x0, y0, z0)] * (1 - tx) * (1 - ty) * (1 - tz) +
            f[idC(x1, y0, z0)] * tx * (1 - ty) * (1 - tz) +
            f[idC(x0, y1, z0)] * (1 - tx) * ty * (1 - tz) +
            f[idC(x1, y1, z0)] * tx * ty * (1 - tz) +
            f[idC(x0, y0, z1)] * (1 - tx) * (1 - ty) * tz +
            f[idC(x1, y0, z1)] * tx * (1 - ty) * tz +
            f[idC(x0, y1, z1)] * (1 - tx) * ty * tz +
            f[idC(x1, y1, z1)] * tx * ty * tz;
    }

    float sampleVxFieldDx(const std::vector<float>& f, Vec3 p, float dx) const {
        const float gx = std::clamp(p.x / dx, 0.f, float(n));
        const float gy = std::clamp(p.y / dx - 0.5f, 0.f, float(n - 1));
        const float gz = std::clamp(p.z / dx - 0.5f, 0.f, float(n - 1));
        const int x0 = std::min(int(gx), n - 1), y0 = int(gy), z0 = int(gz);
        const int x1 = x0 + 1, y1 = std::min(y0 + 1, n - 1), z1 = std::min(z0 + 1, n - 1);
        const float tx = gx - x0, ty = gy - y0, tz = gz - z0;
        return f[idX(x0, y0, z0)] * (1 - tx) * (1 - ty) * (1 - tz) +
            f[idX(x1, y0, z0)] * tx * (1 - ty) * (1 - tz) +
            f[idX(x0, y1, z0)] * (1 - tx) * ty * (1 - tz) +
            f[idX(x1, y1, z0)] * tx * ty * (1 - tz) +
            f[idX(x0, y0, z1)] * (1 - tx) * (1 - ty) * tz +
            f[idX(x1, y0, z1)] * tx * (1 - ty) * tz +
            f[idX(x0, y1, z1)] * (1 - tx) * ty * tz +
            f[idX(x1, y1, z1)] * tx * ty * tz;
    }

    float sampleVyFieldDx(const std::vector<float>& f, Vec3 p, float dx) const {
        const float gx = std::clamp(p.x / dx - 0.5f, 0.f, float(n - 1));
        const float gy = std::clamp(p.y / dx, 0.f, float(n));
        const float gz = std::clamp(p.z / dx - 0.5f, 0.f, float(n - 1));
        const int x0 = int(gx), y0 = std::min(int(gy), n - 1), z0 = int(gz);
        const int x1 = std::min(x0 + 1, n - 1), y1 = y0 + 1, z1 = std::min(z0 + 1, n - 1);
        const float tx = gx - x0, ty = gy - y0, tz = gz - z0;
        return f[idY(x0, y0, z0)] * (1 - tx) * (1 - ty) * (1 - tz) +
            f[idY(x1, y0, z0)] * tx * (1 - ty) * (1 - tz) +
            f[idY(x0, y1, z0)] * (1 - tx) * ty * (1 - tz) +
            f[idY(x1, y1, z0)] * tx * ty * (1 - tz) +
            f[idY(x0, y0, z1)] * (1 - tx) * (1 - ty) * tz +
            f[idY(x1, y0, z1)] * tx * (1 - ty) * tz +
            f[idY(x0, y1, z1)] * (1 - tx) * ty * tz +
            f[idY(x1, y1, z1)] * tx * ty * tz;
    }

    float sampleVzFieldDx(const std::vector<float>& f, Vec3 p, float dx) const {
        const float gx = std::clamp(p.x / dx - 0.5f, 0.f, float(n - 1));
        const float gy = std::clamp(p.y / dx - 0.5f, 0.f, float(n - 1));
        const float gz = std::clamp(p.z / dx, 0.f, float(n));
        const int x0 = int(gx), y0 = int(gy), z0 = std::min(int(gz), n - 1);
        const int x1 = std::min(x0 + 1, n - 1), y1 = std::min(y0 + 1, n - 1), z1 = z0 + 1;
        const float tx = gx - x0, ty = gy - y0, tz = gz - z0;
        return f[idZ(x0, y0, z0)] * (1 - tx) * (1 - ty) * (1 - tz) +
            f[idZ(x1, y0, z0)] * tx * (1 - ty) * (1 - tz) +
            f[idZ(x0, y1, z0)] * (1 - tx) * ty * (1 - tz) +
            f[idZ(x1, y1, z0)] * tx * ty * (1 - tz) +
            f[idZ(x0, y0, z1)] * (1 - tx) * (1 - ty) * tz +
            f[idZ(x1, y0, z1)] * tx * (1 - ty) * tz +
            f[idZ(x0, y1, z1)] * (1 - tx) * ty * tz +
            f[idZ(x1, y1, z1)] * tx * ty * tz;
    }

    void advectScalar(std::vector<float>& field, float dt) {
        const float dx = h();
        const int N = n * n * n;
        std::vector<float> fwd(N);
        for (int z = 0; z < n; ++z)
            for (int y = 0; y < n; ++y)
                for (int x = 0; x < n; ++x)
                    fwd[idC(x, y, z)] = sampleCell(field, Vec3{ (x + .5f) * dx, (y + .5f) * dx, (z + .5f) * dx } - velocity(Vec3{ (x + .5f) * dx, (y + .5f) * dx, (z + .5f) * dx }) * dt);

        if (!macCormackSmoke) { field.swap(fwd); return; }

        std::vector<float> back(N);
        for (int z = 0; z < n; ++z)
            for (int y = 0; y < n; ++y)
                for (int x = 0; x < n; ++x) {
                    Vec3 p{ (x + .5f) * dx, (y + .5f) * dx, (z + .5f) * dx };
                    back[idC(x, y, z)] = sampleCell(fwd, p + velocity(p) * dt);
                }

        for (int z = 0; z < n; ++z)
            for (int y = 0; y < n; ++y)
                for (int x = 0; x < n; ++x) {
                    const int i = idC(x, y, z);
                    Vec3 p{ (x + .5f) * dx, (y + .5f) * dx, (z + .5f) * dx };
                    Vec3 posFwd = p - velocity(p) * dt;
                    Vec3 posBwd = p + velocity(p) * dt;
                    if (posFwd.x < 0 || posFwd.x > boxSize || posFwd.y < 0 || posFwd.y > boxSize || posFwd.z < 0 || posFwd.z > boxSize ||
                        posBwd.x < 0 || posBwd.x > boxSize || posBwd.y < 0 || posBwd.y > boxSize || posBwd.z < 0 || posBwd.z > boxSize) {
                        field[i] = fwd[i]; continue;
                    }
                    float corrected = fwd[i] + 0.5f * (field[i] - back[i]);
                    int gx = std::clamp(int(posFwd.x / dx - 0.5f), 0, n - 1);
                    int gy = std::clamp(int(posFwd.y / dx - 0.5f), 0, n - 1);
                    int gz = std::clamp(int(posFwd.z / dx - 0.5f), 0, n - 1);
                    float fMin = field[idC(gx, gy, gz)], fMax = fMin;
                    for (int dz = 0; dz <= 1; ++dz)
                        for (int dy = 0; dy <= 1; ++dy)
                            for (int dx2 = 0; dx2 <= 1; ++dx2) {
                                int a = std::clamp(gx + dx2, 0, n - 1);
                                int b = std::clamp(gy + dy, 0, n - 1);
                                int c = std::clamp(gz + dz, 0, n - 1);
                                float val = field[idC(a, b, c)];
                                fMin = std::min(fMin, val);
                                fMax = std::max(fMax, val);
                            }
                    if (corrected < fMin || corrected > fMax) field[i] = fwd[i];
                    else field[i] = corrected;
                }
    }

    void advectVx(float dt, const std::vector<float>& vxSrc,
        const std::vector<float>& vxVel, const std::vector<float>& vyVel, const std::vector<float>& vzVel,
        std::vector<float>& out, bool useMC = true) {
        const float dx = h();
        const int N = (int)out.size();
        std::vector<float> fwd(N);
        for (int z = 0; z < n; ++z)
            for (int y = 0; y < n; ++y)
                for (int x = 0; x <= n; ++x) {
                    Vec3 face{ x * dx, (y + .5f) * dx, (z + .5f) * dx };
                    fwd[idX(x, y, z)] = sampleVxField(vxSrc, face - velocity(face, vxVel, vyVel, vzVel) * dt);
                }

        if (!useMC || !macCormackVel) { out.swap(fwd); return; }

        std::vector<float> back(N);
        for (int z = 0; z < n; ++z)
            for (int y = 0; y < n; ++y)
                for (int x = 0; x <= n; ++x) {
                    Vec3 face{ x * dx, (y + .5f) * dx, (z + .5f) * dx };
                    back[idX(x, y, z)] = sampleVxField(fwd, face + velocity(face, vxVel, vyVel, vzVel) * dt);
                }

        for (int z = 0; z < n; ++z)
            for (int y = 0; y < n; ++y)
                for (int x = 0; x <= n; ++x) {
                    const int i = idX(x, y, z);
                    Vec3 face{ x * dx, (y + .5f) * dx, (z + .5f) * dx };
                    Vec3 v = velocity(face, vxVel, vyVel, vzVel);
                    Vec3 posFwd = face - v * dt, posBwd = face + v * dt;
                    float cfl = std::sqrt(v.x * v.x + v.y * v.y + v.z * v.z) * dt / dx;
                    if (posFwd.x < 0 || posFwd.x > boxSize || posFwd.y < 0 || posFwd.y > boxSize || posFwd.z < 0 || posFwd.z > boxSize ||
                        posBwd.x < 0 || posBwd.x > boxSize || posBwd.y < 0 || posBwd.y > boxSize || posBwd.z < 0 || posBwd.z > boxSize ||
                        cfl > cflMc) {
                        out[i] = fwd[i]; continue;
                    }
                    float corrected = fwd[i] + 0.5f * (vxSrc[i] - back[i]);
                    int gx = std::clamp(int(posFwd.x / dx), 0, n);
                    int gy = std::clamp(int(posFwd.y / dx - 0.5f), 0, n - 1);
                    int gz = std::clamp(int(posFwd.z / dx - 0.5f), 0, n - 1);
                    float fMin = vxSrc[idX(gx, gy, gz)], fMax = fMin;
                    for (int dz = 0; dz <= 1; ++dz)
                        for (int dy = 0; dy <= 1; ++dy)
                            for (int dx2 = 0; dx2 <= 1; ++dx2) {
                                int a = std::clamp(gx + dx2, 0, n);
                                int b = std::clamp(gy + dy, 0, n - 1);
                                int c = std::clamp(gz + dz, 0, n - 1);
                                float val = vxSrc[idX(a, b, c)];
                                fMin = std::min(fMin, val);
                                fMax = std::max(fMax, val);
                            }
                    if (corrected < fMin || corrected > fMax) out[i] = fwd[i];
                    else out[i] = corrected;
                }
    }

    void advectVy(float dt, const std::vector<float>& vySrc,
        const std::vector<float>& vxVel, const std::vector<float>& vyVel, const std::vector<float>& vzVel,
        std::vector<float>& out, bool useMC = true) {
        const float dx = h();
        const int N = (int)out.size();
        std::vector<float> fwd(N);
        for (int z = 0; z < n; ++z)
            for (int y = 0; y <= n; ++y)
                for (int x = 0; x < n; ++x) {
                    Vec3 face{ (x + .5f) * dx, y * dx, (z + .5f) * dx };
                    fwd[idY(x, y, z)] = sampleVyField(vySrc, face - velocity(face, vxVel, vyVel, vzVel) * dt);
                }

        if (!useMC || !macCormackVel) { out.swap(fwd); return; }

        std::vector<float> back(N);
        for (int z = 0; z < n; ++z)
            for (int y = 0; y <= n; ++y)
                for (int x = 0; x < n; ++x) {
                    Vec3 face{ (x + .5f) * dx, y * dx, (z + .5f) * dx };
                    back[idY(x, y, z)] = sampleVyField(fwd, face + velocity(face, vxVel, vyVel, vzVel) * dt);
                }

        for (int z = 0; z < n; ++z)
            for (int y = 0; y <= n; ++y)
                for (int x = 0; x < n; ++x) {
                    const int i = idY(x, y, z);
                    Vec3 face{ (x + .5f) * dx, y * dx, (z + .5f) * dx };
                    Vec3 v = velocity(face, vxVel, vyVel, vzVel);
                    Vec3 posFwd = face - v * dt, posBwd = face + v * dt;
                    float cfl = std::sqrt(v.x * v.x + v.y * v.y + v.z * v.z) * dt / dx;
                    if (posFwd.x < 0 || posFwd.x > boxSize || posFwd.y < 0 || posFwd.y > boxSize || posFwd.z < 0 || posFwd.z > boxSize ||
                        posBwd.x < 0 || posBwd.x > boxSize || posBwd.y < 0 || posBwd.y > boxSize || posBwd.z < 0 || posBwd.z > boxSize ||
                        cfl > cflMc) {
                        out[i] = fwd[i]; continue;
                    }
                    float corrected = fwd[i] + 0.5f * (vySrc[i] - back[i]);
                    int gx = std::clamp(int(posFwd.x / dx - 0.5f), 0, n - 1);
                    int gy = std::clamp(int(posFwd.y / dx), 0, n);
                    int gz = std::clamp(int(posFwd.z / dx - 0.5f), 0, n - 1);
                    float fMin = vySrc[idY(gx, gy, gz)], fMax = fMin;
                    for (int dz = 0; dz <= 1; ++dz)
                        for (int dy = 0; dy <= 1; ++dy)
                            for (int dx2 = 0; dx2 <= 1; ++dx2) {
                                int a = std::clamp(gx + dx2, 0, n - 1);
                                int b = std::clamp(gy + dy, 0, n);
                                int c = std::clamp(gz + dz, 0, n - 1);
                                float val = vySrc[idY(a, b, c)];
                                fMin = std::min(fMin, val);
                                fMax = std::max(fMax, val);
                            }
                    if (corrected < fMin || corrected > fMax) out[i] = fwd[i];
                    else out[i] = corrected;
                }
    }

    void advectVz(float dt, const std::vector<float>& vzSrc,
        const std::vector<float>& vxVel, const std::vector<float>& vyVel, const std::vector<float>& vzVel,
        std::vector<float>& out, bool useMC = true) {
        const float dx = h();
        const int N = (int)out.size();
        std::vector<float> fwd(N);
        for (int z = 0; z <= n; ++z)
            for (int y = 0; y < n; ++y)
                for (int x = 0; x < n; ++x) {
                    Vec3 face{ (x + .5f) * dx, (y + .5f) * dx, z * dx };
                    fwd[idZ(x, y, z)] = sampleVzField(vzSrc, face - velocity(face, vxVel, vyVel, vzVel) * dt);
                }

        if (!useMC || !macCormackVel) { out.swap(fwd); return; }

        std::vector<float> back(N);
        for (int z = 0; z <= n; ++z)
            for (int y = 0; y < n; ++y)
                for (int x = 0; x < n; ++x) {
                    Vec3 face{ (x + .5f) * dx, (y + .5f) * dx, z * dx };
                    back[idZ(x, y, z)] = sampleVzField(fwd, face + velocity(face, vxVel, vyVel, vzVel) * dt);
                }

        for (int z = 0; z <= n; ++z)
            for (int y = 0; y < n; ++y)
                for (int x = 0; x < n; ++x) {
                    const int i = idZ(x, y, z);
                    Vec3 face{ (x + .5f) * dx, (y + .5f) * dx, z * dx };
                    Vec3 v = velocity(face, vxVel, vyVel, vzVel);
                    Vec3 posFwd = face - v * dt, posBwd = face + v * dt;
                    float cfl = std::sqrt(v.x * v.x + v.y * v.y + v.z * v.z) * dt / dx;
                    if (posFwd.x < 0 || posFwd.x > boxSize || posFwd.y < 0 || posFwd.y > boxSize || posFwd.z < 0 || posFwd.z > boxSize ||
                        posBwd.x < 0 || posBwd.x > boxSize || posBwd.y < 0 || posBwd.y > boxSize || posBwd.z < 0 || posBwd.z > boxSize ||
                        cfl > cflMc) {
                        out[i] = fwd[i]; continue;
                    }
                    float corrected = fwd[i] + 0.5f * (vzSrc[i] - back[i]);
                    int gx = std::clamp(int(posFwd.x / dx - 0.5f), 0, n - 1);
                    int gy = std::clamp(int(posFwd.y / dx - 0.5f), 0, n - 1);
                    int gz = std::clamp(int(posFwd.z / dx), 0, n);
                    float fMin = vzSrc[idZ(gx, gy, gz)], fMax = fMin;
                    for (int dz = 0; dz <= 1; ++dz)
                        for (int dy = 0; dy <= 1; ++dy)
                            for (int dx2 = 0; dx2 <= 1; ++dx2) {
                                int a = std::clamp(gx + dx2, 0, n - 1);
                                int b = std::clamp(gy + dy, 0, n - 1);
                                int c = std::clamp(gz + dz, 0, n);
                                float val = vzSrc[idZ(a, b, c)];
                                fMin = std::min(fMin, val);
                                fMax = std::max(fMax, val);
                            }
                    if (corrected < fMin || corrected > fMax) out[i] = fwd[i];
                    else out[i] = corrected;
                }
    }

    void applyBoundary() {
        if (!openXm) for (int z = 0; z < n; ++z) for (int y = 0; y < n; ++y) vx[idX(0, y, z)] = 0.f;
        if (!openXp) for (int z = 0; z < n; ++z) for (int y = 0; y < n; ++y) vx[idX(n, y, z)] = 0.f;
        if (!openYm) for (int z = 0; z < n; ++z) for (int x = 0; x < n; ++x) vy[idY(x, 0, z)] = 0.f;
        if (!openYp) for (int z = 0; z < n; ++z) for (int x = 0; x < n; ++x) vy[idY(x, n, z)] = 0.f;
        if (!openZm) for (int y = 0; y < n; ++y) for (int x = 0; x < n; ++x) vz[idZ(x, y, 0)] = 0.f;
        if (!openZp) for (int y = 0; y < n; ++y) for (int x = 0; x < n; ++x) vz[idZ(x, y, n)] = 0.f;
    }

    void project(float dt, int iterations) {
        const float hInv = 1.0f / h(), h2 = h() * h();
        applyBoundary();
        std::fill(divergence.begin(), divergence.end(), 0.f);
        for (int z = 0; z < n; ++z)
            for (int y = 0; y < n; ++y)
                for (int x = 0; x < n; ++x)
                    divergence[idC(x, y, z)] = (vx[idX(x + 1, y, z)] - vx[idX(x, y, z)] +
                        vy[idY(x, y + 1, z)] - vy[idY(x, y, z)] +
                        vz[idZ(x, y, z + 1)] - vz[idZ(x, y, z)]) * hInv;

        std::fill(pressure.begin(), pressure.end(), 0.f);
        const int dx[6] = { -1, 1, 0, 0, 0, 0 }, dy[6] = { 0, 0, -1, 1, 0, 0 }, dz[6] = { 0, 0, 0, 0, -1, 1 };
        const float omega = sorOmega;

        for (int iter = 0; iter < iterations; ++iter) {
            for (int parity = 0; parity < 2; ++parity) {
                for (int z = 0; z < n; ++z)
                    for (int y = 0; y < n; ++y)
                        for (int x = 0; x < n; ++x) {
                            if (((x + y + z) & 1) != parity) continue;
                            const int i = idC(x, y, z);
                            float sum = 0.f;
                            int terms = 0;
                            for (int k = 0; k < 6; ++k) {
                                int a = x + dx[k], b = y + dy[k], c = z + dz[k];
                                if (!inside(a, b, c)) continue;
                                ++terms;
                                sum += pressure[idC(a, b, c)];
                            }
                            pressure[i] = (1.f - omega) * pressure[i] +
                                omega * (sum - divergence[i] * h2 / dt) / std::max(terms, 1);
                        }
            }
        }

        for (int z = 0; z < n; ++z)
            for (int y = 0; y < n; ++y)
                for (int x = 0; x <= n; ++x) {
                    if ((x == 0 && !openXm) || (x == n && !openXp)) continue;
                    float pRight = (x == n) ? 0.f : pressure[idC(x, y, z)];
                    float pLeft = (x == 0) ? 0.f : pressure[idC(x - 1, y, z)];
                    vx[idX(x, y, z)] -= dt * (pRight - pLeft) * hInv;
                }
        for (int z = 0; z < n; ++z)
            for (int y = 0; y <= n; ++y)
                for (int x = 0; x < n; ++x) {
                    if ((y == 0 && !openYm) || (y == n && !openYp)) continue;
                    float pTop = (y == n) ? 0.f : pressure[idC(x, y, z)];
                    float pBottom = (y == 0) ? 0.f : pressure[idC(x, y - 1, z)];
                    vy[idY(x, y, z)] -= dt * (pTop - pBottom) * hInv;
                }
        for (int z = 0; z <= n; ++z)
            for (int y = 0; y < n; ++y)
                for (int x = 0; x < n; ++x) {
                    if ((z == 0 && !openZm) || (z == n && !openZp)) continue;
                    float pFront = (z == n) ? 0.f : pressure[idC(x, y, z)];
                    float pBack = (z == 0) ? 0.f : pressure[idC(x, y, z - 1)];
                    vz[idZ(x, y, z)] -= dt * (pFront - pBack) * hInv;
                }
        applyBoundary();

        if (openXm) for (int z = 0; z < n; ++z) for (int y = 0; y < n; ++y) vx[idX(0, y, z)] = std::min(vx[idX(0, y, z)], 0.f);
        if (openXp) for (int z = 0; z < n; ++z) for (int y = 0; y < n; ++y) vx[idX(n, y, z)] = std::max(vx[idX(n, y, z)], 0.f);
        if (openYm) for (int z = 0; z < n; ++z) for (int x = 0; x < n; ++x) vy[idY(x, 0, z)] = std::min(vy[idY(x, 0, z)], 0.f);
        if (openYp) for (int z = 0; z < n; ++z) for (int x = 0; x < n; ++x) vy[idY(x, n, z)] = std::max(vy[idY(x, n, z)], 0.f);
        if (openZm) for (int y = 0; y < n; ++y) for (int x = 0; x < n; ++x) vz[idZ(x, y, 0)] = std::min(vz[idZ(x, y, 0)], 0.f);
        if (openZp) for (int y = 0; y < n; ++y) for (int x = 0; x < n; ++x) vz[idZ(x, y, n)] = std::max(vz[idZ(x, y, n)], 0.f);
    }

    void emit(float strength, float dt) {
        if (strength <= 0.f) return;
        const float dx = h();
        const float dtScale = dt * 60.0f;
        const float ecx = emitterCenterX * boxSize, ecy = emitterCenterY * boxSize, ecz = emitterCenterZ * boxSize;
        const float erad = emitterRadius * boxSize;
        const int emitR = int(erad / dx + 1);
        const int cx = int(ecx / dx), cy = int(ecy / dx), cz = int(ecz / dx);
        for (int z = std::max(0, cz - emitR); z < std::min(n, cz + emitR); ++z)
            for (int y = std::max(0, cy - emitR); y < std::min(n, cy + emitR); ++y)
                for (int x = std::max(0, cx - emitR); x < std::min(n, cx + emitR); ++x) {
                    Vec3 p{ (x + .5f) * dx, (y + .5f) * dx, (z + .5f) * dx };
                    float qx = (p.x - ecx) / erad, qy = (p.y - ecy) / erad, qz = (p.z - ecz) / erad;
                    float q = qx * qx + qy * qy + qz * qz;
                    if (q > 1.0f) continue;
                    const int i = idC(x, y, z);
                    smoke[i] = 1.0f;
                    const float kick = emitterKickBase + emitterKickScale * strength;
                    vz[idZ(x, y, z)] = std::max(vz[idZ(x, y, z)], kick);
                    vz[idZ(x, y, z + 1)] = std::max(vz[idZ(x, y, z + 1)], kick);
                    const float j = (frand() - 0.5f) * stirStrength * dtScale;
                    vx[idX(x, y, z)] += j; vx[idX(x + 1, y, z)] += j;
                    vy[idY(x, y, z)] += j; vy[idY(x, y + 1, z)] += j;
                }
    }

    void applyBuoyancy(float dtStep, float buoyancyVal) {
        for (int z = 0; z < n; ++z)
            for (int y = 0; y < n; ++y)
                for (int x = 0; x < n; ++x) {
                    const float f = buoyancyVal * smoke[idC(x, y, z)] * dtStep;
                    vz[idZ(x, y, z)] += buoyancySplit * f;
                    vz[idZ(x, y, z + 1)] += buoyancySplit * f;
                }
    }

    void step(float dt, float buoyancy, float sourceStrength, int projectIterations) {
        const float halfDt = 0.5f * dt;
        const int halfIters = std::max(1, projectIterations / 2);

        applyBuoyancy(halfDt, buoyancy);
        emit(sourceStrength, halfDt);
        advectScalar(smoke, halfDt);

        std::vector<float> vx0 = vx, vy0 = vy, vz0 = vz;
        std::vector<float> vxOld = vx, vyOld = vy, vzOld = vz;
        advectVx(halfDt, vxOld, vxOld, vyOld, vzOld, vx);
        std::vector<float> vxTilde = vx;
        vx = vx0; vy = vy0; vz = vz0;
        vxOld = vx; vyOld = vy; vzOld = vz;
        advectVy(halfDt, vyOld, vxOld, vyOld, vzOld, vy);
        std::vector<float> vyTilde = vy;
        vx = vx0; vy = vy0; vz = vz0;
        vxOld = vx; vyOld = vy; vzOld = vz;
        advectVz(halfDt, vzOld, vxOld, vyOld, vzOld, vz);
        std::vector<float> vzTilde = vz;
        vx = vxTilde; vy = vyTilde;

        project(halfDt, halfIters);

        std::vector<float> vxMid = vx, vyMid = vy, vzMid = vz;
        std::vector<float> vxHat = vxMid, vyHat = vyMid, vzHat = vzMid;
        for (size_t i = 0; i < vx.size(); ++i) vxHat[i] = 2.f * vxMid[i] - vxTilde[i];
        for (size_t i = 0; i < vy.size(); ++i) vyHat[i] = 2.f * vyMid[i] - vyTilde[i];
        for (size_t i = 0; i < vz.size(); ++i) vzHat[i] = 2.f * vzMid[i] - vzTilde[i];

        advectScalar(smoke, halfDt);
        advectVx(halfDt, vxHat, vxMid, vyMid, vzMid, vx, macCormackVel);
        advectVy(halfDt, vyHat, vxMid, vyMid, vzMid, vy, macCormackVel);
        advectVz(halfDt, vzHat, vxMid, vyMid, vzMid, vz, macCormackVel);

        applyBuoyancy(halfDt, buoyancy);
        emit(sourceStrength, halfDt);
        project(halfDt, halfIters);
    }

    void stepNoReflect(float dt, float buoyancy, float sourceStrength, int projectIterations) {
        advectScalar(smoke, dt);
        {
            std::vector<float> vxOld = vx, vyOld = vy, vzOld = vz;
            advectVx(dt, vxOld, vxOld, vyOld, vzOld, vx);
            advectVy(dt, vyOld, vxOld, vyOld, vzOld, vy);
            advectVz(dt, vzOld, vxOld, vyOld, vzOld, vz);
        }
        applyBuoyancy(dt, buoyancy);
        emit(sourceStrength, dt);
        project(dt, projectIterations);
    }

    void rebox(float oldSz, float newSz) {
        if (newSz == oldSz) return;
        boxSize = newSz;
        std::fill(smoke.begin(), smoke.end(), 0.f);
        std::fill(vx.begin(), vx.end(), 0.f);
        std::fill(vy.begin(), vy.end(), 0.f);
        std::fill(vz.begin(), vz.end(), 0.f);
        std::fill(pressure.begin(), pressure.end(), 0.f);
        std::fill(divergence.begin(), divergence.end(), 0.f);
    }
};

struct VolumeRenderer {
    GLuint prog = 0, vao = 0, vbo = 0, volumeTex = 0;
    int n = 0;

    void init(int res, float boxSz) {
        n = res;
        const float s = boxSz;
        float cube[] = {
            // -Z face
            0,0,0,  0,s,0,  s,s,0,  0,0,0,  s,s,0,  s,0,0,
            // +Z face
            0,0,s,  s,0,s,  s,s,s,  0,0,s,  s,s,s,  0,s,s,
            // -Y face
            0,0,0,  s,0,0,  s,0,s,  0,0,0,  s,0,s,  0,0,s,
            // +Y face
            0,s,0,  0,s,s,  s,s,s,  0,s,0,  s,s,s,  s,s,0,
            // -X face
            0,0,0,  0,0,s,  0,s,s,  0,0,0,  0,s,s,  0,s,0,
            // +X face
            s,0,0,  s,s,0,  s,s,s,  s,0,0,  s,s,s,  s,0,s,
        };
        glGenVertexArrays(1, &vao);
        glGenBuffers(1, &vbo);
        glBindVertexArray(vao);
        glBindBuffer(GL_ARRAY_BUFFER, vbo);
        glBufferData(GL_ARRAY_BUFFER, sizeof(cube), cube, GL_STATIC_DRAW);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
        glEnableVertexAttribArray(0);

        glGenTextures(1, &volumeTex);
        glBindTexture(GL_TEXTURE_3D, volumeTex);
        glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);

        const char* vert = R"glsl(#version 330 core
            layout(location=0) in vec3 aPos;
            uniform mat4 u_MVP;
            out vec3 vWorldPos;
            void main() {
                vWorldPos = aPos;
                gl_Position = u_MVP * vec4(aPos, 1.0);
            }
        )glsl";
        const char* frag = R"glsl(#version 330 core
            in vec3 vWorldPos;
            out vec4 FragColor;
            uniform sampler3D u_Volume;
            uniform float u_StepScale;
            uniform float u_BoxSize;
            uniform float u_GridRes;
            uniform float u_AlphaMul;
            uniform vec3 u_CamPos;
            void main() {
                vec3 ro = u_CamPos;
                vec3 rd = normalize(vWorldPos - ro);
                vec3 t0 = (vec3(0.0) - ro) / rd;
                vec3 t1 = (vec3(u_BoxSize) - ro) / rd;
                vec3 tn = min(t0, t1), tf = max(t0, t1);
                float tnear = max(max(tn.x, tn.y), tn.z);
                float tfar  = min(min(tf.x, tf.y), tf.z);
                if (tnear < 0.0) tnear = 0.0;
                if (tnear >= tfar) discard;
                float step = u_StepScale / u_GridRes;
                vec3 bg = vec3(0.03, 0.04, 0.07);
                vec3 smokeCol = vec3(0.9, 0.85, 0.75);
                vec4 col = vec4(0.0);
                bool hitSmoke = false;
                for (float t = tnear; t < tfar; t += step) {
                    float d = texture(u_Volume, (ro + rd * t) / u_BoxSize).r;
                    if (d > 0.001) hitSmoke = true;
                    float alpha = clamp(d * step * u_AlphaMul, 0.0, 1.0);
                    col.rgb += (1.0 - col.a) * smokeCol * alpha;
                    col.a += (1.0 - col.a) * alpha;
                    if (col.a > 0.99) break;
                }
                if (hitSmoke) {
                    FragColor = vec4(bg * (1.0 - col.a) + col.rgb, 1.0);
                } else {
                    discard;
                }
            }
        )glsl";

        auto compile = [](GLenum type, const char* s) {
            GLuint sh = glCreateShader(type);
            glShaderSource(sh, 1, &s, nullptr); glCompileShader(sh);
            int ok = 0; glGetShaderiv(sh, GL_COMPILE_STATUS, &ok);
            if (!ok) { char log[1024]; glGetShaderInfoLog(sh, 1024, nullptr, log); std::fprintf(stderr, "GLSL %s compile: %s\n", type == GL_VERTEX_SHADER ? "VS" : "FS", log); }
            return sh;
            };
        GLuint vs = compile(GL_VERTEX_SHADER, vert), fs = compile(GL_FRAGMENT_SHADER, frag);
        prog = glCreateProgram();
        glAttachShader(prog, vs); glAttachShader(prog, fs); glLinkProgram(prog);
        int linkOk = 0; glGetProgramiv(prog, GL_LINK_STATUS, &linkOk);
        if (!linkOk) { char log[1024]; glGetProgramInfoLog(prog, 1024, nullptr, log); std::fprintf(stderr, "Link: %s\n", log); }
        glDeleteShader(vs); glDeleteShader(fs);
    }

    void setBoxSize(float boxSz) {
        const float s = boxSz;
        float cube[] = {
            // -Z face
            0,0,0,  0,s,0,  s,s,0,  0,0,0,  s,s,0,  s,0,0,
            // +Z face
            0,0,s,  s,0,s,  s,s,s,  0,0,s,  s,s,s,  0,s,s,
            // -Y face
            0,0,0,  s,0,0,  s,0,s,  0,0,0,  s,0,s,  0,0,s,
            // +Y face
            0,s,0,  0,s,s,  s,s,s,  0,s,0,  s,s,s,  s,s,0,
            // -X face
            0,0,0,  0,0,s,  0,s,s,  0,0,0,  0,s,s,  0,s,0,
            // +X face
            s,0,0,  s,s,0,  s,s,s,  s,0,0,  s,s,s,  s,0,s,
        };
        glBindBuffer(GL_ARRAY_BUFFER, vbo);
        glBufferData(GL_ARRAY_BUFFER, sizeof(cube), cube, GL_STATIC_DRAW);
    }

    void upload(const std::vector<float>& density) {
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_3D, volumeTex);
        glTexImage3D(GL_TEXTURE_3D, 0, GL_R16F, n, n, n, 0, GL_RED, GL_FLOAT, density.data());
    }

    void render(int w, int h, const float* mvp, float cx, float cy, float cz, float stepScale, float boxSize, float alphaMul, FILE* logFile) {
        glViewport(0, 0, w, h);
        glDisable(GL_DEPTH_TEST);
        glEnable(GL_CULL_FACE);
        glCullFace(GL_FRONT);
        glUseProgram(prog);
        glUniformMatrix4fv(glGetUniformLocation(prog, "u_MVP"), 1, GL_FALSE, mvp);
        glUniform3f(glGetUniformLocation(prog, "u_CamPos"), cx, cy, cz);
        glUniform1f(glGetUniformLocation(prog, "u_StepScale"), stepScale);
        glUniform1f(glGetUniformLocation(prog, "u_BoxSize"), boxSize);
        glUniform1f(glGetUniformLocation(prog, "u_GridRes"), float(n));
        glUniform1f(glGetUniformLocation(prog, "u_AlphaMul"), alphaMul);
        glUniform1i(glGetUniformLocation(prog, "u_Volume"), 0);
        glBindTexture(GL_TEXTURE_3D, volumeTex);
        glBindVertexArray(vao);
        glDrawArrays(GL_TRIANGLES, 0, 36);
        glBindVertexArray(0);
        glUseProgram(0);
        glDisable(GL_CULL_FACE);
        glCullFace(GL_BACK);
    }

    void shutdown() {
        if (prog) glDeleteProgram(prog);
        if (vao) glDeleteVertexArrays(1, &vao);
        if (vbo) glDeleteBuffers(1, &vbo);
        if (volumeTex) glDeleteTextures(1, &volumeTex);
    }
};

static void buildMVP(float* m, float& cx, float& cy, float& cz, float azimuth, float elevation, float dist, float aspect, float boxSize) {
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

int main() {
    if (!glfwInit()) return 1;
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    GLFWwindow* window = glfwCreateWindow(960, 960, "vapor", nullptr, nullptr);
    if (!window) return 1;
    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGui::StyleColorsDark();
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init("#version 330");
    if (!loadOpenGLFunctions()) { glfwDestroyWindow(window); glfwTerminate(); return 1; }

    SmokeSim3D sim(32);
    VolumeRenderer renderer;

    FILE* logFile = std::fopen("vapor_console.log", "w");
    auto logPrint = [&](const char* fmt, ...) {
        va_list args1, args2;
        va_start(args1, fmt); va_start(args2, fmt);
        std::vprintf(fmt, args1);
        if (logFile) { std::vfprintf(logFile, fmt, args2); std::fflush(logFile); }
        va_end(args2); va_end(args1);
        };

    renderer.init(sim.n, sim.boxSize);

    bool paused = false, useReflection = false, debugPrintOn = false;
    float buoyancy = 5.0f, sourceStrength = 1.0f;
    int projectIterations = sim.n * 2;
    float alphaMul = 30.0f;
    int debugFrame = 0;
    const float maxDt = 0.033f;
    float azimuth = -1.2f, elevation = 0.3f, dist = sim.boxSize * 2.5f, stepScale = 1.5f;
    float prevBoxSize = sim.boxSize;
    bool dragging = false;
    double lastX = 0, lastY = 0;

    while (!glfwWindowShouldClose(window)) {
        glfwPollEvents();
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        if (!paused) {
            float dt = std::min(ImGui::GetIO().DeltaTime, maxDt);
            if (useReflection) sim.step(dt, buoyancy, sourceStrength, projectIterations);
            else sim.stepNoReflect(dt, buoyancy, sourceStrength, projectIterations);
            if (debugPrintOn) {
                ++debugFrame;
                float dmax = *std::max_element(sim.smoke.begin(), sim.smoke.end());
                float vmax = 0.f;
                for (float v : sim.vz) vmax = std::max(vmax, std::abs(v));
                logPrint("f %d maxSmoke=%.3f maxVz=%.2f\n", debugFrame, dmax, vmax);
            }
        }

        int w, h; glfwGetFramebufferSize(window, &w, &h);
        glBindFramebuffer(GL_FRAMEBUFFER, 0);
        glClearColor(0.03f, 0.04f, 0.07f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        renderer.upload(sim.smoke);
        float mvp[16], cx, cy, cz;
        buildMVP(mvp, cx, cy, cz, azimuth, elevation, dist, float(w) / float(h), sim.boxSize);
        renderer.render(w, h, mvp, cx, cy, cz, stepScale, sim.boxSize, alphaMul, logFile);

        if (debugPrintOn && debugFrame % 5 == 0) {
            int ascW = 80, ascH = 24;
            std::vector<unsigned char> full(w * h * 3);
            glReadBuffer(GL_BACK);
            glPixelStorei(GL_PACK_ALIGNMENT, 1);
            glReadPixels(0, 0, w, h, GL_RGB, GL_UNSIGNED_BYTE, full.data());
            glPixelStorei(GL_PACK_ALIGNMENT, 4);
            logPrint("  framebuffer %dx%d:\n", w, h);
            for (int ry = ascH - 1; ry >= 0; --ry) {
                logPrint("  ");
                for (int rx = 0; rx < ascW; ++rx) {
                    int sx = rx * w / ascW;
                    int sy = ry * h / ascH;
                    int i = (sy * w + sx) * 3;
                    float lum = 0.2126f * full[i] + 0.7152f * full[i + 1] + 0.0722f * full[i + 2];
                    const char* ramp = " .:-=+*#%@";
                    int idx = int(lum / 255.0f * 9.0f);
                    idx = idx < 0 ? 0 : idx > 9 ? 9 : idx;
                    logPrint("%c", ramp[idx]);
                }
                logPrint("\n");
            }
        }

        ImGui::SetNextWindowPos({ 16, 16 }, ImGuiCond_FirstUseEver);
        ImGui::Begin("Controls", nullptr, ImGuiWindowFlags_AlwaysAutoResize);
        if (ImGui::Button(paused ? "Resume" : "Pause")) paused = !paused;
        ImGui::SameLine();
        if (ImGui::Button("Reset")) {
            sim = SmokeSim3D(sim.n);
            prevBoxSize = sim.boxSize;
            renderer.setBoxSize(sim.boxSize);
            dist = sim.boxSize * 2.5f;
        }
        ImGui::Checkbox("Reflection", &useReflection);
        ImGui::Checkbox("MC smoke", &sim.macCormackSmoke);
        ImGui::Checkbox("MC vel", &sim.macCormackVel);
        ImGui::Checkbox("Debug print", &debugPrintOn);
        ImGui::SliderFloat("Buoyancy", &buoyancy, 0.f, 10.f);
        ImGui::SliderFloat("Source", &sourceStrength, 0.f, 3.f);
        ImGui::SliderFloat("Box size", &sim.boxSize, 0.5f, 5.f);
        ImGui::SliderFloat("Emitt radius", &sim.emitterRadius, 0.05f, 0.5f);
        ImGui::SliderFloat("Stir", &sim.stirStrength, 0.f, 2.f);
        ImGui::SliderFloat("Alpha mul", &alphaMul, 5.f, 100.f);
        ImGui::SliderInt("SOR its", &projectIterations, 10, 500);
        ImGui::SliderFloat("MC CFL", &sim.cflMc, 0.5f, 10.f);
        ImGui::SliderFloat("Step", &stepScale, 0.1f, 3.f);
        ImGui::Text("Drag mouse to orbit, scroll to zoom");
        ImGui::Text("FPS: %.1f", ImGui::GetIO().Framerate);
        ImGui::End();

        if (sim.boxSize != prevBoxSize) {
            sim.rebox(prevBoxSize, sim.boxSize);
            renderer.setBoxSize(sim.boxSize);
            dist = sim.boxSize * 2.5f;
            prevBoxSize = sim.boxSize;
        }

        ImGui::Render();
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

        if (!ImGui::GetIO().WantCaptureMouse && ImGui::GetIO().MouseWheel != 0.f) {
            dist = std::clamp(dist - ImGui::GetIO().MouseWheel * dist * 0.1f, sim.boxSize * 0.3f, sim.boxSize * 8.f);
        }

        if (ImGui::IsMouseDown(0) && !ImGui::GetIO().WantCaptureMouse) {
            double mx, my; glfwGetCursorPos(window, &mx, &my);
            if (!dragging) { dragging = true; lastX = mx; lastY = my; }
            else { azimuth += float((mx - lastX) * 0.005f); elevation = std::clamp(elevation + float((lastY - my) * 0.005f), -1.5f, 1.5f); lastX = mx; lastY = my; }
        }
        else dragging = false;

        glfwSwapBuffers(window);
    }
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
    renderer.shutdown();
    glfwDestroyWindow(window);
    glfwTerminate();
    if (logFile) std::fclose(logFile);
}