#include <algorithm>
#include <cmath>
#include <cstdarg>
#include <cstdint>
#include <cstdio>
#include <numeric>
#include <vector>

#include "gl_loader.h"
#include <GLFW/glfw3.h>
#include <imgui.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_opengl3.h>

struct Vec2 { float x, y; };
static Vec2 operator-(Vec2 a, Vec2 b) { return { a.x - b.x, a.y - b.y }; }
static Vec2 operator+(Vec2 a, Vec2 b) { return { a.x + b.x, a.y + b.y }; }
static Vec2 operator*(Vec2 a, float s) { return { a.x * s, a.y * s }; }

static float frand() {
    static uint32_t state = 0x9E3779B9u;
    state = state * 1664525u + 1013904223u;
    return (state >> 8) * (1.0f / 16777216.0f);
}

struct SmokeSim2D {
    int n = 256;
    int displayMargin = 32;
    float emitterRadius = 0.08f;
    float emitterFalloff = 3.5f;
    float emitterCenterX = 0.5f, emitterCenterY = 0.2f;
    float emitterSmokeScale = 3.0f;
    float emitterKickBase = 0.6f;
    float emitterKickScale = 1.0f;
    float stirStrength = 0.6f;
    float sorOmega = 1.95f;
    float cflMc = 3.0f;
    float buoyancySplit = 0.5f;
    bool macCormackSmoke = true, macCormackVel = true;
    bool openLeft = true, openRight = true, openTop = true, openBottom = false;
    std::vector<float> smoke, vx, vy, pressure, divergence;

    explicit SmokeSim2D(int res) {
        n = res;
        smoke.assign(n * n, 0);
        vx.assign((n + 1) * n, 0);
        vy.assign(n * (n + 1), 0);
        pressure.assign(n * n, 0);
        divergence.assign(n * n, 0);
    }

    int idC(int x, int y) const { return x + n * y; }
    int idX(int x, int y) const { return x + (n + 1) * y; }
    int idY(int x, int y) const { return x + n * y; }
    float h() const { return 1.0f / n; }
    bool inside(int x, int y) const { return x >= 0 && y >= 0 && x < n && y < n; }
    bool insideVx(int x, int y) const { return x >= 0 && y >= 0 && x <= n && y < n; }
    bool insideVy(int x, int y) const { return x >= 0 && y >= 0 && x < n && y <= n; }
    static float mix(float a, float b, float t) { return a + (b - a) * t; }

    float sampleCell(const std::vector<float>& f, Vec2 p) const {
        const float gx = std::clamp(p.x / h() - 0.5f, 0.f, float(n - 1));
        const float gy = std::clamp(p.y / h() - 0.5f, 0.f, float(n - 1));
        const int x0 = int(gx), y0 = int(gy);
        const int x1 = std::min(x0 + 1, n - 1), y1 = std::min(y0 + 1, n - 1);
        const float tx = gx - x0, ty = gy - y0;
        return mix(mix(f[idC(x0, y0)], f[idC(x1, y0)], tx),
                   mix(f[idC(x0, y1)], f[idC(x1, y1)], tx), ty);
    }

    float sampleVxField(const std::vector<float>& f, Vec2 p) const {
        const float gx = std::clamp(p.x / h(), 0.f, float(n));
        const float gy = std::clamp(p.y / h() - 0.5f, 0.f, float(n - 1));
        const int x0 = std::min(int(gx), n - 1), y0 = int(gy);
        const int x1 = x0 + 1, y1 = std::min(y0 + 1, n - 1);
        const float tx = gx - x0, ty = gy - y0;
        return mix(mix(f[idX(x0, y0)], f[idX(x1, y0)], tx),
                   mix(f[idX(x0, y1)], f[idX(x1, y1)], tx), ty);
    }

    float sampleVx(Vec2 p) const { return sampleVxField(vx, p); }

    float sampleVyField(const std::vector<float>& f, Vec2 p) const {
        const float gx = std::clamp(p.x / h() - 0.5f, 0.f, float(n - 1));
        const float gy = std::clamp(p.y / h(), 0.f, float(n));
        const int x0 = int(gx), y0 = std::min(int(gy), n - 1);
        const int x1 = std::min(x0 + 1, n - 1), y1 = y0 + 1;
        const float tx = gx - x0, ty = gy - y0;
        return mix(mix(f[idY(x0, y0)], f[idY(x1, y0)], tx),
                   mix(f[idY(x0, y1)], f[idY(x1, y1)], tx), ty);
    }

    float sampleVy(Vec2 p) const { return sampleVyField(vy, p); }

    Vec2 velocity(Vec2 p) const { return { sampleVx(p), sampleVy(p) }; }
    Vec2 velocity(Vec2 p, const std::vector<float>& vxSrc, const std::vector<float>& vySrc) const {
        return { sampleVxField(vxSrc, p), sampleVyField(vySrc, p) };
    }

    float semiLagrangian(const std::vector<float>& source, Vec2 p, float dt) const {
        return sampleCell(source, p - velocity(p) * dt);
    }

    void advectScalar(std::vector<float>& field, float dt) {
        std::vector<float> fwd(n * n);
        const float dx = h();
        for (int y = 0; y < n; ++y)
            for (int x = 0; x < n; ++x)
                fwd[idC(x, y)] = semiLagrangian(field, { (x + .5f) * dx, (y + .5f) * dx }, dt);

        if (!macCormackSmoke) { field.swap(fwd); return; }

        std::vector<float> back(n * n);
        for (int y = 0; y < n; ++y)
            for (int x = 0; x < n; ++x) {
                Vec2 p{ (x + .5f) * dx, (y + .5f) * dx };
                back[idC(x, y)] = sampleCell(fwd, p + velocity(p) * dt);
            }

        for (int y = 0; y < n; ++y)
            for (int x = 0; x < n; ++x) {
                const int i = idC(x, y);
                const Vec2 p{ (x + .5f) * dx, (y + .5f) * dx };
                const Vec2 posFwd = p - velocity(p) * dt;
                const Vec2 posBwd = p + velocity(p) * dt;
                if (posFwd.x < 0.f || posFwd.x > 1.f || posFwd.y < 0.f || posFwd.y > 1.f ||
                    posBwd.x < 0.f || posBwd.x > 1.f || posBwd.y < 0.f || posBwd.y > 1.f) {
                    field[i] = fwd[i];
                    continue;
                }
                float corrected = fwd[i] + 0.5f * (field[i] - back[i]);
                float fMin = field[i], fMax = field[i];
                for (int dy = -1; dy <= 1; ++dy)
                    for (int dx2 = -1; dx2 <= 1; ++dx2) {
                        const int a = x + dx2, b = y + dy;
                        if (inside(a, b)) {
                            fMin = std::min(fMin, field[idC(a, b)]);
                            fMax = std::max(fMax, field[idC(a, b)]);
                        }
                    }
                if (corrected < fMin || corrected > fMax)
                    field[i] = fwd[i];
                else
                    field[i] = corrected;
            }
    }

    void advectVx(float dt, const std::vector<float>& vxSrc,
                  const std::vector<float>& vxVel, const std::vector<float>& vyVel,
                  std::vector<float>& out, bool useMC = true) {
        std::vector<float> fwd(out.size());
        const float dx = h();
        for (int y = 0; y < n; ++y)
            for (int x = 0; x <= n; ++x) {
                const Vec2 face{ x * dx, (y + .5f) * dx };
                fwd[idX(x, y)] = sampleVxField(vxSrc, face - velocity(face, vxVel, vyVel) * dt);
            }

        if (!useMC || !macCormackVel) { out.swap(fwd); return; }

        std::vector<float> back(out.size());
        for (int y = 0; y < n; ++y)
            for (int x = 0; x <= n; ++x) {
                const Vec2 face{ x * dx, (y + .5f) * dx };
                back[idX(x, y)] = sampleVxField(fwd, face + velocity(face, vxVel, vyVel) * dt);
            }

        for (int y = 0; y < n; ++y)
            for (int x = 0; x <= n; ++x) {
                const int i = idX(x, y);
                const Vec2 face{ x * dx, (y + .5f) * dx };
                const Vec2 v = velocity(face, vxVel, vyVel);
                const Vec2 posFwd = face - v * dt;
                const Vec2 posBwd = face + v * dt;
                if (posFwd.x < 0.f || posFwd.x > 1.f || posFwd.y < 0.f || posFwd.y > 1.f ||
                    posBwd.x < 0.f || posBwd.x > 1.f || posBwd.y < 0.f || posBwd.y > 1.f ||
                    std::sqrt(v.x * v.x + v.y * v.y) * dt / dx > cflMc) {
                    out[i] = fwd[i];
                    continue;
                }
                float corrected = fwd[i] + 0.5f * (vxSrc[i] - back[i]);
                float fMin = vxSrc[i], fMax = vxSrc[i];
                for (int dy = -1; dy <= 1; ++dy)
                    for (int dx2 = -1; dx2 <= 1; ++dx2) {
                        const int a = x + dx2, b = y + dy;
                        if (insideVx(a, b)) {
                            fMin = std::min(fMin, vxSrc[idX(a, b)]);
                            fMax = std::max(fMax, vxSrc[idX(a, b)]);
                        }
                    }
                if (corrected < fMin || corrected > fMax)
                    out[i] = fwd[i];
                else
                    out[i] = corrected;
            }
    }

    void advectVy(float dt, const std::vector<float>& vySrc,
                  const std::vector<float>& vxVel, const std::vector<float>& vyVel,
                  std::vector<float>& out, bool useMC = true) {
        std::vector<float> fwd(out.size());
        const float dx = h();
        for (int y = 0; y <= n; ++y)
            for (int x = 0; x < n; ++x) {
                const Vec2 face{ (x + .5f) * dx, y * dx };
                fwd[idY(x, y)] = sampleVyField(vySrc, face - velocity(face, vxVel, vyVel) * dt);
            }

        if (!useMC || !macCormackVel) { out.swap(fwd); return; }

        std::vector<float> back(out.size());
        for (int y = 0; y <= n; ++y)
            for (int x = 0; x < n; ++x) {
                const Vec2 face{ (x + .5f) * dx, y * dx };
                back[idY(x, y)] = sampleVyField(fwd, face + velocity(face, vxVel, vyVel) * dt);
            }

        for (int y = 0; y <= n; ++y)
            for (int x = 0; x < n; ++x) {
                const int i = idY(x, y);
                const Vec2 face{ (x + .5f) * dx, y * dx };
                const Vec2 v = velocity(face, vxVel, vyVel);
                const Vec2 posFwd = face - v * dt;
                const Vec2 posBwd = face + v * dt;
                if (posFwd.x < 0.f || posFwd.x > 1.f || posFwd.y < 0.f || posFwd.y > 1.f ||
                    posBwd.x < 0.f || posBwd.x > 1.f || posBwd.y < 0.f || posBwd.y > 1.f ||
                    std::sqrt(v.x * v.x + v.y * v.y) * dt / dx > cflMc) {
                    out[i] = fwd[i];
                    continue;
                }
                float corrected = fwd[i] + 0.5f * (vySrc[i] - back[i]);
                float fMin = vySrc[i], fMax = vySrc[i];
                for (int dy = -1; dy <= 1; ++dy)
                    for (int dx2 = -1; dx2 <= 1; ++dx2) {
                        const int a = x + dx2, b = y + dy;
                        if (insideVy(a, b)) {
                            fMin = std::min(fMin, vySrc[idY(a, b)]);
                            fMax = std::max(fMax, vySrc[idY(a, b)]);
                        }
                    }
                if (corrected < fMin || corrected > fMax)
                    out[i] = fwd[i];
                else
                    out[i] = corrected;
            }
    }

    void applyBoundary() {
        if (!openLeft)
            for (int y = 0; y < n; ++y) vx[idX(0, y)] = 0.f;
        if (!openRight)
            for (int y = 0; y < n; ++y) vx[idX(n, y)] = 0.f;
        if (!openBottom)
            for (int x = 0; x < n; ++x) vy[idY(x, 0)] = 0.f;
        if (!openTop)
            for (int x = 0; x < n; ++x) vy[idY(x, n)] = 0.f;
    }

    void project(float dt, int iterations) {
        const float hInv = 1.0f / h(), h2 = h() * h();

        applyBoundary();
        std::fill(divergence.begin(), divergence.end(), 0.f);
        for (int y = 0; y < n; ++y)
            for (int x = 0; x < n; ++x)
                divergence[idC(x, y)] = (vx[idX(x + 1, y)] - vx[idX(x, y)] +
                                         vy[idY(x, y + 1)] - vy[idY(x, y)]) * hInv;

        std::fill(pressure.begin(), pressure.end(), 0.f);
        const int dx[4] = { -1, 1, 0, 0 }, dy[4] = { 0, 0, -1, 1 };
        const float omega = sorOmega;

        for (int iter = 0; iter < iterations; ++iter) {
            for (int parity = 0; parity < 2; ++parity) {
                for (int y = 0; y < n; ++y)
                    for (int x = 0; x < n; ++x) {
                        if (((x + y) & 1) != parity) continue;
                        const int i = idC(x, y);
                        float sum = 0.f;
                        int terms = 0;
                        for (int k = 0; k < 4; ++k) {
                            const int a = x + dx[k], b = y + dy[k];
                            if (!inside(a, b)) continue;
                            ++terms;
                            sum += pressure[idC(a, b)];
                        }
                        pressure[i] = (1.f - omega) * pressure[i] +
                                      omega * (sum - divergence[i] * h2 / dt) / std::max(terms, 1);
                    }
            }
        }

        for (int y = 0; y < n; ++y)
            for (int x = 1; x < n; ++x)
                vx[idX(x, y)] -= dt * (pressure[idC(x, y)] - pressure[idC(x - 1, y)]) * hInv;
        for (int y = 1; y < n; ++y)
            for (int x = 0; x < n; ++x)
                vy[idY(x, y)] -= dt * (pressure[idC(x, y)] - pressure[idC(x, y - 1)]) * hInv;
        applyBoundary();

        if (openLeft)
            for (int y = 0; y < n; ++y) vx[idX(0, y)] = std::min(vx[idX(0, y)], 0.f);
        if (openRight)
            for (int y = 0; y < n; ++y) vx[idX(n, y)] = std::max(vx[idX(n, y)], 0.f);
        if (openBottom)
            for (int x = 0; x < n; ++x) vy[idY(x, 0)] = std::min(vy[idY(x, 0)], 0.f);
        if (openTop)
            for (int x = 0; x < n; ++x) vy[idY(x, n)] = std::max(vy[idY(x, n)], 0.f);
    }

    void emit(float strength, float dt) {
        if (strength <= 0.f) return;
        const float r2 = emitterRadius * emitterRadius, dx = h();
        const float dtScale = dt * 60.0f;
        const int emitH = std::min(int(emitterCenterY * n + emitterRadius * n + 2), n);
        for (int y = std::max(0, int(emitterCenterY * n - emitterRadius * n - 1)); y < emitH; ++y)
            for (int x = 0; x < n; ++x) {
                Vec2 p{ (x + .5f) * dx, (y + .5f) * dx };
                const float dx2 = (p.x - emitterCenterX) * (p.x - emitterCenterX);
                const float dy2 = (p.y - emitterCenterY) * (p.y - emitterCenterY);
                const float q = (dx2 + dy2) / r2;
                if (q > 1.0f) continue;
                const float w = std::exp(-emitterFalloff * q);
                const int i = idC(x, y);
                smoke[i] = std::max(smoke[i], std::min(1.0f, strength * emitterSmokeScale * w));
                const float kick = (emitterKickBase + emitterKickScale * strength) * w;
                vy[idY(x, y)] = std::max(vy[idY(x, y)], kick);
                vy[idY(x, y + 1)] = std::max(vy[idY(x, y + 1)], kick);
                const float j = (frand() - 0.5f) * stirStrength * w * dtScale;
                vx[idX(x, y)] += j;
                vx[idX(x + 1, y)] += j;
            }
    }

    void step(float dt, float buoyancy, float sourceStrength, int projectIterations) {
        const float halfDt = 0.5f * dt;
        const int halfIters = std::max(1, projectIterations / 2);
        const float dxInv = 1.0f / h();

        // ---- Stage 1: forward half-step of smoke using u₀ ----
        advectScalar(smoke, halfDt);

        for (int y = 0; y < n; ++y) {
            if (openLeft) {
                const float vf = -std::min(vx[idX(0, y)], 0.f);
                smoke[idC(0, y)] = std::max(0.f, smoke[idC(0, y)] * (1.f - vf * halfDt * dxInv));
            }
            if (openRight) {
                const float vf = std::max(vx[idX(n, y)], 0.f);
                smoke[idC(n - 1, y)] = std::max(0.f, smoke[idC(n - 1, y)] * (1.f - vf * halfDt * dxInv));
            }
        }
        for (int x = 0; x < n; ++x) {
            if (openBottom) {
                const float vf = -std::min(vy[idY(x, 0)], 0.f);
                smoke[idC(x, 0)] = std::max(0.f, smoke[idC(x, 0)] * (1.f - vf * halfDt * dxInv));
            }
            if (openTop) {
                const float vf = std::max(vy[idY(x, n)], 0.f);
                smoke[idC(x, n - 1)] = std::max(0.f, smoke[idC(x, n - 1)] * (1.f - vf * halfDt * dxInv));
            }
        }

        // ---- Stage 2: forward advection of velocity, forces, project → u₁/₂ ----
        std::vector<float> vx0 = vx, vy0 = vy;

        std::vector<float> vxOld = vx, vyOld = vy; advectVx(halfDt, vxOld, vxOld, vyOld, vx);
        std::vector<float> vxTilde = vx;

        vx = vx0; vy = vy0;
        vxOld = vx; vyOld = vy; advectVy(halfDt, vyOld, vxOld, vyOld, vy);
        std::vector<float> vyTilde = vy;

        vx = vxTilde;  // now vx,vy = ũ₁/₂ (forward-advected)

        for (int y = 0; y < n; ++y)
            for (int x = 0; x < n; ++x) {
                const float f = buoyancy * smoke[idC(x, y)] * halfDt;
                vy[idY(x, y)] += buoyancySplit * f;
                vy[idY(x, y + 1)] += buoyancySplit * f;
            }

        emit(sourceStrength, dt);
        project(halfDt, halfIters);

        std::vector<float> vxMid = vx, vyMid = vy;  // u₁/₂ (projected midpoint)

        // ---- Stage 3: reflect  û = 2·u₁/₂ − ũ₁/₂ ----
        std::vector<float> vxHat = vxMid, vyHat = vyMid;
        for (size_t i = 0; i < vx.size(); ++i) vxHat[i] = 2.f * vxMid[i] - vxTilde[i];
        for (size_t i = 0; i < vy.size(); ++i) vyHat[i] = 2.f * vyMid[i] - vyTilde[i];

        // ---- Stage 4: advect reflected fields using midpoint velocity, final project ----
        advectScalar(smoke, halfDt);

        for (int y = 0; y < n; ++y) {
            if (openLeft) {
                const float vf = -std::min(vxMid[idX(0, y)], 0.f);
                smoke[idC(0, y)] = std::max(0.f, smoke[idC(0, y)] * (1.f - vf * halfDt * dxInv));
            }
            if (openRight) {
                const float vf = std::max(vxMid[idX(n, y)], 0.f);
                smoke[idC(n - 1, y)] = std::max(0.f, smoke[idC(n - 1, y)] * (1.f - vf * halfDt * dxInv));
            }
        }
        for (int x = 0; x < n; ++x) {
            if (openBottom) {
                const float vf = -std::min(vyMid[idY(x, 0)], 0.f);
                smoke[idC(x, 0)] = std::max(0.f, smoke[idC(x, 0)] * (1.f - vf * halfDt * dxInv));
            }
            if (openTop) {
                const float vf = std::max(vyMid[idY(x, n)], 0.f);
                smoke[idC(x, n - 1)] = std::max(0.f, smoke[idC(x, n - 1)] * (1.f - vf * halfDt * dxInv));
            }
        }

        advectVx(halfDt, vxHat, vxMid, vyMid, vx, macCormackVel);
        advectVy(halfDt, vyHat, vxMid, vyMid, vy, macCormackVel);

        project(halfDt, halfIters);
    }

    void stepNoReflect(float dt, float buoyancy, float sourceStrength, int projectIterations) {
        const float dxInv = 1.0f / h();

        // ---- Advect smoke full dt ----
        advectScalar(smoke, dt);

        for (int y = 0; y < n; ++y) {
            if (openLeft) {
                const float vf = -std::min(vx[idX(0, y)], 0.f);
                smoke[idC(0, y)] = std::max(0.f, smoke[idC(0, y)] * (1.f - vf * dt * dxInv));
            }
            if (openRight) {
                const float vf = std::max(vx[idX(n, y)], 0.f);
                smoke[idC(n - 1, y)] = std::max(0.f, smoke[idC(n - 1, y)] * (1.f - vf * dt * dxInv));
            }
        }
        for (int x = 0; x < n; ++x) {
            if (openBottom) {
                const float vf = -std::min(vy[idY(x, 0)], 0.f);
                smoke[idC(x, 0)] = std::max(0.f, smoke[idC(x, 0)] * (1.f - vf * dt * dxInv));
            }
            if (openTop) {
                const float vf = std::max(vy[idY(x, n)], 0.f);
                smoke[idC(x, n - 1)] = std::max(0.f, smoke[idC(x, n - 1)] * (1.f - vf * dt * dxInv));
            }
        }

        // ---- Advect velocity full dt (MacCormack self-advection) ----
        {
            std::vector<float> vxOld = vx, vyOld = vy;
            advectVx(dt, vxOld, vxOld, vyOld, vx);
            advectVy(dt, vyOld, vxOld, vyOld, vy);
        }

        // ---- Forces ----
        for (int y = 0; y < n; ++y)
            for (int x = 0; x < n; ++x) {
                const float f = buoyancy * smoke[idC(x, y)] * dt;
                vy[idY(x, y)] += buoyancySplit * f;
                vy[idY(x, y + 1)] += buoyancySplit * f;
            }

        emit(sourceStrength, dt);
        project(dt, projectIterations);
    }

    void debugPrint(int frame) const {
        static FILE* logFile = nullptr;
        if (!logFile) { fopen_s(&logFile, "vapor_debug.log", "w"); }
        auto log = [&](const char* fmt, ...) {
            va_list args;
            va_start(args, fmt);
            std::vprintf(fmt, args);
            va_end(args);
            va_start(args, fmt);
            if (logFile) vfprintf(logFile, fmt, args);
            va_end(args);
        };
        log("--- frame %d ---\n", frame);
        log("density  max=%+.4f  min=%+.4f  total=%+.4f\n",
            *std::max_element(smoke.begin(), smoke.end()),
            *std::min_element(smoke.begin(), smoke.end()),
            std::accumulate(smoke.begin(), smoke.end(), 0.0f) / (n * n));

        float divMax = 0.f;
        for (int y = 0; y < n; ++y)
            for (int x = 0; x < n; ++x) {
                float d = (vx[idX(x + 1, y)] - vx[idX(x, y)] +
                           vy[idY(x, y + 1)] - vy[idY(x, y)]) / h();
                if (std::fabs(d) > std::fabs(divMax)) divMax = d;
            }
        log("div max=%+.6f\n", divMax);

        float vMax = 0.f;
        for (float v : vx) if (std::fabs(v) > vMax) vMax = std::fabs(v);
        for (float v : vy) if (std::fabs(v) > vMax) vMax = std::fabs(v);
        log("vel max=%.4f\n", vMax);

        int zeroInPlume = 0, holes = 0;
        for (int y = 0; y < n; ++y)
            for (int x = 0; x < n; ++x) {
                const int i = idC(x, y);
                if (smoke[i] <= 0.f) {
                    bool nearSmoke = false;
                    for (int dy = -2; dy <= 2 && !nearSmoke; ++dy)
                        for (int dx2 = -2; dx2 <= 2; ++dx2)
                            if (inside(x + dx2, y + dy) && smoke[idC(x + dx2, y + dy)] > 0.01f)
                                { nearSmoke = true; break; }
                    if (nearSmoke) ++zeroInPlume;
                }
                if (smoke[i] <= 0.f) {
                    int nc = 0;
                    for (int dy = -1; dy <= 1; ++dy)
                        for (int dx2 = -1; dx2 <= 1; ++dx2)
                            if ((dx2 || dy) && inside(x + dx2, y + dy) && smoke[idC(x + dx2, y + dy)] > 0.1f)
                                ++nc;
                    if (nc >= 5) ++holes;
                }
            }
        log("holes(surrounded)=%d  edges=%d\n\n", holes, zeroInPlume);

        const int cw = 80, ch = 60;
        const char map[] = " .-~=+*#@%";
        const int yLo = 0, yHi = n - 1;
        for (int cy = ch - 1; cy >= 0; --cy) {
            int y = yLo + cy * (yHi - yLo) / (ch - 1);
            for (int cx = 0; cx < cw; ++cx) {
                int x = cx * (n - 1) / (cw - 1);
                float d = smoke[idC(x, y)];
                int idx = std::max(0, std::min(9, int(d * 10.f)));
                std::putchar(map[idx]);
                if (logFile) fputc(map[idx], logFile);
            }
            std::putchar('\n');
            if (logFile) fputc('\n', logFile);
        }
        std::fflush(stdout);
        if (logFile) std::fflush(logFile);
    }
};

struct FieldRenderer {
    GLuint prog = 0, vao = 0, vbo = 0, fieldTex = 0, fbo = 0, renderTex = 0;
    int texWidth = 0, texHeight = 0;

    void init() {
        const float quadv[] = { -1.f, -1.f, 1.f, -1.f, -1.f, 1.f, -1.f, 1.f, 1.f, -1.f, 1.f, 1.f };
        glGenVertexArrays(1, &vao);
        glGenBuffers(1, &vbo);
        glBindVertexArray(vao);
        glBindBuffer(GL_ARRAY_BUFFER, vbo);
        glBufferData(GL_ARRAY_BUFFER, sizeof(quadv), quadv, GL_STATIC_DRAW);
        glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 2 * sizeof(float), (void*)0);
        glEnableVertexAttribArray(0);

        glGenTextures(1, &fieldTex);
        glBindTexture(GL_TEXTURE_2D, fieldTex);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

        const char* vert = "#version 330 core\nlayout(location=0) in vec2 aPos; out vec2 vUV; void main(){ vUV=aPos*0.5+0.5; gl_Position=vec4(aPos,0,1); }";
        const char* frag = R"glsl(
            #version 330 core
            in vec2 vUV; out vec4 FragColor;
            uniform sampler2D u_Field;
            uniform int u_ShowSpeed;
            uniform float u_SpeedMax;
            void main() {
                vec3 bg = vec3(0.03, 0.04, 0.07);
                if (u_ShowSpeed == 0) {
                    float d = texture(u_Field, vUV).r;
                    float a = clamp(d * 2.5, 0.0, 1.0);
                    vec3 col = mix(bg, vec3(0.78, 0.80, 0.87), a);
                    FragColor = vec4(col, 1.0);
                } else {
                    float s = clamp(texture(u_Field, vUV).r / max(u_SpeedMax, 1e-4), 0.0, 1.0);
                    vec3 c0 = vec3(0.05, 0.08, 0.15);
                    vec3 c1 = vec3(0.0, 0.55, 0.85);
                    vec3 c2 = vec3(0.95, 0.75, 0.15);
                    vec3 col = s < 0.5 ? mix(c0, c1, s * 2.0) : mix(c1, c2, (s - 0.5) * 2.0);
                    FragColor = vec4(col, 1.0);
                }
            }
        )glsl";

        auto compile = [](GLenum type, const char* s) {
            GLuint sh = glCreateShader(type);
            glShaderSource(sh, 1, &s, nullptr);
            glCompileShader(sh);
            int ok = 0;
            glGetShaderiv(sh, GL_COMPILE_STATUS, &ok);
            if (!ok) {
                char log[1024];
                glGetShaderInfoLog(sh, 1024, nullptr, log);
                std::fprintf(stderr, "Shader compile error (%s): %s\n", type == GL_VERTEX_SHADER ? "VS" : "FS", log);
            }
            return sh;
        };
        GLuint vs = compile(GL_VERTEX_SHADER, vert), fs = compile(GL_FRAGMENT_SHADER, frag);
        prog = glCreateProgram();
        glAttachShader(prog, vs);
        glAttachShader(prog, fs);
        glLinkProgram(prog);
        int linkOk = 0;
        glGetProgramiv(prog, GL_LINK_STATUS, &linkOk);
        if (!linkOk) {
            char log[1024];
            glGetProgramInfoLog(prog, 1024, nullptr, log);
            std::fprintf(stderr, "Shader link error: %s\n", log);
        }
        glDeleteShader(vs);
        glDeleteShader(fs);
    }

    void shutdown() {
        if (prog)  glDeleteProgram(prog);
        if (vao)   glDeleteVertexArrays(1, &vao);
        if (vbo)   glDeleteBuffers(1, &vbo);
        if (fieldTex) glDeleteTextures(1, &fieldTex);
        if (renderTex) glDeleteTextures(1, &renderTex);
        if (fbo)   glDeleteFramebuffers(1, &fbo);
        prog = vao = vbo = fieldTex = fbo = renderTex = 0;
    }

    void upload(const std::vector<float>& field, int n) {
        glBindTexture(GL_TEXTURE_2D, fieldTex);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_R16F, n, n, 0, GL_RED, GL_FLOAT, field.data());
    }

    void render(int w, int h, bool showSpeed, float speedMax) {
        if (w != texWidth || h != texHeight) {
            texWidth = w;
            texHeight = h;
            if (!fbo) {
                glGenFramebuffers(1, &fbo);
                glGenTextures(1, &renderTex);
            }
            glBindTexture(GL_TEXTURE_2D, renderTex);
            glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, w, h, 0, GL_RGBA, GL_UNSIGNED_BYTE, nullptr);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
            glBindFramebuffer(GL_FRAMEBUFFER, fbo);
            glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, renderTex, 0);
        }
        glBindFramebuffer(GL_FRAMEBUFFER, fbo);
        glViewport(0, 0, w, h);
        glClear(GL_COLOR_BUFFER_BIT);
        glUseProgram(prog);
        glUniform1i(glGetUniformLocation(prog, "u_Field"), 0);
        glUniform1i(glGetUniformLocation(prog, "u_ShowSpeed"), showSpeed ? 1 : 0);
        glUniform1f(glGetUniformLocation(prog, "u_SpeedMax"), speedMax);
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, fieldTex);
        glBindVertexArray(vao);
        glDrawArrays(GL_TRIANGLES, 0, 6);
        glBindFramebuffer(GL_FRAMEBUFFER, 0);
    }
};

int main() {
    if (!glfwInit()) return 1;
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    GLFWwindow* window = glfwCreateWindow(960, 960, "Vapor - 2D smoke plume", nullptr, nullptr);
    if (!window) return 1;
    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGui::StyleColorsDark();
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init("#version 330");
    if (!loadOpenGLFunctions()) {
        std::fprintf(stderr, "Failed to load required OpenGL 3.3 functions.\n");
        glfwDestroyWindow(window);
        glfwTerminate();
        return 1;
    }

    SmokeSim2D sim(256);
    FieldRenderer renderer;
    renderer.init();

    bool paused = false, showSpeed = false, useReflection = false, debugPrintOn = false;
    float buoyancy = 0.8f, sourceStrength = 1.0f;
    int projectIterations = 400;
    int debugFrame = 0;
    const float maxDt = 0.033f;

    while (!glfwWindowShouldClose(window)) {
        glfwPollEvents();
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        ImGui::SetNextWindowPos({ 0, 0 });
        ImGui::SetNextWindowSize(ImGui::GetIO().DisplaySize);
        ImGui::Begin("Sim", nullptr, ImGuiWindowFlags_NoDecoration | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoBringToFrontOnFocus);

        if (!paused) {
            float dt = std::min(ImGui::GetIO().DeltaTime, maxDt);
            if (useReflection)
                sim.step(dt, buoyancy, sourceStrength, projectIterations);
            else
                sim.stepNoReflect(dt, buoyancy, sourceStrength, projectIterations);
            if (debugPrintOn)
                sim.debugPrint(++debugFrame);
        }

        ImVec2 avail = ImGui::GetContentRegionAvail();
        const float s = std::max(1.f, std::min(avail.x, avail.y));
        ImVec2 size(s, s);
        ImGui::SetCursorPos({ std::max(0.f, (avail.x - s) * 0.5f), std::max(0.f, (avail.y - s) * 0.5f) });

        const int m = sim.displayMargin, dispN = sim.n - 2 * m;
        std::vector<float> view(dispN * dispN);
        float speedMax = 1.f;
        if (showSpeed) {
            speedMax = 0.f;
            for (int dy = 0; dy < dispN; ++dy)
                for (int dx = 0; dx < dispN; ++dx) {
                    const int x = dx + m, y = dy + m;
                    const float ux = 0.5f * (sim.vx[sim.idX(x, y)] + sim.vx[sim.idX(x + 1, y)]);
                    const float uy = 0.5f * (sim.vy[sim.idY(x, y)] + sim.vy[sim.idY(x, y + 1)]);
                    const float sp = std::sqrt(ux * ux + uy * uy);
                    view[dy * dispN + dx] = sp;
                    speedMax = std::max(speedMax, sp);
                }
        } else {
            for (int dy = 0; dy < dispN; ++dy)
                for (int dx = 0; dx < dispN; ++dx)
                    view[dy * dispN + dx] = sim.smoke[sim.idC(dx + m, dy + m)];
        }

        renderer.upload(view, dispN);
        renderer.render((int)size.x, (int)size.y, showSpeed, speedMax);
        ImGui::Image((ImTextureID)(intptr_t)renderer.renderTex, size, { 0, 1 }, { 1, 0 });
        ImGui::End();

        ImGui::SetNextWindowPos({ 16, 16 }, ImGuiCond_FirstUseEver);
        ImGui::Begin("Controls", nullptr, ImGuiWindowFlags_AlwaysAutoResize);
        if (ImGui::Button(paused ? "Resume" : "Pause")) paused = !paused;
        ImGui::SameLine();
        if (ImGui::Button("Reset")) {
            bool ol = sim.openLeft, orr = sim.openRight, ot = sim.openTop, ob = sim.openBottom;
            bool ms = sim.macCormackSmoke, mv = sim.macCormackVel;
            sim = SmokeSim2D(sim.n);
            sim.openLeft = ol; sim.openRight = orr; sim.openTop = ot; sim.openBottom = ob;
            sim.macCormackSmoke = ms; sim.macCormackVel = mv;
        }
        ImGui::Checkbox("Show velocity", &showSpeed);
        ImGui::Checkbox("Reflection", &useReflection);
        ImGui::Checkbox("MacCormack smoke", &sim.macCormackSmoke);
        ImGui::Checkbox("MacCormack vel", &sim.macCormackVel);
        ImGui::Checkbox("Debug print", &debugPrintOn);
        ImGui::SliderFloat("Buoyancy", &buoyancy, 0.0f, 10.0f);
        ImGui::SliderFloat("Source", &sourceStrength, 0.0f, 3.0f);
        ImGui::SliderInt("SOR iterations", &projectIterations, 10, 800);
        ImGui::SliderFloat("MC CFL limit", &sim.cflMc, 0.5f, 10.0f);
        ImGui::SeparatorText("Boundaries");
        ImGui::Checkbox("Open left", &sim.openLeft);
        ImGui::SameLine();
        ImGui::Checkbox("Open right", &sim.openRight);
        ImGui::Checkbox("Open top", &sim.openTop);
        ImGui::SameLine();
        ImGui::Checkbox("Open bottom", &sim.openBottom);
        ImGui::End();

        ImGui::Render();
        int w, h;
        glfwGetFramebufferSize(window, &w, &h);
        glViewport(0, 0, w, h);
        glClear(GL_COLOR_BUFFER_BIT);
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
        glfwSwapBuffers(window);
    }
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
    renderer.shutdown();
    glfwDestroyWindow(window);
    glfwTerminate();
}
