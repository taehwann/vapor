#pragma once
#include "common/MacAlgorithm.hpp"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <vector>
#include <stdexcept>
#include "common/FluidDiagnostics.hpp"
// Compatibility implementation of ba20987. Preserve operation order for replay.
// Random inlet stirring is opt-in (stirStrength > 0); 0 keeps the historical replay identical.
struct ReflectionPoint2D { float x, y; };
static ReflectionPoint2D operator-(ReflectionPoint2D a, ReflectionPoint2D b) { return { a.x - b.x, a.y - b.y }; }
static ReflectionPoint2D operator+(ReflectionPoint2D a, ReflectionPoint2D b) { return { a.x + b.x, a.y + b.y }; }
static ReflectionPoint2D operator*(ReflectionPoint2D a, float s) { return { a.x * s, a.y * s }; }


struct ReflectionMacFluidSolver2D {
    bool macCormackSmoke = true;
    float sourceStrength = 1.f, buoyancy = .8f;
    float smokeDecay = .3007525f; // exp(-rate/60) rounds to historical .995f.
    int projectIterations = 60;
    ReflectionMacFluidSolver2D& parameters() noexcept { return *this; }
    float boxSize() const noexcept { return 1.f; }
    void setBoxSize(float size) {
        if (size != 1.f) throw std::invalid_argument("2D reflection requires box_size=1");
    }
    void resetState() {
        for (auto* field : {&smoke,&vx,&vy,&pressure,&divergence}) std::fill(field->begin(),field->end(),0.f);
        randomState = randomSeed;
    }
    void advance(float dt) {
        if (!std::isfinite(emitterKickMultiplier) || emitterKickMultiplier < 0.f)
            throw std::invalid_argument("Invalid emitter kick multiplier");
        if (!std::isfinite(dt) || dt<=0.f || !std::isfinite(buoyancy) ||
            !std::isfinite(smokeDecay) || smokeDecay<0.f || !std::isfinite(sourceStrength) || sourceStrength<0.f ||
            !std::isfinite(emitterRadius) || emitterRadius<=0.f || emitterRadius>1.f ||
            !std::isfinite(stirStrength) || stirStrength<0.f ||
            !std::isfinite(sorOmega) || sorOmega<=0.f || sorOmega>=2.f || projectIterations<0)
            throw std::invalid_argument("Invalid 2D reflection step or parameters");
        applyBoundary();
        step(dt, buoyancy, std::exp(-smokeDecay/60.f), sourceStrength, projectIterations);
    }
    FluidDiagnostics diagnostics() const {
        FluidDiagnostics result;
        double energy=0.;
        for (const auto* field : {&vx,&vy}) for (float v : *field) {
            energy+=double(v)*v;
            result.maxVelocity=std::max(result.maxVelocity,std::abs(v));
        }
        result.kineticEnergy=float(.5*energy*h()*h());
        for (int y=0;y<n;++y) for(int x=0;x<n;++x)
            result.maxDivergence=std::max(result.maxDivergence,std::abs(
                (vx[idX(x+1,y)]-vx[idX(x,y)]+vy[idY(x,y+1)]-vy[idY(x,y)])/h()));
        return result;
    }
    int n = 256;
    float emitterRadius = 0.08f;
    float emitterFalloff = 3.5f;
    float emitterCenterX = 0.5f, emitterCenterY = 0.2f;
    float emitterSmokeScale = 3.0f;
    float emitterKickMultiplier = 1.f;
    float emitterKickBase = 0.6f;
    float emitterKickScale = 1.0f;
    int emitterHeight = 4;
    float sorOmega = 1.9f;
    float buoyancySplit = 0.5f;
    float buoyancyExp = 0.25f;
    float stirStrength = 0.f; // Random inlet stirring velocity; 0 keeps the historical replay identical.
    std::uint32_t randomSeed = 0x9E3779B9u;
    std::uint32_t randomState = 0x9E3779B9u;
    static constexpr bool openLeft = false, openRight = false, openTop = false, openBottom = false;
    static constexpr float decayFramesPerSec = 60.0f;
    std::vector<float> smoke, vx, vy, pressure, divergence;

    // Ready-to-run historical plume; velocity transport is always semi-Lagrangian.
    explicit ReflectionMacFluidSolver2D(AdvectionMethod smokeMethod = AdvectionMethod::SemiLagrangian)
        : ReflectionMacFluidSolver2D(128) {
        if (smokeMethod != AdvectionMethod::SemiLagrangian && smokeMethod != AdvectionMethod::MacCormack)
            throw std::invalid_argument("Unknown smoke advection method");
        macCormackSmoke = smokeMethod == AdvectionMethod::MacCormack;
        emitterRadius = .04f;
        projectIterations = 120;
    }
    explicit ReflectionMacFluidSolver2D(int res) {
        if (res<2 || res>512) throw std::invalid_argument("2D reflection resolution must be 2..512");
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
    static float mix(float a, float b, float t) { return a + (b - a) * t; }

    float randomUnit() {
        randomState = randomState * 1664525u + 1013904223u;
        return (randomState >> 8) * (1.0f / 16777216.0f);
    }

    float sampleCell(const std::vector<float>& f, ReflectionPoint2D p) const {
        const float gx = std::clamp(p.x / h() - 0.5f, 0.f, float(n - 1));
        const float gy = std::clamp(p.y / h() - 0.5f, 0.f, float(n - 1));
        const int x0 = int(gx), y0 = int(gy);
        const int x1 = std::min(x0 + 1, n - 1), y1 = std::min(y0 + 1, n - 1);
        const float tx = gx - x0, ty = gy - y0;
        return mix(mix(f[idC(x0, y0)], f[idC(x1, y0)], tx),
                   mix(f[idC(x0, y1)], f[idC(x1, y1)], tx), ty);
    }

    float sampleVxField(const std::vector<float>& f, ReflectionPoint2D p) const {
        const float gx = std::clamp(p.x / h(), 0.f, float(n));
        const float gy = std::clamp(p.y / h() - 0.5f, 0.f, float(n - 1));
        const int x0 = std::min(int(gx), n - 1), y0 = int(gy);
        const int x1 = x0 + 1, y1 = std::min(y0 + 1, n - 1);
        const float tx = gx - x0, ty = gy - y0;
        return mix(mix(f[idX(x0, y0)], f[idX(x1, y0)], tx),
                   mix(f[idX(x0, y1)], f[idX(x1, y1)], tx), ty);
    }

    float sampleVx(ReflectionPoint2D p) const { return sampleVxField(vx, p); }

    float sampleVyField(const std::vector<float>& f, ReflectionPoint2D p) const {
        const float gx = std::clamp(p.x / h() - 0.5f, 0.f, float(n - 1));
        const float gy = std::clamp(p.y / h(), 0.f, float(n));
        const int x0 = int(gx), y0 = std::min(int(gy), n - 1);
        const int x1 = std::min(x0 + 1, n - 1), y1 = y0 + 1;
        const float tx = gx - x0, ty = gy - y0;
        return mix(mix(f[idY(x0, y0)], f[idY(x1, y0)], tx),
                   mix(f[idY(x0, y1)], f[idY(x1, y1)], tx), ty);
    }

    float sampleVy(ReflectionPoint2D p) const { return sampleVyField(vy, p); }

    ReflectionPoint2D velocity(ReflectionPoint2D p) const { return { sampleVx(p), sampleVy(p) }; }
    ReflectionPoint2D velocity(ReflectionPoint2D p, const std::vector<float>& vxSrc, const std::vector<float>& vySrc) const {
        return { sampleVxField(vxSrc, p), sampleVyField(vySrc, p) };
    }

    float semiLagrangian(const std::vector<float>& source, ReflectionPoint2D p, float dt) const {
        ReflectionPoint2D q = p - velocity(p) * dt;
        const float pad = 0.5f * h();
        q.x = std::clamp(q.x, pad, 1.0f - pad);
        q.y = std::clamp(q.y, pad, 1.0f - pad);
        return sampleCell(source, q);
    }

    void advectScalar(std::vector<float>& field, float dt) {
        std::vector<float> fwd(n * n);
        const float dx = h();
        for (int y = 0; y < n; ++y)
            for (int x = 0; x < n; ++x)
                fwd[idC(x, y)] = semiLagrangian(field, { (x + .5f) * dx, (y + .5f) * dx }, dt);

        if (!macCormackSmoke) {
            field.swap(fwd);
            return;
        }

        std::vector<float> back(n * n);
        for (int y = 0; y < n; ++y)
            for (int x = 0; x < n; ++x) {
                ReflectionPoint2D p{ (x + .5f) * dx, (y + .5f) * dx };
                ReflectionPoint2D q = p + velocity(p) * dt;
                const float pad = 0.5f * h();
                q.x = std::clamp(q.x, pad, 1.0f - pad);
                q.y = std::clamp(q.y, pad, 1.0f - pad);
                back[idC(x, y)] = sampleCell(fwd, q);
            }

        for (int y = 0; y < n; ++y)
            for (int x = 0; x < n; ++x) {
                const int i = idC(x, y);
                float corrected = fwd[i] + 0.5f * (field[i] - back[i]);
                float fMin = std::min(field[i], fwd[i]), fMax = std::max(field[i], fwd[i]);
                for (int dy = -1; dy <= 1; ++dy)
                    for (int dx2 = -1; dx2 <= 1; ++dx2) {
                        const int a = x + dx2, b = y + dy;
                        if (inside(a, b)) {
                            fMin = std::min({ fMin, fwd[idC(a, b)], field[idC(a, b)] });
                            fMax = std::max({ fMax, fwd[idC(a, b)], field[idC(a, b)] });
                        }
                    }
                field[i] = std::clamp(corrected, fMin, fMax);
            }
    }

    void advectVx(float dt, const std::vector<float>& vxSrc, const std::vector<float>& vySrc) {
        std::vector<float> result(vx.size());
        const float dx = h();
        for (int y = 0; y < n; ++y)
            for (int x = 0; x <= n; ++x) {
                const ReflectionPoint2D face{ x * dx, (y + .5f) * dx };
                ReflectionPoint2D q = face - velocity(face, vxSrc, vySrc) * dt;
                q.x = std::clamp(q.x, 0.f, 1.0f);
                q.y = std::clamp(q.y, 0.5f * h(), 1.0f - 0.5f * h());
                result[idX(x, y)] = sampleVxField(vxSrc, q);
            }
        vx.swap(result);
    }

    void advectVy(float dt, const std::vector<float>& vxSrc, const std::vector<float>& vySrc) {
        std::vector<float> result(vy.size());
        const float dx = h();
        for (int y = 0; y <= n; ++y)
            for (int x = 0; x < n; ++x) {
                const ReflectionPoint2D face{ (x + .5f) * dx, y * dx };
                ReflectionPoint2D q = face - velocity(face, vxSrc, vySrc) * dt;
                q.x = std::clamp(q.x, 0.5f * h(), 1.0f - 0.5f * h());
                q.y = std::clamp(q.y, 0.f, 1.0f);
                result[idY(x, y)] = sampleVyField(vySrc, q);
            }
        vy.swap(result);
    }

    // Free slip: zero normal face velocity; tangential samples stay untouched.
    // Clamped sampling extends tangential velocity constantly beyond the wall.
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
        const int emitH = std::min(int(emitterCenterY * n + emitterRadius * n + 2), n);
        for (int y = std::max(0, int(emitterCenterY * n - emitterRadius * n - 1)); y < emitH; ++y)
            for (int x = 0; x < n; ++x) {
                ReflectionPoint2D p{ (x + .5f) * dx, (y + .5f) * dx };
                const float dx2 = (p.x - emitterCenterX) * (p.x - emitterCenterX);
                const float dy2 = (p.y - emitterCenterY) * (p.y - emitterCenterY);
                const float q = (dx2 + dy2) / r2;
                if (q > 1.0f) continue;
                const float w = std::exp(-emitterFalloff * q);
                const int i = idC(x, y);
                smoke[i] = std::max(smoke[i], std::min(1.0f, strength * emitterSmokeScale * w));
                if (emitterKickMultiplier > 0.f) vy[idY(x, y)] = std::max(vy[idY(x, y)], emitterKickMultiplier * (emitterKickBase + emitterKickScale * strength) * w);
                if (stirStrength > 0.f) {
                    vx[idX(x, y)] += (randomUnit() - 0.5f) * stirStrength * w;
                    vx[idX(x + 1, y)] += (randomUnit() - 0.5f) * stirStrength * w;
                }
            }
    }

    void step(float dt, float buoyancy, float dissipation, float sourceStrength, int projectIterations) {
        const float halfDt = 0.5f * dt;
        const int halfIters = std::max(1, projectIterations / 2);
        const float halfDecay = std::pow(dissipation, halfDt * decayFramesPerSec);
        const float dxInv = 1.0f / h(), dx = h();

        // ---- Stage 1: forward half-step of smoke using u₀ ----
        advectScalar(smoke, halfDt);

        for (int i = 0; i < n * n; ++i) smoke[i] *= halfDecay;

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

        { std::vector<float> vxOld = vx, vyOld = vy; advectVx(halfDt, vxOld, vyOld); }
        std::vector<float> vxTilde = vx;

        vx = vx0; vy = vy0;
        { std::vector<float> vxOld = vx, vyOld = vy; advectVy(halfDt, vxOld, vyOld); }
        std::vector<float> vyTilde = vy;

        vx = vxTilde;  // now vx,vy = ũ₁/₂ (forward-advected)

        for (int y = 0; y < n; ++y)
            for (int x = 0; x < n; ++x) {
                const float f = buoyancy * std::pow(smoke[idC(x, y)], buoyancyExp) * halfDt;
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

        for (int i = 0; i < n * n; ++i) smoke[i] *= halfDecay;

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

        // Advect reflected vx using midpoint velocity
        {
            std::vector<float> result(vx.size());
            for (int y = 0; y < n; ++y)
                for (int x = 0; x <= n; ++x) {
                    ReflectionPoint2D face{ x * dx, (y + .5f) * dx };
                    ReflectionPoint2D q = face - velocity(face, vxMid, vyMid) * halfDt;
                    q.x = std::clamp(q.x, 0.f, 1.f);
                    q.y = std::clamp(q.y, 0.5f * dx, 1.f - 0.5f * dx);
                    result[idX(x, y)] = sampleVxField(vxHat, q);
                }
            vx.swap(result);
        }
        // Advect reflected vy using midpoint velocity
        {
            std::vector<float> result(vy.size());
            for (int y = 0; y <= n; ++y)
                for (int x = 0; x < n; ++x) {
                    ReflectionPoint2D face{ (x + .5f) * dx, y * dx };
                    ReflectionPoint2D q = face - velocity(face, vxMid, vyMid) * halfDt;
                    q.x = std::clamp(q.x, 0.5f * dx, 1.f - 0.5f * dx);
                    q.y = std::clamp(q.y, 0.f, 1.f);
                    result[idY(x, y)] = sampleVyField(vyHat, q);
                }
            vy.swap(result);
        }

        project(halfDt, halfIters);
    }
};
