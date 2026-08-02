#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
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
    int emitterHeight = 4;
    float stirStrength = 0.6f;
    float sorOmega = 1.9f;
    float buoyancySplit = 0.5f;
    float buoyancyExp = 0.25f;
    bool openLeft = true, openRight = true, openTop = true, openBottom = false;
    static constexpr float decayFramesPerSec = 60.0f;
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
        Vec2 q = p - velocity(p) * dt;
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

        std::vector<float> back(n * n);
        for (int y = 0; y < n; ++y)
            for (int x = 0; x < n; ++x) {
                Vec2 p{ (x + .5f) * dx, (y + .5f) * dx };
                Vec2 q = p + velocity(p) * dt;
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
                const Vec2 face{ x * dx, (y + .5f) * dx };
                Vec2 q = face - velocity(face, vxSrc, vySrc) * dt;
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
                const Vec2 face{ (x + .5f) * dx, y * dx };
                Vec2 q = face - velocity(face, vxSrc, vySrc) * dt;
                q.x = std::clamp(q.x, 0.5f * h(), 1.0f - 0.5f * h());
                q.y = std::clamp(q.y, 0.f, 1.0f);
                result[idY(x, y)] = sampleVyField(vySrc, q);
            }
        vy.swap(result);
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
                vy[idY(x, y)] = std::max(vy[idY(x, y)], (emitterKickBase + emitterKickScale * strength) * w);
                vx[idX(x, y)] += (frand() - 0.5f) * stirStrength * w;
                vx[idX(x + 1, y)] += (frand() - 0.5f) * stirStrength * w;
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
                    Vec2 face{ x * dx, (y + .5f) * dx };
                    Vec2 q = face - velocity(face, vxMid, vyMid) * halfDt;
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
                    Vec2 face{ (x + .5f) * dx, y * dx };
                    Vec2 q = face - velocity(face, vxMid, vyMid) * halfDt;
                    q.x = std::clamp(q.x, 0.5f * dx, 1.f - 0.5f * dx);
                    q.y = std::clamp(q.y, 0.f, 1.f);
                    result[idY(x, y)] = sampleVyField(vyHat, q);
                }
            vy.swap(result);
        }

        project(halfDt, halfIters);
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

    bool paused = false, showSpeed = false;
    float buoyancy = 0.8f, dissipation = 0.995f, sourceStrength = 1.0f;
    int projectIterations = 60;
    const float maxDt = 0.033f;

    while (!glfwWindowShouldClose(window)) {
        glfwPollEvents();
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        ImGui::SetNextWindowPos({ 0, 0 });
        ImGui::SetNextWindowSize(ImGui::GetIO().DisplaySize);
        ImGui::Begin("Sim", nullptr, ImGuiWindowFlags_NoDecoration | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoBringToFrontOnFocus);

        if (!paused)
            sim.step(std::min(ImGui::GetIO().DeltaTime, maxDt), buoyancy, dissipation, sourceStrength, projectIterations);

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
            sim = SmokeSim2D(sim.n);
            sim.openLeft = ol; sim.openRight = orr; sim.openTop = ot; sim.openBottom = ob;
        }
        ImGui::Checkbox("Show velocity", &showSpeed);
        ImGui::SliderFloat("Buoyancy", &buoyancy, 0.0f, 10.0f);
        ImGui::SliderFloat("Buoyancy exp", &sim.buoyancyExp, 0.1f, 1.0f, "%.2f");
        ImGui::SliderFloat("Dissipation", &dissipation, 0.970f, 1.000f, "%.3f");
        ImGui::SliderFloat("Source", &sourceStrength, 0.0f, 3.0f);
        ImGui::SliderInt("SOR iterations", &projectIterations, 10, 200);
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
