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

static GLuint compileComputeShader(const char* src) {
    GLuint sh = glCreateShader(GL_COMPUTE_SHADER);
    glShaderSource(sh, 1, &src, nullptr);
    glCompileShader(sh);
    int ok = 0; glGetShaderiv(sh, GL_COMPILE_STATUS, &ok);
    if (!ok) { char log[1024]; glGetShaderInfoLog(sh, 1024, nullptr, log); std::fprintf(stderr, "Compute compile: %s\n", log); }
    GLuint prog = glCreateProgram();
    glAttachShader(prog, sh);
    glLinkProgram(prog);
    int linkOk = 0; glGetProgramiv(prog, GL_LINK_STATUS, &linkOk);
    if (!linkOk) { char log[1024]; glGetProgramInfoLog(prog, 1024, nullptr, log); std::fprintf(stderr, "Compute link: %s\n", log); }
    glDeleteShader(sh);
    return prog;
}

struct GpuSORSolver {
    GLuint progDiv = 0, progRBGS = 0, progCorrect = 0;
    GLuint progAdvectSL = 0, progAdvectMC = 0, progCopy = 0;
    GLuint ssboP = 0, ssboDiv = 0, ssboVx = 0, ssboVy = 0, ssboVz = 0;
    GLuint ssboSmoke = 0, ssboFwd = 0, ssboAdvectSrc = 0, ssboAdvectOut = 0;
    int n = 0, vxSz = 0, vySz = 0, vzSz = 0;
    int maxSz = 0;
    bool enabled = false, advectionEnabled = false;

    void init(int res, int vxSize, int vySize, int vzSize) {
        n = res; vxSz = vxSize; vySz = vySize; vzSz = vzSize;
        if (!glDispatchCompute) { enabled = false; return; }

        const char* divSrc = R"glsl(#version 430 core
layout(local_size_x = 8, local_size_y = 8, local_size_z = 4) in;
layout(std430, binding = 0) buffer VxBuf { float vx[]; };
layout(std430, binding = 1) buffer VyBuf { float vy[]; };
layout(std430, binding = 2) buffer VzBuf { float vz[]; };
layout(std430, binding = 3) buffer DivBuf { float div[]; };
uniform int u_n, u_n1; uniform float u_hInv;
void main() {
    int x = int(gl_GlobalInvocationID.x), y = int(gl_GlobalInvocationID.y), z = int(gl_GlobalInvocationID.z);
    if (x >= u_n || y >= u_n || z >= u_n) return;
    int idxC = x + u_n * (y + u_n * z);
    int idxX0 = x + u_n1 * (y + u_n * z), idxX1 = (x+1) + u_n1 * (y + u_n * z);
    int idxY0 = x + u_n * (y + u_n1 * z), idxY1 = x + u_n * ((y+1) + u_n1 * z);
    int idxZ0 = x + u_n * (y + u_n * z), idxZ1 = x + u_n * (y + u_n * (z+1));
    div[idxC] = (vx[idxX1] - vx[idxX0] + vy[idxY1] - vy[idxY0] + vz[idxZ1] - vz[idxZ0]) * u_hInv;
}
)glsl";

        const char* rbgsSrc = R"glsl(#version 430 core
layout(local_size_x = 4, local_size_y = 4, local_size_z = 4) in;
layout(std430, binding = 0) buffer PBuf { float p[]; };
layout(std430, binding = 1) buffer DBuf { float div[]; };
uniform int u_n, u_parity; uniform float u_omega, u_h2_div_dt;
void main() {
    int x = int(gl_GlobalInvocationID.x), y = int(gl_GlobalInvocationID.y), z = int(gl_GlobalInvocationID.z);
    if (x >= u_n || y >= u_n || z >= u_n) return;
    if (((x + y + z) & 1) != u_parity) return;
    const ivec3 nb[6] = { ivec3(-1,0,0), ivec3(1,0,0), ivec3(0,-1,0), ivec3(0,1,0), ivec3(0,0,-1), ivec3(0,0,1) };
    int i = x + u_n * (y + u_n * z);
    float sum = 0.0; int terms = 0;
    for (int k = 0; k < 6; ++k) {
        int nx = x + nb[k].x, ny = y + nb[k].y, nz = z + nb[k].z;
        if (nx < 0 || nx >= u_n || ny < 0 || ny >= u_n || nz < 0 || nz >= u_n) { sum += p[i]; ++terms; }
        else { sum += p[nx + u_n * (ny + u_n * nz)]; ++terms; }
    }
    p[i] = (1.0 - u_omega) * p[i] + u_omega * (sum - div[i] * u_h2_div_dt) / float(max(terms, 1));
}
)glsl";

        const char* correctSrc = R"glsl(#version 430 core
layout(local_size_x = 4, local_size_y = 4, local_size_z = 4) in;
layout(std430, binding = 0) buffer VxBuf { float vx[]; };
layout(std430, binding = 1) buffer VyBuf { float vy[]; };
layout(std430, binding = 2) buffer VzBuf { float vz[]; };
layout(std430, binding = 3) buffer PBuf { float p[]; };
uniform int u_n, u_mode, u_n1; uniform float u_dt, u_hInv;
int idC(int x, int y, int z) { return x + u_n * (y + u_n * z); }
int idX(int x, int y, int z) { return x + u_n1 * (y + u_n * z); }
int idY(int x, int y, int z) { return x + u_n * (y + u_n1 * z); }
int idZ(int x, int y, int z) { return x + u_n * (y + u_n * z); }
void main() {
    int x = int(gl_GlobalInvocationID.x), y = int(gl_GlobalInvocationID.y), z = int(gl_GlobalInvocationID.z);
    if (u_mode == 0) { if (x < 1 || x >= u_n || y >= u_n || z >= u_n) return; vx[idX(x,y,z)] -= u_dt * (p[idC(x,y,z)] - p[idC(x-1,y,z)]) * u_hInv; }
    else if (u_mode == 1) { if (x >= u_n || y < 1 || y >= u_n || z >= u_n) return; vy[idY(x,y,z)] -= u_dt * (p[idC(x,y,z)] - p[idC(x,y-1,z)]) * u_hInv; }
    else if (u_mode == 2) { if (x >= u_n || y >= u_n || z < 1 || z >= u_n) return; vz[idZ(x,y,z)] -= u_dt * (p[idC(x,y,z)] - p[idC(x,y,z-1)]) * u_hInv; }
}
)glsl";

        progDiv = compileComputeShader(divSrc);
        progRBGS = compileComputeShader(rbgsSrc);
        progCorrect = compileComputeShader(correctSrc);

        const char* advectSlSrc = R"glsl(#version 430 core
layout(local_size_x = 8, local_size_y = 8, local_size_z = 4) in;
layout(std430, binding = 0) buffer VxBuf { float vx[]; };
layout(std430, binding = 1) buffer VyBuf { float vy[]; };
layout(std430, binding = 2) buffer VzBuf { float vz[]; };
layout(std430, binding = 3) buffer SrcBuf { float src[]; };
layout(std430, binding = 4) buffer FwdBuf { float fwd[]; };
uniform int u_n, u_mode; uniform float u_dx, u_dt;
int idC(int x, int y, int z) { return x + u_n * (y + u_n * z); }
int idX(int x, int y, int z) { return x + (u_n + 1) * (y + u_n * z); }
int idY(int x, int y, int z) { return x + u_n * (y + (u_n + 1) * z); }
int idZ(int x, int y, int z) { return x + u_n * (y + u_n * z); }
vec3 velocity(vec3 p) {
    float gx = clamp(p.x / u_dx, 0.0f, float(u_n));
    float gy = clamp(p.y / u_dx - 0.5f, 0.0f, float(u_n - 1));
    float gz = clamp(p.z / u_dx - 0.5f, 0.0f, float(u_n - 1));
    int x0 = min(int(gx), u_n - 1), y0 = int(gy), z0 = int(gz);
    int x1 = x0 + 1, y1 = min(y0 + 1, u_n - 1), z1 = min(z0 + 1, u_n - 1);
    float tx = gx - float(x0), ty = gy - float(y0), tz = gz - float(z0);
    float vxVal = vx[idX(x0,y0,z0)]*(1-tx)*(1-ty)*(1-tz) + vx[idX(x1,y0,z0)]*tx*(1-ty)*(1-tz) +
                  vx[idX(x0,y1,z0)]*(1-tx)*ty*(1-tz) + vx[idX(x1,y1,z0)]*tx*ty*(1-tz) +
                  vx[idX(x0,y0,z1)]*(1-tx)*(1-ty)*tz + vx[idX(x1,y0,z1)]*tx*(1-ty)*tz +
                  vx[idX(x0,y1,z1)]*(1-tx)*ty*tz + vx[idX(x1,y1,z1)]*tx*ty*tz;
    gx = clamp(p.x / u_dx - 0.5f, 0.0f, float(u_n - 1));
    gy = clamp(p.y / u_dx, 0.0f, float(u_n));
    x0 = int(gx); y0 = min(int(gy), u_n - 1); z0 = int(gz);
    x1 = min(x0 + 1, u_n - 1); y1 = y0 + 1; z1 = min(z0 + 1, u_n - 1);
    tx = gx - float(x0); ty = gy - float(y0);
    float vyVal = vy[idY(x0,y0,z0)]*(1-tx)*(1-ty)*(1-tz) + vy[idY(x1,y0,z0)]*tx*(1-ty)*(1-tz) +
                  vy[idY(x0,y1,z0)]*(1-tx)*ty*(1-tz) + vy[idY(x1,y1,z0)]*tx*ty*(1-tz) +
                  vy[idY(x0,y0,z1)]*(1-tx)*(1-ty)*tz + vy[idY(x1,y0,z1)]*tx*(1-ty)*tz +
                  vy[idY(x0,y1,z1)]*(1-tx)*ty*tz + vy[idY(x1,y1,z1)]*tx*ty*tz;
    gy = clamp(p.y / u_dx - 0.5f, 0.0f, float(u_n - 1));
    gz = clamp(p.z / u_dx, 0.0f, float(u_n));
    y0 = int(gy); z0 = min(int(gz), u_n - 1);
    y1 = min(y0 + 1, u_n - 1); z1 = z0 + 1;
    tx = gx - float(x0); ty = gy - float(y0); tz = gz - float(z0);
    float vzVal = vz[idZ(x0,y0,z0)]*(1-tx)*(1-ty)*(1-tz) + vz[idZ(x1,y0,z0)]*tx*(1-ty)*(1-tz) +
                  vz[idZ(x0,y1,z0)]*(1-tx)*ty*(1-tz) + vz[idZ(x1,y1,z0)]*tx*ty*(1-tz) +
                  vz[idZ(x0,y0,z1)]*(1-tx)*(1-ty)*tz + vz[idZ(x1,y0,z1)]*tx*(1-ty)*tz +
                  vz[idZ(x0,y1,z1)]*(1-tx)*ty*tz + vz[idZ(x1,y1,z1)]*tx*ty*tz;
    return vec3(vxVal, vyVal, vzVal);
}
float sampleSrc(vec3 p) {
    float dx = u_dx; int n = u_n;
    if (u_mode == 0) {
        float gx = clamp(p.x/dx-0.5f, 0.0f, float(n-1));
        float gy = clamp(p.y/dx-0.5f, 0.0f, float(n-1));
        float gz = clamp(p.z/dx-0.5f, 0.0f, float(n-1));
        int x0=int(gx), y0=int(gy), z0=int(gz);
        int x1=min(x0+1,n-1), y1=min(y0+1,n-1), z1=min(z0+1,n-1);
        float tx=gx-float(x0), ty=gy-float(y0), tz=gz-float(z0);
        return src[idC(x0,y0,z0)]*(1-tx)*(1-ty)*(1-tz)+src[idC(x1,y0,z0)]*tx*(1-ty)*(1-tz)+
               src[idC(x0,y1,z0)]*(1-tx)*ty*(1-tz)+src[idC(x1,y1,z0)]*tx*ty*(1-tz)+
               src[idC(x0,y0,z1)]*(1-tx)*(1-ty)*tz+src[idC(x1,y0,z1)]*tx*(1-ty)*tz+
               src[idC(x0,y1,z1)]*(1-tx)*ty*tz+src[idC(x1,y1,z1)]*tx*ty*tz;
    } else if (u_mode == 1) {
        float gx = clamp(p.x/dx, 0.0f, float(n));
        float gy = clamp(p.y/dx-0.5f, 0.0f, float(n-1));
        float gz = clamp(p.z/dx-0.5f, 0.0f, float(n-1));
        int x0=min(int(gx),n-1), y0=int(gy), z0=int(gz);
        int x1=x0+1, y1=min(y0+1,n-1), z1=min(z0+1,n-1);
        float tx=gx-float(x0), ty=gy-float(y0), tz=gz-float(z0);
        return src[idX(x0,y0,z0)]*(1-tx)*(1-ty)*(1-tz)+src[idX(x1,y0,z0)]*tx*(1-ty)*(1-tz)+
               src[idX(x0,y1,z0)]*(1-tx)*ty*(1-tz)+src[idX(x1,y1,z0)]*tx*ty*(1-tz)+
               src[idX(x0,y0,z1)]*(1-tx)*(1-ty)*tz+src[idX(x1,y0,z1)]*tx*(1-ty)*tz+
               src[idX(x0,y1,z1)]*(1-tx)*ty*tz+src[idX(x1,y1,z1)]*tx*ty*tz;
    } else if (u_mode == 2) {
        float gx = clamp(p.x/dx-0.5f, 0.0f, float(n-1));
        float gy = clamp(p.y/dx, 0.0f, float(n));
        float gz = clamp(p.z/dx-0.5f, 0.0f, float(n-1));
        int x0=int(gx), y0=min(int(gy),n-1), z0=int(gz);
        int x1=min(x0+1,n-1), y1=y0+1, z1=min(z0+1,n-1);
        float tx=gx-float(x0), ty=gy-float(y0), tz=gz-float(z0);
        return src[idY(x0,y0,z0)]*(1-tx)*(1-ty)*(1-tz)+src[idY(x1,y0,z0)]*tx*(1-ty)*(1-tz)+
               src[idY(x0,y1,z0)]*(1-tx)*ty*(1-tz)+src[idY(x1,y1,z0)]*tx*ty*(1-tz)+
               src[idY(x0,y0,z1)]*(1-tx)*(1-ty)*tz+src[idY(x1,y0,z1)]*tx*(1-ty)*tz+
               src[idY(x0,y1,z1)]*(1-tx)*ty*tz+src[idY(x1,y1,z1)]*tx*ty*tz;
    } else {
        float gx = clamp(p.x/dx-0.5f, 0.0f, float(n-1));
        float gy = clamp(p.y/dx-0.5f, 0.0f, float(n-1));
        float gz = clamp(p.z/dx, 0.0f, float(n));
        int x0=int(gx), y0=int(gy), z0=min(int(gz),n-1);
        int x1=min(x0+1,n-1), y1=min(y0+1,n-1), z1=z0+1;
        float tx=gx-float(x0), ty=gy-float(y0), tz=gz-float(z0);
        return src[idZ(x0,y0,z0)]*(1-tx)*(1-ty)*(1-tz)+src[idZ(x1,y0,z0)]*tx*(1-ty)*(1-tz)+
               src[idZ(x0,y1,z0)]*(1-tx)*ty*(1-tz)+src[idZ(x1,y1,z0)]*tx*ty*(1-tz)+
               src[idZ(x0,y0,z1)]*(1-tx)*(1-ty)*tz+src[idZ(x1,y0,z1)]*tx*(1-ty)*tz+
               src[idZ(x0,y1,z1)]*(1-tx)*ty*tz+src[idZ(x1,y1,z1)]*tx*ty*tz;
    }
}
int srcIdx(int i) { if (u_mode == 1) return i; if (u_mode == 2) return i; return i; }
vec3 facePos(int x, int y, int z) {
    float dx = u_dx;
    if (u_mode == 0) return vec3((x+0.5f)*dx, (y+0.5f)*dx, (z+0.5f)*dx);
    if (u_mode == 1) return vec3(float(x)*dx, (y+0.5f)*dx, (z+0.5f)*dx);
    if (u_mode == 2) return vec3((x+0.5f)*dx, float(y)*dx, (z+0.5f)*dx);
    return vec3((x+0.5f)*dx, (y+0.5f)*dx, float(z)*dx);
}
void main() {
    int x = int(gl_GlobalInvocationID.x), y = int(gl_GlobalInvocationID.y), z = int(gl_GlobalInvocationID.z);
    if (u_mode == 0 && (x >= u_n || y >= u_n || z >= u_n)) return;
    if (u_mode == 1 && (x > u_n || y >= u_n || z >= u_n)) return;
    if (u_mode == 2 && (x >= u_n || y > u_n || z >= u_n)) return;
    if (u_mode == 3 && (x >= u_n || y >= u_n || z > u_n)) return;
    vec3 face = facePos(x, y, z);
    vec3 v = velocity(face);
    int i = (u_mode == 0) ? idC(x,y,z) : (u_mode == 1) ? idX(x,y,z) : (u_mode == 2) ? idY(x,y,z) : idZ(x,y,z);
    fwd[i] = sampleSrc(face - v * u_dt);
}
)glsl";

        const char* advectMcSrc = R"glsl(#version 430 core
layout(local_size_x = 8, local_size_y = 8, local_size_z = 4) in;
layout(std430, binding = 0) buffer VxBuf { float vx[]; };
layout(std430, binding = 1) buffer VyBuf { float vy[]; };
layout(std430, binding = 2) buffer VzBuf { float vz[]; };
layout(std430, binding = 3) buffer SrcBuf { float src[]; };
layout(std430, binding = 4) buffer FwdBuf { float fwd[]; };
layout(std430, binding = 5) buffer OutBuf { float outBuf[]; };
uniform int u_n, u_mode; uniform float u_dx, u_dt, u_boxSize, u_cflMc;
int idC(int x, int y, int z) { return x + u_n * (y + u_n * z); }
int idX(int x, int y, int z) { return x + (u_n + 1) * (y + u_n * z); }
int idY(int x, int y, int z) { return x + u_n * (y + (u_n + 1) * z); }
int idZ(int x, int y, int z) { return x + u_n * (y + u_n * z); }
vec3 velocity(vec3 p) {
    float gx = clamp(p.x / u_dx, 0.0f, float(u_n));
    float gy = clamp(p.y / u_dx - 0.5f, 0.0f, float(u_n - 1));
    float gz = clamp(p.z / u_dx - 0.5f, 0.0f, float(u_n - 1));
    int x0 = min(int(gx), u_n - 1), y0 = int(gy), z0 = int(gz);
    int x1 = x0 + 1, y1 = min(y0 + 1, u_n - 1), z1 = min(z0 + 1, u_n - 1);
    float tx = gx - float(x0), ty = gy - float(y0), tz = gz - float(z0);
    float vxVal = vx[idX(x0,y0,z0)]*(1-tx)*(1-ty)*(1-tz) + vx[idX(x1,y0,z0)]*tx*(1-ty)*(1-tz) +
                  vx[idX(x0,y1,z0)]*(1-tx)*ty*(1-tz) + vx[idX(x1,y1,z0)]*tx*ty*(1-tz) +
                  vx[idX(x0,y0,z1)]*(1-tx)*(1-ty)*tz + vx[idX(x1,y0,z1)]*tx*(1-ty)*tz +
                  vx[idX(x0,y1,z1)]*(1-tx)*ty*tz + vx[idX(x1,y1,z1)]*tx*ty*tz;
    gx = clamp(p.x / u_dx - 0.5f, 0.0f, float(u_n - 1));
    gy = clamp(p.y / u_dx, 0.0f, float(u_n));
    x0 = int(gx); y0 = min(int(gy), u_n - 1); z0 = int(gz);
    x1 = min(x0 + 1, u_n - 1); y1 = y0 + 1; z1 = min(z0 + 1, u_n - 1);
    tx = gx - float(x0); ty = gy - float(y0);
    float vyVal = vy[idY(x0,y0,z0)]*(1-tx)*(1-ty)*(1-tz) + vy[idY(x1,y0,z0)]*tx*(1-ty)*(1-tz) +
                  vy[idY(x0,y1,z0)]*(1-tx)*ty*(1-tz) + vy[idY(x1,y1,z0)]*tx*ty*(1-tz) +
                  vy[idY(x0,y0,z1)]*(1-tx)*(1-ty)*tz + vy[idY(x1,y0,z1)]*tx*(1-ty)*tz +
                  vy[idY(x0,y1,z1)]*(1-tx)*ty*tz + vy[idY(x1,y1,z1)]*tx*ty*tz;
    gy = clamp(p.y / u_dx - 0.5f, 0.0f, float(u_n - 1));
    gz = clamp(p.z / u_dx, 0.0f, float(u_n));
    y0 = int(gy); z0 = min(int(gz), u_n - 1);
    y1 = min(y0 + 1, u_n - 1); z1 = z0 + 1;
    tx = gx - float(x0); ty = gy - float(y0); tz = gz - float(z0);
    float vzVal = vz[idZ(x0,y0,z0)]*(1-tx)*(1-ty)*(1-tz) + vz[idZ(x1,y0,z0)]*tx*(1-ty)*(1-tz) +
                  vz[idZ(x0,y1,z0)]*(1-tx)*ty*(1-tz) + vz[idZ(x1,y1,z0)]*tx*ty*(1-tz) +
                  vz[idZ(x0,y0,z1)]*(1-tx)*(1-ty)*tz + vz[idZ(x1,y0,z1)]*tx*(1-ty)*tz +
                  vz[idZ(x0,y1,z1)]*(1-tx)*ty*tz + vz[idZ(x1,y1,z1)]*tx*ty*tz;
    return vec3(vxVal, vyVal, vzVal);
}
float sampleSrc(vec3 p) {
    float dx = u_dx; int n = u_n;
    if (u_mode == 0) {
        float gx = clamp(p.x/dx-0.5f, 0.0f, float(n-1));
        float gy = clamp(p.y/dx-0.5f, 0.0f, float(n-1));
        float gz = clamp(p.z/dx-0.5f, 0.0f, float(n-1));
        int x0=int(gx), y0=int(gy), z0=int(gz);
        int x1=min(x0+1,n-1), y1=min(y0+1,n-1), z1=min(z0+1,n-1);
        float tx=gx-float(x0), ty=gy-float(y0), tz=gz-float(z0);
        return src[idC(x0,y0,z0)]*(1-tx)*(1-ty)*(1-tz)+src[idC(x1,y0,z0)]*tx*(1-ty)*(1-tz)+
               src[idC(x0,y1,z0)]*(1-tx)*ty*(1-tz)+src[idC(x1,y1,z0)]*tx*ty*(1-tz)+
               src[idC(x0,y0,z1)]*(1-tx)*(1-ty)*tz+src[idC(x1,y0,z1)]*tx*(1-ty)*tz+
               src[idC(x0,y1,z1)]*(1-tx)*ty*tz+src[idC(x1,y1,z1)]*tx*ty*tz;
    } else if (u_mode == 1) {
        float gx = clamp(p.x/dx, 0.0f, float(n));
        float gy = clamp(p.y/dx-0.5f, 0.0f, float(n-1));
        float gz = clamp(p.z/dx-0.5f, 0.0f, float(n-1));
        int x0=min(int(gx),n-1), y0=int(gy), z0=int(gz);
        int x1=x0+1, y1=min(y0+1,n-1), z1=min(z0+1,n-1);
        float tx=gx-float(x0), ty=gy-float(y0), tz=gz-float(z0);
        return src[idX(x0,y0,z0)]*(1-tx)*(1-ty)*(1-tz)+src[idX(x1,y0,z0)]*tx*(1-ty)*(1-tz)+
               src[idX(x0,y1,z0)]*(1-tx)*ty*(1-tz)+src[idX(x1,y1,z0)]*tx*ty*(1-tz)+
               src[idX(x0,y0,z1)]*(1-tx)*(1-ty)*tz+src[idX(x1,y0,z1)]*tx*(1-ty)*tz+
               src[idX(x0,y1,z1)]*(1-tx)*ty*tz+src[idX(x1,y1,z1)]*tx*ty*tz;
    } else if (u_mode == 2) {
        float gx = clamp(p.x/dx-0.5f, 0.0f, float(n-1));
        float gy = clamp(p.y/dx, 0.0f, float(n));
        float gz = clamp(p.z/dx-0.5f, 0.0f, float(n-1));
        int x0=int(gx), y0=min(int(gy),n-1), z0=int(gz);
        int x1=min(x0+1,n-1), y1=y0+1, z1=min(z0+1,n-1);
        float tx=gx-float(x0), ty=gy-float(y0), tz=gz-float(z0);
        return src[idY(x0,y0,z0)]*(1-tx)*(1-ty)*(1-tz)+src[idY(x1,y0,z0)]*tx*(1-ty)*(1-tz)+
               src[idY(x0,y1,z0)]*(1-tx)*ty*(1-tz)+src[idY(x1,y1,z0)]*tx*ty*(1-tz)+
               src[idY(x0,y0,z1)]*(1-tx)*(1-ty)*tz+src[idY(x1,y0,z1)]*tx*(1-ty)*tz+
               src[idY(x0,y1,z1)]*(1-tx)*ty*tz+src[idY(x1,y1,z1)]*tx*ty*tz;
    } else {
        float gx = clamp(p.x/dx-0.5f, 0.0f, float(n-1));
        float gy = clamp(p.y/dx-0.5f, 0.0f, float(n-1));
        float gz = clamp(p.z/dx, 0.0f, float(n));
        int x0=int(gx), y0=int(gy), z0=min(int(gz),n-1);
        int x1=min(x0+1,n-1), y1=min(y0+1,n-1), z1=z0+1;
        float tx=gx-float(x0), ty=gy-float(y0), tz=gz-float(z0);
        return src[idZ(x0,y0,z0)]*(1-tx)*(1-ty)*(1-tz)+src[idZ(x1,y0,z0)]*tx*(1-ty)*(1-tz)+
               src[idZ(x0,y1,z0)]*(1-tx)*ty*(1-tz)+src[idZ(x1,y1,z0)]*tx*ty*(1-tz)+
               src[idZ(x0,y0,z1)]*(1-tx)*(1-ty)*tz+src[idZ(x1,y0,z1)]*tx*(1-ty)*tz+
               src[idZ(x0,y1,z1)]*(1-tx)*ty*tz+src[idZ(x1,y1,z1)]*tx*ty*tz;
    }
}
float sampleFwdBuf(vec3 p) {
    float dx = u_dx; int n = u_n;
    if (u_mode == 0) {
        float gx = clamp(p.x/dx-0.5f, 0.0f, float(n-1));
        float gy = clamp(p.y/dx-0.5f, 0.0f, float(n-1));
        float gz = clamp(p.z/dx-0.5f, 0.0f, float(n-1));
        int x0=int(gx), y0=int(gy), z0=int(gz);
        int x1=min(x0+1,n-1), y1=min(y0+1,n-1), z1=min(z0+1,n-1);
        float tx=gx-float(x0), ty=gy-float(y0), tz=gz-float(z0);
        return fwd[idC(x0,y0,z0)]*(1-tx)*(1-ty)*(1-tz)+fwd[idC(x1,y0,z0)]*tx*(1-ty)*(1-tz)+
               fwd[idC(x0,y1,z0)]*(1-tx)*ty*(1-tz)+fwd[idC(x1,y1,z0)]*tx*ty*(1-tz)+
               fwd[idC(x0,y0,z1)]*(1-tx)*(1-ty)*tz+fwd[idC(x1,y0,z1)]*tx*(1-ty)*tz+
               fwd[idC(x0,y1,z1)]*(1-tx)*ty*tz+fwd[idC(x1,y1,z1)]*tx*ty*tz;
    } else if (u_mode == 1) {
        float gx = clamp(p.x/dx, 0.0f, float(n));
        float gy = clamp(p.y/dx-0.5f, 0.0f, float(n-1));
        float gz = clamp(p.z/dx-0.5f, 0.0f, float(n-1));
        int x0=min(int(gx),n-1), y0=int(gy), z0=int(gz);
        int x1=x0+1, y1=min(y0+1,n-1), z1=min(z0+1,n-1);
        float tx=gx-float(x0), ty=gy-float(y0), tz=gz-float(z0);
        return fwd[idX(x0,y0,z0)]*(1-tx)*(1-ty)*(1-tz)+fwd[idX(x1,y0,z0)]*tx*(1-ty)*(1-tz)+
               fwd[idX(x0,y1,z0)]*(1-tx)*ty*(1-tz)+fwd[idX(x1,y1,z0)]*tx*ty*(1-tz)+
               fwd[idX(x0,y0,z1)]*(1-tx)*(1-ty)*tz+fwd[idX(x1,y0,z1)]*tx*(1-ty)*tz+
               fwd[idX(x0,y1,z1)]*(1-tx)*ty*tz+fwd[idX(x1,y1,z1)]*tx*ty*tz;
    } else if (u_mode == 2) {
        float gx = clamp(p.x/dx-0.5f, 0.0f, float(n-1));
        float gy = clamp(p.y/dx, 0.0f, float(n));
        float gz = clamp(p.z/dx-0.5f, 0.0f, float(n-1));
        int x0=int(gx), y0=min(int(gy),n-1), z0=int(gz);
        int x1=min(x0+1,n-1), y1=y0+1, z1=min(z0+1,n-1);
        float tx=gx-float(x0), ty=gy-float(y0), tz=gz-float(z0);
        return fwd[idY(x0,y0,z0)]*(1-tx)*(1-ty)*(1-tz)+fwd[idY(x1,y0,z0)]*tx*(1-ty)*(1-tz)+
               fwd[idY(x0,y1,z0)]*(1-tx)*ty*(1-tz)+fwd[idY(x1,y1,z0)]*tx*ty*(1-tz)+
               fwd[idY(x0,y0,z1)]*(1-tx)*(1-ty)*tz+fwd[idY(x1,y0,z1)]*tx*(1-ty)*tz+
               fwd[idY(x0,y1,z1)]*(1-tx)*ty*tz+fwd[idY(x1,y1,z1)]*tx*ty*tz;
    } else {
        float gx = clamp(p.x/dx-0.5f, 0.0f, float(n-1));
        float gy = clamp(p.y/dx-0.5f, 0.0f, float(n-1));
        float gz = clamp(p.z/dx, 0.0f, float(n));
        int x0=int(gx), y0=int(gy), z0=min(int(gz),n-1);
        int x1=min(x0+1,n-1), y1=min(y0+1,n-1), z1=z0+1;
        float tx=gx-float(x0), ty=gy-float(y0), tz=gz-float(z0);
        return fwd[idZ(x0,y0,z0)]*(1-tx)*(1-ty)*(1-tz)+fwd[idZ(x1,y0,z0)]*tx*(1-ty)*(1-tz)+
               fwd[idZ(x0,y1,z0)]*(1-tx)*ty*(1-tz)+fwd[idZ(x1,y1,z0)]*tx*ty*(1-tz)+
               fwd[idZ(x0,y0,z1)]*(1-tx)*(1-ty)*tz+fwd[idZ(x1,y0,z1)]*tx*(1-ty)*tz+
               fwd[idZ(x0,y1,z1)]*(1-tx)*ty*tz+fwd[idZ(x1,y1,z1)]*tx*ty*tz;
    }
}
vec3 facePos(int x, int y, int z) {
    float dx = u_dx;
    if (u_mode == 0) return vec3((x+0.5f)*dx, (y+0.5f)*dx, (z+0.5f)*dx);
    if (u_mode == 1) return vec3(float(x)*dx, (y+0.5f)*dx, (z+0.5f)*dx);
    if (u_mode == 2) return vec3((x+0.5f)*dx, float(y)*dx, (z+0.5f)*dx);
    return vec3((x+0.5f)*dx, (y+0.5f)*dx, float(z)*dx);
}
int clampX(int v) { return (u_mode == 1) ? clamp(v, 0, u_n) : clamp(v, 0, u_n - 1); }
int clampY(int v) { return (u_mode == 2) ? clamp(v, 0, u_n) : clamp(v, 0, u_n - 1); }
int clampZ(int v) { return (u_mode == 3) ? clamp(v, 0, u_n) : clamp(v, 0, u_n - 1); }
void main() {
    int x = int(gl_GlobalInvocationID.x), y = int(gl_GlobalInvocationID.y), z = int(gl_GlobalInvocationID.z);
    if (u_mode == 0 && (x >= u_n || y >= u_n || z >= u_n)) return;
    if (u_mode == 1 && (x > u_n || y >= u_n || z >= u_n)) return;
    if (u_mode == 2 && (x >= u_n || y > u_n || z >= u_n)) return;
    if (u_mode == 3 && (x >= u_n || y >= u_n || z > u_n)) return;
    int i = (u_mode == 0) ? idC(x,y,z) : (u_mode == 1) ? idX(x,y,z) : (u_mode == 2) ? idY(x,y,z) : idZ(x,y,z);
    vec3 face = facePos(x, y, z);
    vec3 v = velocity(face);
    vec3 posFwd = face - v * u_dt;
    if (posFwd.x < 0.0f || posFwd.x > u_boxSize || posFwd.y < 0.0f || posFwd.y > u_boxSize || posFwd.z < 0.0f || posFwd.z > u_boxSize ||
        face.x + v.x * u_dt < 0.0f || face.x + v.x * u_dt > u_boxSize ||
        face.y + v.y * u_dt < 0.0f || face.y + v.y * u_dt > u_boxSize ||
        face.z + v.z * u_dt < 0.0f || face.z + v.z * u_dt > u_boxSize ||
        length(v) * u_dt / u_dx > u_cflMc) {
        outBuf[i] = fwd[i];
        return;
    }
    float backVal = sampleFwdBuf(face + v * u_dt);
    float corrected = fwd[i] + 0.5f * (src[i] - backVal);
    int gx, gy, gz;
    float dx = u_dx;
    if (u_mode == 0) {
        gx = clamp(int(posFwd.x/dx-0.5f), 0, u_n-1);
        gy = clamp(int(posFwd.y/dx-0.5f), 0, u_n-1);
        gz = clamp(int(posFwd.z/dx-0.5f), 0, u_n-1);
    } else if (u_mode == 1) {
        gx = clamp(int(posFwd.x/dx), 0, u_n);
        gy = clamp(int(posFwd.y/dx-0.5f), 0, u_n-1);
        gz = clamp(int(posFwd.z/dx-0.5f), 0, u_n-1);
    } else if (u_mode == 2) {
        gx = clamp(int(posFwd.x/dx-0.5f), 0, u_n-1);
        gy = clamp(int(posFwd.y/dx), 0, u_n);
        gz = clamp(int(posFwd.z/dx-0.5f), 0, u_n-1);
    } else {
        gx = clamp(int(posFwd.x/dx-0.5f), 0, u_n-1);
        gy = clamp(int(posFwd.y/dx-0.5f), 0, u_n-1);
        gz = clamp(int(posFwd.z/dx), 0, u_n);
    }
    float fMin = src[i], fMax = fMin;
    for (int dz = 0; dz <= 1; ++dz) {
        for (int dy = 0; dy <= 1; ++dy) {
            for (int dx2 = 0; dx2 <= 1; ++dx2) {
                int nx = clampX(gx + dx2), ny = clampY(gy + dy), nz = clampZ(gz + dz);
                int ni = (u_mode == 0) ? idC(nx,ny,nz) : (u_mode == 1) ? idX(nx,ny,nz) : (u_mode == 2) ? idY(nx,ny,nz) : idZ(nx,ny,nz);
                float val = src[ni];
                fMin = min(fMin, val);
                fMax = max(fMax, val);
            }
        }
    }
    outBuf[i] = clamp(corrected, fMin, fMax);
}
)glsl";

        const char* copySrc = R"glsl(#version 430 core
layout(local_size_x = 64, local_size_y = 1, local_size_z = 1) in;
layout(std430, binding = 0) buffer SrcBuf { float src[]; };
layout(std430, binding = 1) buffer DstBuf { float dst[]; };
uniform int u_count;
void main() {
    int i = int(gl_GlobalInvocationID.x);
    if (i >= u_count) return;
    dst[i] = src[i];
}
)glsl";

        progAdvectSL = compileComputeShader(advectSlSrc);
        progAdvectMC = compileComputeShader(advectMcSrc);
        progCopy = compileComputeShader(copySrc);

        auto newSSBO = [](GLuint& handle, int count) {
            glGenBuffers(1, &handle);
            glBindBuffer(GL_SHADER_STORAGE_BUFFER, handle);
            glBufferData(GL_SHADER_STORAGE_BUFFER, (GLsizeiptr)count * sizeof(float), nullptr, GL_DYNAMIC_DRAW);
            };

        newSSBO(ssboP, n * n * n);
        newSSBO(ssboDiv, n * n * n);
        newSSBO(ssboVx, vxSz);
        newSSBO(ssboVy, vySz);
        newSSBO(ssboVz, vzSz);

        maxSz = std::max({ vxSz, vySz, vzSz, n * n * n });
        newSSBO(ssboSmoke, n * n * n);
        newSSBO(ssboFwd, maxSz);
        newSSBO(ssboAdvectSrc, maxSz);
        newSSBO(ssboAdvectOut, maxSz);

        enabled = true;
    }

    void project(float dt, int iterations, float omega, float hInv, float h,
        const float* vxData, const float* vyData, const float* vzData,
        float* vxOut, float* vyOut, float* vzOut) {
        const int N = n * n * n;

        glBindBuffer(GL_SHADER_STORAGE_BUFFER, ssboVx);
        glBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, (GLsizeiptr)vxSz * sizeof(float), vxData);
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, ssboVy);
        glBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, (GLsizeiptr)vySz * sizeof(float), vyData);
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, ssboVz);
        glBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, (GLsizeiptr)vzSz * sizeof(float), vzData);

        std::vector<float> zeroP(N, 0.f);
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, ssboP);
        glBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, (GLsizeiptr)N * sizeof(float), zeroP.data());

        glUseProgram(progDiv);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, ssboVx);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, ssboVy);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, ssboVz);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 3, ssboDiv);
        glUniform1i(glGetUniformLocation(progDiv, "u_n"), n);
        glUniform1i(glGetUniformLocation(progDiv, "u_n1"), n + 1);
        glUniform1f(glGetUniformLocation(progDiv, "u_hInv"), hInv);
        glDispatchCompute((n + 7) / 8, (n + 7) / 8, (n + 3) / 4);
        glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);

        glUseProgram(progRBGS);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, ssboP);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, ssboDiv);
        glUniform1i(glGetUniformLocation(progRBGS, "u_n"), n);
        glUniform1f(glGetUniformLocation(progRBGS, "u_omega"), omega);
        glUniform1f(glGetUniformLocation(progRBGS, "u_h2_div_dt"), h * h / dt);
        int gX = (n + 3) / 4, gY = (n + 3) / 4, gZ = (n + 3) / 4;
        for (int iter = 0; iter < iterations; ++iter) {
            glUniform1i(glGetUniformLocation(progRBGS, "u_parity"), 0);
            glDispatchCompute(gX, gY, gZ);
            glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
            glUniform1i(glGetUniformLocation(progRBGS, "u_parity"), 1);
            glDispatchCompute(gX, gY, gZ);
            glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
        }

        glUseProgram(progCorrect);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, ssboVx);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, ssboVy);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, ssboVz);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 3, ssboP);
        glUniform1i(glGetUniformLocation(progCorrect, "u_n"), n);
        glUniform1i(glGetUniformLocation(progCorrect, "u_n1"), n + 1);
        glUniform1f(glGetUniformLocation(progCorrect, "u_dt"), dt);
        glUniform1f(glGetUniformLocation(progCorrect, "u_hInv"), hInv);

        glUniform1i(glGetUniformLocation(progCorrect, "u_mode"), 0);
        glDispatchCompute(gX, gY, gZ);
        glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
        glUniform1i(glGetUniformLocation(progCorrect, "u_mode"), 1);
        glDispatchCompute(gX, gY, gZ);
        glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
        glUniform1i(glGetUniformLocation(progCorrect, "u_mode"), 2);
        glDispatchCompute(gX, gY, gZ);
        glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);

        glBindBuffer(GL_SHADER_STORAGE_BUFFER, ssboVx);
        glGetBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, (GLsizeiptr)vxSz * sizeof(float), vxOut);
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, ssboVy);
        glGetBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, (GLsizeiptr)vySz * sizeof(float), vyOut);
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, ssboVz);
        glGetBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, (GLsizeiptr)vzSz * sizeof(float), vzOut);
    }

    void uploadVelocity(const float* vxData, const float* vyData, const float* vzData) {
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, ssboVx);
        glBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, (GLsizeiptr)vxSz * sizeof(float), vxData);
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, ssboVy);
        glBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, (GLsizeiptr)vySz * sizeof(float), vyData);
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, ssboVz);
        glBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, (GLsizeiptr)vzSz * sizeof(float), vzData);
    }

    void dispatchSL(int mode, int count, float dx, float dt) {
        glUseProgram(progAdvectSL);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, ssboVx);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, ssboVy);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, ssboVz);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 3, ssboAdvectSrc);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 4, ssboFwd);
        glUniform1i(glGetUniformLocation(progAdvectSL, "u_n"), n);
        glUniform1i(glGetUniformLocation(progAdvectSL, "u_mode"), mode);
        glUniform1f(glGetUniformLocation(progAdvectSL, "u_dx"), dx);
        glUniform1f(glGetUniformLocation(progAdvectSL, "u_dt"), dt);
        int gx = (n + 7) / 8, gy = (n + 7) / 8, gz = (n + 3) / 4;
        if (mode == 1) gx = (n + 1 + 7) / 8;
        else if (mode == 2) gy = (n + 1 + 7) / 8;
        else if (mode == 3) gz = (n + 1 + 3) / 4;
        glDispatchCompute(gx, gy, gz);
        glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
    }

    void dispatchMC(int mode, float dx, float dt, float boxSize, float cflMc, GLuint srcBo, GLuint outBo) {
        glUseProgram(progAdvectMC);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, ssboVx);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, ssboVy);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, ssboVz);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 3, srcBo);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 4, ssboFwd);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 5, outBo);
        glUniform1i(glGetUniformLocation(progAdvectMC, "u_n"), n);
        glUniform1i(glGetUniformLocation(progAdvectMC, "u_mode"), mode);
        glUniform1f(glGetUniformLocation(progAdvectMC, "u_dx"), dx);
        glUniform1f(glGetUniformLocation(progAdvectMC, "u_dt"), dt);
        glUniform1f(glGetUniformLocation(progAdvectMC, "u_boxSize"), boxSize);
        glUniform1f(glGetUniformLocation(progAdvectMC, "u_cflMc"), cflMc);
        int gx = (n + 7) / 8, gy = (n + 7) / 8, gz = (n + 3) / 4;
        if (mode == 1) gx = (n + 1 + 7) / 8;
        else if (mode == 2) gy = (n + 1 + 7) / 8;
        else if (mode == 3) gz = (n + 1 + 3) / 4;
        glDispatchCompute(gx, gy, gz);
        glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
    }

    void copySSBO(GLuint dstBo, GLuint srcBo, int count) {
        glUseProgram(progCopy);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, srcBo);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, dstBo);
        glUniform1i(glGetUniformLocation(progCopy, "u_count"), count);
        glDispatchCompute((count + 63) / 64, 1, 1);
        glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
    }

    void advectSmokeGPU(float dt, float* smokeData, const float* vxVel, const float* vyVel, const float* vzVel, float boxSize, float cflMc, bool useMC) {
        const float dx = boxSize / n;
        const int N = n * n * n;

        glBindBuffer(GL_SHADER_STORAGE_BUFFER, ssboSmoke);
        glBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, (GLsizeiptr)N * sizeof(float), smokeData);
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, ssboAdvectSrc);
        glBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, (GLsizeiptr)N * sizeof(float), smokeData);
        uploadVelocity(vxVel, vyVel, vzVel);

        dispatchSL(0, N, dx, dt);

        if (!useMC) {
            copySSBO(ssboSmoke, ssboFwd, N);
        } else {
            dispatchMC(0, dx, dt, boxSize, cflMc, ssboSmoke, ssboAdvectOut);
            copySSBO(ssboSmoke, ssboAdvectOut, N);
        }

        glBindBuffer(GL_SHADER_STORAGE_BUFFER, ssboSmoke);
        glGetBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, (GLsizeiptr)N * sizeof(float), smokeData);
    }

    void advectVxGPU(float dt, const float* vxSrc, int srcSz, const float* vxVel, const float* vyVel, const float* vzVel, float* vxOut, float boxSize, float cflMc, bool useMC) {
        const float dx = boxSize / n;
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, ssboAdvectSrc);
        glBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, (GLsizeiptr)srcSz * sizeof(float), vxSrc);
        uploadVelocity(vxVel, vyVel, vzVel);

        dispatchSL(1, srcSz, dx, dt);

        if (!useMC) {
            glBindBuffer(GL_SHADER_STORAGE_BUFFER, ssboAdvectOut);
            copySSBO(ssboAdvectOut, ssboFwd, srcSz);
        } else {
            dispatchMC(1, dx, dt, boxSize, cflMc, ssboAdvectSrc, ssboAdvectOut);
        }

        glBindBuffer(GL_SHADER_STORAGE_BUFFER, ssboAdvectOut);
        glGetBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, (GLsizeiptr)srcSz * sizeof(float), vxOut);
    }

    void advectVyGPU(float dt, const float* vySrc, int srcSz, const float* vxVel, const float* vyVel, const float* vzVel, float* vyOut, float boxSize, float cflMc, bool useMC) {
        const float dx = boxSize / n;
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, ssboAdvectSrc);
        glBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, (GLsizeiptr)srcSz * sizeof(float), vySrc);
        uploadVelocity(vxVel, vyVel, vzVel);

        dispatchSL(2, srcSz, dx, dt);

        if (!useMC) {
            glBindBuffer(GL_SHADER_STORAGE_BUFFER, ssboAdvectOut);
            copySSBO(ssboAdvectOut, ssboFwd, srcSz);
        } else {
            dispatchMC(2, dx, dt, boxSize, cflMc, ssboAdvectSrc, ssboAdvectOut);
        }

        glBindBuffer(GL_SHADER_STORAGE_BUFFER, ssboAdvectOut);
        glGetBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, (GLsizeiptr)srcSz * sizeof(float), vyOut);
    }

    void advectVzGPU(float dt, const float* vzSrc, int srcSz, const float* vxVel, const float* vyVel, const float* vzVel, float* vzOut, float boxSize, float cflMc, bool useMC) {
        const float dx = boxSize / n;
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, ssboAdvectSrc);
        glBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, (GLsizeiptr)srcSz * sizeof(float), vzSrc);
        uploadVelocity(vxVel, vyVel, vzVel);

        dispatchSL(3, srcSz, dx, dt);

        if (!useMC) {
            glBindBuffer(GL_SHADER_STORAGE_BUFFER, ssboAdvectOut);
            copySSBO(ssboAdvectOut, ssboFwd, srcSz);
        } else {
            dispatchMC(3, dx, dt, boxSize, cflMc, ssboAdvectSrc, ssboAdvectOut);
        }

        glBindBuffer(GL_SHADER_STORAGE_BUFFER, ssboAdvectOut);
        glGetBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, (GLsizeiptr)srcSz * sizeof(float), vzOut);
    }

    void shutdown() {
        if (progDiv) glDeleteProgram(progDiv);
        if (progRBGS) glDeleteProgram(progRBGS);
        if (progCorrect) glDeleteProgram(progCorrect);
        if (progAdvectSL) glDeleteProgram(progAdvectSL);
        if (progAdvectMC) glDeleteProgram(progAdvectMC);
        if (progCopy) glDeleteProgram(progCopy);
        if (ssboP) glDeleteBuffers(1, &ssboP);
        if (ssboDiv) glDeleteBuffers(1, &ssboDiv);
        if (ssboVx) glDeleteBuffers(1, &ssboVx);
        if (ssboVy) glDeleteBuffers(1, &ssboVy);
        if (ssboVz) glDeleteBuffers(1, &ssboVz);
        if (ssboSmoke) glDeleteBuffers(1, &ssboSmoke);
        if (ssboFwd) glDeleteBuffers(1, &ssboFwd);
        if (ssboAdvectSrc) glDeleteBuffers(1, &ssboAdvectSrc);
        if (ssboAdvectOut) glDeleteBuffers(1, &ssboAdvectOut);
    }
};

struct SmokeSim3D {
    int n = 64;
    // Spacious room defaults so smoke dissipates before hitting walls
    float boxSize = 4.5f;
    float emitterRadius = 0.035f;
    float emitterCenterX = 0.5f, emitterCenterY = 0.5f, emitterCenterZ = 0.1f;
    float emitterKickBase = 0.8f;
    float emitterKickScale = 0.8f;
    float stirStrength = 0.5f;
    float sorOmega = 1.95f;
    float buoyancySplit = 0.5f;
    float cflMc = 3.0f;
    bool macCormackSmoke = true, macCormackVel = true;

    // Dissipation tuned to fade out before reaching the top wall
    float smokeDecay = 0.06f;

    std::vector<float> smoke, vx, vy, vz, pressure, divergence;

    // Persistent scratch buffers to eliminate per-frame heap allocations
    std::vector<float> fwdScratch, backScratch;
    std::vector<float> vx0, vy0, vz0;
    std::vector<float> vxTilde, vyTilde, vzTilde;
    std::vector<float> vxHat, vyHat, vzHat;
    std::vector<float> vxDiv0, vyDiv0, vzDiv0;

    GpuSORSolver gpuSolver;

    explicit SmokeSim3D(int res) {
        n = res;
        allocateBuffers();
    }

    void allocateBuffers() {
        const int N = n * n * n;
        const int NX = (n + 1) * n * n;
        const int NY = n * (n + 1) * n;
        const int NZ = n * n * (n + 1);

        smoke.assign(N, 0.f);
        pressure.assign(N, 0.f);
        divergence.assign(N, 0.f);
        fwdScratch.assign(std::max({ NX, NY, NZ }), 0.f);
        backScratch.assign(std::max({ NX, NY, NZ }), 0.f);

        vx.assign(NX, 0.f);
        vy.assign(NY, 0.f);
        vz.assign(NZ, 0.f);

        vx0.assign(NX, 0.f);
        vy0.assign(NY, 0.f);
        vz0.assign(NZ, 0.f);

        vxTilde.assign(NX, 0.f);
        vyTilde.assign(NY, 0.f);
        vzTilde.assign(NZ, 0.f);

        vxHat.assign(NX, 0.f);
        vyHat.assign(NY, 0.f);
        vzHat.assign(NZ, 0.f);

        gpuSolver.init(n, (int)vx.size(), (int)vy.size(), (int)vz.size());
    }

    inline int idC(int x, int y, int z) const { return x + n * (y + n * z); }
    inline int idX(int x, int y, int z) const { return x + (n + 1) * (y + n * z); }
    inline int idY(int x, int y, int z) const { return x + n * (y + (n + 1) * z); }
    inline int idZ(int x, int y, int z) const { return x + n * (y + n * z); }
    inline float h() const { return boxSize / n; }
    inline bool inside(int x, int y, int z) const { return x >= 0 && y >= 0 && z >= 0 && x < n && y < n && z < n; }
    inline bool insideVx(int x, int y, int z) const { return x >= 0 && y >= 0 && z >= 0 && x <= n && y < n && z < n; }
    inline bool insideVy(int x, int y, int z) const { return x >= 0 && y >= 0 && z >= 0 && x < n && y <= n && z < n; }
    inline bool insideVz(int x, int y, int z) const { return x >= 0 && y >= 0 && z >= 0 && x < n && y < n && z <= n; }

    inline float sampleCell(const std::vector<float>& f, Vec3 p) const {
        float gx = std::clamp(p.x / h() - 0.5f, 0.f, float(n - 1));
        float gy = std::clamp(p.y / h() - 0.5f, 0.f, float(n - 1));
        float gz = std::clamp(p.z / h() - 0.5f, 0.f, float(n - 1));
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

    inline float sampleVxField(const std::vector<float>& f, Vec3 p) const {
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

    inline float sampleVyField(const std::vector<float>& f, Vec3 p) const {
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

    inline float sampleVzField(const std::vector<float>& f, Vec3 p) const {
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

    inline float sampleVx(Vec3 p) const { return sampleVxField(vx, p); }
    inline float sampleVy(Vec3 p) const { return sampleVyField(vy, p); }
    inline float sampleVz(Vec3 p) const { return sampleVzField(vz, p); }
    inline Vec3 velocity(Vec3 p) const { return { sampleVx(p), sampleVy(p), sampleVz(p) }; }
    inline Vec3 velocity(Vec3 p, const std::vector<float>& vxSrc, const std::vector<float>& vySrc, const std::vector<float>& vzSrc) const {
        return { sampleVxField(vxSrc, p), sampleVyField(vySrc, p), sampleVzField(vzSrc, p) };
    }

    void advectScalar(std::vector<float>& field, float dt) {
        const float dx = h();
        const int N = n * n * n;

#pragma omp parallel for
        for (int z = 0; z < n; ++z) {
            for (int y = 0; y < n; ++y) {
                for (int x = 0; x < n; ++x) {
                    Vec3 p{ (x + .5f) * dx, (y + .5f) * dx, (z + .5f) * dx };
                    fwdScratch[idC(x, y, z)] = sampleCell(field, p - velocity(p) * dt);
                }
            }
        }

        if (!macCormackSmoke) {
            std::copy(fwdScratch.begin(), fwdScratch.begin() + N, field.begin());
            return;
        }

#pragma omp parallel for
        for (int z = 0; z < n; ++z) {
            for (int y = 0; y < n; ++y) {
                for (int x = 0; x < n; ++x) {
                    Vec3 p{ (x + .5f) * dx, (y + .5f) * dx, (z + .5f) * dx };
                    backScratch[idC(x, y, z)] = sampleCell(fwdScratch, p + velocity(p) * dt);
                }
            }
        }

#pragma omp parallel for
        for (int z = 0; z < n; ++z) {
            for (int y = 0; y < n; ++y) {
                for (int x = 0; x < n; ++x) {
                    const int i = idC(x, y, z);
                    Vec3 p{ (x + .5f) * dx, (y + .5f) * dx, (z + .5f) * dx };
                    Vec3 v = velocity(p);
                    Vec3 posFwd = p - v * dt;
                    Vec3 posBwd = p + v * dt;
                    if (posFwd.x < 0 || posFwd.x > boxSize || posFwd.y < 0 || posFwd.y > boxSize || posFwd.z < 0 || posFwd.z > boxSize ||
                        posBwd.x < 0 || posBwd.x > boxSize || posBwd.y < 0 || posBwd.y > boxSize || posBwd.z < 0 || posBwd.z > boxSize) {
                        field[i] = fwdScratch[i];
                        continue;
                    }
                    float corrected = fwdScratch[i] + 0.5f * (field[i] - backScratch[i]);
                    int gx = std::clamp(int(posFwd.x / dx - 0.5f), 0, n - 1);
                    int gy = std::clamp(int(posFwd.y / dx - 0.5f), 0, n - 1);
                    int gz = std::clamp(int(posFwd.z / dx - 0.5f), 0, n - 1);
                    float fMin = field[idC(gx, gy, gz)], fMax = fMin;
                    for (int dz = 0; dz <= 1; ++dz)
                        for (int dy = 0; dy <= 1; ++dy)
                            for (int dx2 = 0; dx2 <= 1; ++dx2) {
                                float val = field[idC(std::clamp(gx + dx2, 0, n - 1),
                                    std::clamp(gy + dy, 0, n - 1),
                                    std::clamp(gz + dz, 0, n - 1))];
                                fMin = std::min(fMin, val);
                                fMax = std::max(fMax, val);
                            }
                    field[i] = std::clamp(corrected, fMin, fMax);
                }
            }
        }
    }

    void advectVx(float dt, const std::vector<float>& vxSrc,
        const std::vector<float>& vxVel, const std::vector<float>& vyVel, const std::vector<float>& vzVel,
        std::vector<float>& out, bool useMC = true) {
        const float dx = h();
        const int N = (int)out.size();

#pragma omp parallel for
        for (int z = 0; z < n; ++z) {
            for (int y = 0; y < n; ++y) {
                for (int x = 0; x <= n; ++x) {
                    Vec3 face{ x * dx, (y + .5f) * dx, (z + .5f) * dx };
                    fwdScratch[idX(x, y, z)] = sampleVxField(vxSrc, face - velocity(face, vxVel, vyVel, vzVel) * dt);
                }
            }
        }

        if (!useMC || !macCormackVel) {
            std::copy(fwdScratch.begin(), fwdScratch.begin() + N, out.begin());
            return;
        }

#pragma omp parallel for
        for (int z = 0; z < n; ++z) {
            for (int y = 0; y < n; ++y) {
                for (int x = 0; x <= n; ++x) {
                    Vec3 face{ x * dx, (y + .5f) * dx, (z + .5f) * dx };
                    backScratch[idX(x, y, z)] = sampleVxField(fwdScratch, face + velocity(face, vxVel, vyVel, vzVel) * dt);
                }
            }
        }

#pragma omp parallel for
        for (int z = 0; z < n; ++z) {
            for (int y = 0; y < n; ++y) {
                for (int x = 0; x <= n; ++x) {
                    const int i = idX(x, y, z);
                    Vec3 face{ x * dx, (y + .5f) * dx, (z + .5f) * dx };
                    Vec3 v = velocity(face, vxVel, vyVel, vzVel);
                    Vec3 posFwd = face - v * dt, posBwd = face + v * dt;
                    float cfl = std::sqrt(v.x * v.x + v.y * v.y + v.z * v.z) * dt / dx;
                    if (posFwd.x < 0 || posFwd.x > boxSize || posFwd.y < 0 || posFwd.y > boxSize || posFwd.z < 0 || posFwd.z > boxSize ||
                        posBwd.x < 0 || posBwd.x > boxSize || posBwd.y < 0 || posBwd.y > boxSize || posBwd.z < 0 || posBwd.z > boxSize ||
                        cfl > cflMc) {
                        out[i] = fwdScratch[i];
                        continue;
                    }
                    float corrected = fwdScratch[i] + 0.5f * (vxSrc[i] - backScratch[i]);
                    int gx = std::clamp(int(posFwd.x / dx), 0, n);
                    int gy = std::clamp(int(posFwd.y / dx - 0.5f), 0, n - 1);
                    int gz = std::clamp(int(posFwd.z / dx - 0.5f), 0, n - 1);
                    float fMin = vxSrc[idX(gx, gy, gz)], fMax = fMin;
                    for (int dz = 0; dz <= 1; ++dz)
                        for (int dy = 0; dy <= 1; ++dy)
                            for (int dx2 = 0; dx2 <= 1; ++dx2) {
                                float val = vxSrc[idX(std::clamp(gx + dx2, 0, n),
                                    std::clamp(gy + dy, 0, n - 1),
                                    std::clamp(gz + dz, 0, n - 1))];
                                fMin = std::min(fMin, val);
                                fMax = std::max(fMax, val);
                            }
                    out[i] = std::clamp(corrected, fMin, fMax);
                }
            }
        }
    }

    void advectVy(float dt, const std::vector<float>& vySrc,
        const std::vector<float>& vxVel, const std::vector<float>& vyVel, const std::vector<float>& vzVel,
        std::vector<float>& out, bool useMC = true) {
        const float dx = h();
        const int N = (int)out.size();

#pragma omp parallel for
        for (int z = 0; z < n; ++z) {
            for (int y = 0; y <= n; ++y) {
                for (int x = 0; x < n; ++x) {
                    Vec3 face{ (x + .5f) * dx, y * dx, (z + .5f) * dx };
                    fwdScratch[idY(x, y, z)] = sampleVyField(vySrc, face - velocity(face, vxVel, vyVel, vzVel) * dt);
                }
            }
        }

        if (!useMC || !macCormackVel) {
            std::copy(fwdScratch.begin(), fwdScratch.begin() + N, out.begin());
            return;
        }

#pragma omp parallel for
        for (int z = 0; z < n; ++z) {
            for (int y = 0; y <= n; ++y) {
                for (int x = 0; x < n; ++x) {
                    Vec3 face{ (x + .5f) * dx, y * dx, (z + .5f) * dx };
                    backScratch[idY(x, y, z)] = sampleVyField(fwdScratch, face + velocity(face, vxVel, vyVel, vzVel) * dt);
                }
            }
        }

#pragma omp parallel for
        for (int z = 0; z < n; ++z) {
            for (int y = 0; y <= n; ++y) {
                for (int x = 0; x < n; ++x) {
                    const int i = idY(x, y, z);
                    Vec3 face{ (x + .5f) * dx, y * dx, (z + .5f) * dx };
                    Vec3 v = velocity(face, vxVel, vyVel, vzVel);
                    Vec3 posFwd = face - v * dt, posBwd = face + v * dt;
                    float cfl = std::sqrt(v.x * v.x + v.y * v.y + v.z * v.z) * dt / dx;
                    if (posFwd.x < 0 || posFwd.x > boxSize || posFwd.y < 0 || posFwd.y > boxSize || posFwd.z < 0 || posFwd.z > boxSize ||
                        posBwd.x < 0 || posBwd.x > boxSize || posBwd.y < 0 || posBwd.y > boxSize || posBwd.z < 0 || posBwd.z > boxSize ||
                        cfl > cflMc) {
                        out[i] = fwdScratch[i];
                        continue;
                    }
                    float corrected = fwdScratch[i] + 0.5f * (vySrc[i] - backScratch[i]);
                    int gx = std::clamp(int(posFwd.x / dx - 0.5f), 0, n - 1);
                    int gy = std::clamp(int(posFwd.y / dx), 0, n);
                    int gz = std::clamp(int(posFwd.z / dx - 0.5f), 0, n - 1);
                    float fMin = vySrc[idY(gx, gy, gz)], fMax = fMin;
                    for (int dz = 0; dz <= 1; ++dz)
                        for (int dy = 0; dy <= 1; ++dy)
                            for (int dx2 = 0; dx2 <= 1; ++dx2) {
                                float val = vySrc[idY(std::clamp(gx + dx2, 0, n - 1),
                                    std::clamp(gy + dy, 0, n),
                                    std::clamp(gz + dz, 0, n - 1))];
                                fMin = std::min(fMin, val);
                                fMax = std::max(fMax, val);
                            }
                    out[i] = std::clamp(corrected, fMin, fMax);
                }
            }
        }
    }

    void advectVz(float dt, const std::vector<float>& vzSrc,
        const std::vector<float>& vxVel, const std::vector<float>& vyVel, const std::vector<float>& vzVel,
        std::vector<float>& out, bool useMC = true) {
        const float dx = h();
        const int N = (int)out.size();

#pragma omp parallel for
        for (int z = 0; z <= n; ++z) {
            for (int y = 0; y < n; ++y) {
                for (int x = 0; x < n; ++x) {
                    Vec3 face{ (x + .5f) * dx, (y + .5f) * dx, z * dx };
                    fwdScratch[idZ(x, y, z)] = sampleVzField(vzSrc, face - velocity(face, vxVel, vyVel, vzVel) * dt);
                }
            }
        }

        if (!useMC || !macCormackVel) {
            std::copy(fwdScratch.begin(), fwdScratch.begin() + N, out.begin());
            return;
        }

#pragma omp parallel for
        for (int z = 0; z <= n; ++z) {
            for (int y = 0; y < n; ++y) {
                for (int x = 0; x < n; ++x) {
                    Vec3 face{ (x + .5f) * dx, (y + .5f) * dx, z * dx };
                    backScratch[idZ(x, y, z)] = sampleVzField(fwdScratch, face + velocity(face, vxVel, vyVel, vzVel) * dt);
                }
            }
        }

#pragma omp parallel for
        for (int z = 0; z <= n; ++z) {
            for (int y = 0; y < n; ++y) {
                for (int x = 0; x < n; ++x) {
                    const int i = idZ(x, y, z);
                    Vec3 face{ (x + .5f) * dx, (y + .5f) * dx, z * dx };
                    Vec3 v = velocity(face, vxVel, vyVel, vzVel);
                    Vec3 posFwd = face - v * dt, posBwd = face + v * dt;
                    float cfl = std::sqrt(v.x * v.x + v.y * v.y + v.z * v.z) * dt / dx;
                    if (posFwd.x < 0 || posFwd.x > boxSize || posFwd.y < 0 || posFwd.y > boxSize || posFwd.z < 0 || posFwd.z > boxSize ||
                        posBwd.x < 0 || posBwd.x > boxSize || posBwd.y < 0 || posBwd.y > boxSize || posBwd.z < 0 || posBwd.z > boxSize ||
                        cfl > cflMc) {
                        out[i] = fwdScratch[i];
                        continue;
                    }
                    float corrected = fwdScratch[i] + 0.5f * (vzSrc[i] - backScratch[i]);
                    int gx = std::clamp(int(posFwd.x / dx - 0.5f), 0, n - 1);
                    int gy = std::clamp(int(posFwd.y / dx - 0.5f), 0, n - 1);
                    int gz = std::clamp(int(posFwd.z / dx), 0, n);
                    float fMin = vzSrc[idZ(gx, gy, gz)], fMax = fMin;
                    for (int dz = 0; dz <= 1; ++dz)
                        for (int dy = 0; dy <= 1; ++dy)
                            for (int dx2 = 0; dx2 <= 1; ++dx2) {
                                float val = vzSrc[idZ(std::clamp(gx + dx2, 0, n - 1),
                                    std::clamp(gy + dy, 0, n - 1),
                                    std::clamp(gz + dz, 0, n))];
                                fMin = std::min(fMin, val);
                                fMax = std::max(fMax, val);
                            }
                    out[i] = std::clamp(corrected, fMin, fMax);
                }
            }
        }
    }

    // 100% Closed solid wall boundaries: enforce zero normal velocity at all 6 outer faces
    void enforceVelocityBoundaries() {
        for (int z = 0; z < n; ++z) {
            for (int y = 0; y < n; ++y) {
                vx[idX(0, y, z)] = 0.f;
                vx[idX(n, y, z)] = 0.f;
            }
        }
        for (int z = 0; z < n; ++z) {
            for (int x = 0; x < n; ++x) {
                vy[idY(x, 0, z)] = 0.f;
                vy[idY(x, n, z)] = 0.f;
            }
        }
        for (int y = 0; y < n; ++y) {
            for (int x = 0; x < n; ++x) {
                vz[idZ(x, y, 0)] = 0.f;
                vz[idZ(x, y, n)] = 0.f;
            }
        }
    }

    void project(float dt, int iterations) {
        const float hInv = 1.0f / h(), h2 = h() * h();

        enforceVelocityBoundaries();

        if (gpuSolver.enabled) {
            gpuSolver.project(dt, iterations, sorOmega, hInv, h(),
                vx.data(), vy.data(), vz.data(),
                vx.data(), vy.data(), vz.data());
        } else {
#pragma omp parallel for
            for (int z = 0; z < n; ++z) {
                for (int y = 0; y < n; ++y) {
                    for (int x = 0; x < n; ++x) {
                        divergence[idC(x, y, z)] = (vx[idX(x + 1, y, z)] - vx[idX(x, y, z)] +
                            vy[idY(x, y + 1, z)] - vy[idY(x, y, z)] +
                            vz[idZ(x, y, z + 1)] - vz[idZ(x, y, z)]) * hInv;
                    }
                }
            }

            std::fill(pressure.begin(), pressure.end(), 0.f);
            const int dx[6] = { -1, 1, 0, 0, 0, 0 }, dy[6] = { 0, 0, -1, 1, 0, 0 }, dz[6] = { 0, 0, 0, 0, -1, 1 };
            const float omega = sorOmega;

            for (int iter = 0; iter < iterations; ++iter) {
                for (int parity = 0; parity < 2; ++parity) {
#pragma omp parallel for
                    for (int z = 0; z < n; ++z) {
                        for (int y = 0; y < n; ++y) {
                            int xStart = (y + z + parity) & 1;
                            for (int x = xStart; x < n; x += 2) {
                                const int i = idC(x, y, z);
                                float sum = 0.f;
                                int terms = 0;
                                for (int k = 0; k < 6; ++k) {
                                    int a = x + dx[k], b = y + dy[k], c = z + dz[k];
                                    if (!inside(a, b, c)) {
                                        sum += pressure[i];
                                        ++terms;
                                        continue;
                                    }
                                    ++terms;
                                    sum += pressure[idC(a, b, c)];
                                }
                                pressure[i] = (1.f - omega) * pressure[i] +
                                    omega * (sum - divergence[i] * h2 / dt) / std::max(terms, 1);
                            }
                        }
                    }
                }
            }

#pragma omp parallel for
            for (int z = 0; z < n; ++z) {
                for (int y = 0; y < n; ++y) {
                    for (int x = 1; x < n; ++x) {
                        vx[idX(x, y, z)] -= dt * (pressure[idC(x, y, z)] - pressure[idC(x - 1, y, z)]) * hInv;
                    }
                }
            }

#pragma omp parallel for
            for (int z = 0; z < n; ++z) {
                for (int y = 1; y < n; ++y) {
                    for (int x = 0; x < n; ++x) {
                        vy[idY(x, y, z)] -= dt * (pressure[idC(x, y, z)] - pressure[idC(x, y - 1, z)]) * hInv;
                    }
                }
            }

#pragma omp parallel for
            for (int z = 1; z < n; ++z) {
                for (int y = 0; y < n; ++y) {
                    for (int x = 0; x < n; ++x) {
                        vz[idZ(x, y, z)] -= dt * (pressure[idC(x, y, z)] - pressure[idC(x, y, z - 1)]) * hInv;
                    }
                }
            }
        }

        enforceVelocityBoundaries();
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
                    if (qx * qx + qy * qy + qz * qz > 1.0f) continue;
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
#pragma omp parallel for
        for (int z = 0; z < n; ++z) {
            for (int y = 0; y < n; ++y) {
                for (int x = 0; x < n; ++x) {
                    const float f = buoyancyVal * smoke[idC(x, y, z)] * dtStep;
                    // Skip bottom (z==0) and top (z==n) solid wall faces
                    if (z > 0)     vz[idZ(x, y, z)] += buoyancySplit * f;
                    if (z < n - 1) vz[idZ(x, y, z + 1)] += buoyancySplit * f;
                }
            }
        }
    }

    // Pure smoke density dissipation (velocity remains divergence-free and inviscid)
    void applyDissipation(float dt) {
        const float sFactor = std::clamp(1.0f - smokeDecay * dt, 0.0f, 1.0f);

#pragma omp parallel for
        for (int i = 0; i < (int)smoke.size(); ++i) {
            smoke[i] *= sFactor;
        }
    }

    void step(float dt, float buoyancy, float sourceStrength, int projectIterations) {
        const float halfDt = 0.5f * dt;
        const int halfIters = std::max(1, projectIterations / 2);

        emit(sourceStrength, dt);
        applyBuoyancy(dt, buoyancy);
        enforceVelocityBoundaries();

        vx0 = vx; vy0 = vy; vz0 = vz;

        if (gpuSolver.advectionEnabled) {
            gpuSolver.advectVxGPU(halfDt, vx0.data(), gpuSolver.vxSz, vx0.data(), vy0.data(), vz0.data(), vx.data(), boxSize, cflMc, macCormackVel);
            vxTilde = vx;
            gpuSolver.advectVyGPU(halfDt, vy0.data(), gpuSolver.vySz, vx0.data(), vy0.data(), vz0.data(), vy.data(), boxSize, cflMc, macCormackVel);
            vyTilde = vy;
            gpuSolver.advectVzGPU(halfDt, vz0.data(), gpuSolver.vzSz, vx0.data(), vy0.data(), vz0.data(), vz.data(), boxSize, cflMc, macCormackVel);
            vzTilde = vz;
        } else {
            advectVx(halfDt, vx0, vx0, vy0, vz0, vx);
            vxTilde = vx;
            advectVy(halfDt, vy0, vx0, vy0, vz0, vy);
            vyTilde = vy;
            advectVz(halfDt, vz0, vx0, vy0, vz0, vz);
            vzTilde = vz;
        }

        project(halfDt, halfIters);

        vxDiv0 = vx; vyDiv0 = vy; vzDiv0 = vz;

#pragma omp parallel for
        for (int i = 0; i < (int)vx.size(); ++i) vxHat[i] = 2.f * vx[i] - vxTilde[i];
#pragma omp parallel for
        for (int i = 0; i < (int)vy.size(); ++i) vyHat[i] = 2.f * vy[i] - vyTilde[i];
#pragma omp parallel for
        for (int i = 0; i < (int)vz.size(); ++i) vzHat[i] = 2.f * vz[i] - vzTilde[i];

        if (gpuSolver.advectionEnabled) {
            gpuSolver.advectVxGPU(halfDt, vxHat.data(), gpuSolver.vxSz, vxDiv0.data(), vyDiv0.data(), vzDiv0.data(), vx.data(), boxSize, cflMc, macCormackVel);
            gpuSolver.advectVyGPU(halfDt, vyHat.data(), gpuSolver.vySz, vxDiv0.data(), vyDiv0.data(), vzDiv0.data(), vy.data(), boxSize, cflMc, macCormackVel);
            gpuSolver.advectVzGPU(halfDt, vzHat.data(), gpuSolver.vzSz, vxDiv0.data(), vyDiv0.data(), vzDiv0.data(), vz.data(), boxSize, cflMc, macCormackVel);
        } else {
            advectVx(halfDt, vxHat, vxDiv0, vyDiv0, vzDiv0, vx, macCormackVel);
            advectVy(halfDt, vyHat, vxDiv0, vyDiv0, vzDiv0, vy, macCormackVel);
            advectVz(halfDt, vzHat, vxDiv0, vyDiv0, vzDiv0, vz, macCormackVel);
        }

        project(halfDt, halfIters);

        if (gpuSolver.advectionEnabled) {
            const int N = n * n * n;
            gpuSolver.advectSmokeGPU(dt, smoke.data(), vx.data(), vy.data(), vz.data(), boxSize, cflMc, macCormackSmoke);
        } else {
            advectScalar(smoke, dt);
        }
        applyDissipation(dt);
    }

    void stepNoReflect(float dt, float buoyancy, float sourceStrength, int projectIterations) {
        emit(sourceStrength, dt);
        applyBuoyancy(dt, buoyancy);
        enforceVelocityBoundaries();

        vx0 = vx; vy0 = vy; vz0 = vz;

        if (gpuSolver.advectionEnabled) {
            gpuSolver.advectVxGPU(dt, vx0.data(), gpuSolver.vxSz, vx0.data(), vy0.data(), vz0.data(), vx.data(), boxSize, cflMc, macCormackVel);
            gpuSolver.advectVyGPU(dt, vy0.data(), gpuSolver.vySz, vx0.data(), vy0.data(), vz0.data(), vy.data(), boxSize, cflMc, macCormackVel);
            gpuSolver.advectVzGPU(dt, vz0.data(), gpuSolver.vzSz, vx0.data(), vy0.data(), vz0.data(), vz.data(), boxSize, cflMc, macCormackVel);
        } else {
            advectVx(dt, vx0, vx0, vy0, vz0, vx);
            advectVy(dt, vy0, vx0, vy0, vz0, vy);
            advectVz(dt, vz0, vx0, vy0, vz0, vz);
        }

        project(dt, projectIterations);

        if (gpuSolver.advectionEnabled) {
            const int N = n * n * n;
            gpuSolver.advectSmokeGPU(dt, smoke.data(), vx.data(), vy.data(), vz.data(), boxSize, cflMc, macCormackSmoke);
        } else {
            advectScalar(smoke, dt);
        }
        applyDissipation(dt);
    }

    float kineticEnergy() const {
        const float dx = h();
        double ke = 0.0;
        for (int z = 0; z < n; ++z) for (int y = 0; y < n; ++y) for (int x = 0; x <= n; ++x)
            ke += (double)vx[idX(x, y, z)] * vx[idX(x, y, z)];
        for (int z = 0; z < n; ++z) for (int y = 0; y <= n; ++y) for (int x = 0; x < n; ++x)
            ke += (double)vy[idY(x, y, z)] * vy[idY(x, y, z)];
        for (int z = 0; z <= n; ++z) for (int y = 0; y < n; ++y) for (int x = 0; x < n; ++x)
            ke += (double)vz[idZ(x, y, z)] * vz[idZ(x, y, z)];
        return float(0.5 * ke * dx * dx * dx);
    }

    float maxVelocity() const {
        float m = 0.f;
        for (float f : vx) m = std::max(m, std::abs(f));
        for (float f : vy) m = std::max(m, std::abs(f));
        for (float f : vz) m = std::max(m, std::abs(f));
        return m;
    }

    float computeDivergenceNorm() {
        const float hInv = 1.0f / h();
        float maxDiv = 0.f;
        for (int z = 0; z < n; ++z) {
            for (int y = 0; y < n; ++y) {
                for (int x = 0; x < n; ++x) {
                    float div = std::abs((vx[idX(x + 1, y, z)] - vx[idX(x, y, z)] +
                        vy[idY(x, y + 1, z)] - vy[idY(x, y, z)] +
                        vz[idZ(x, y, z + 1)] - vz[idZ(x, y, z)]) * hInv);
                    maxDiv = std::max(maxDiv, div);
                }
            }
        }
        return maxDiv;
    }

    void rebox(float oldSz, float newSz) {
        if (newSz == oldSz) return;
        boxSize = newSz;
        resetState();
    }

    void resetState() {
        std::fill(smoke.begin(), smoke.end(), 0.f);
        std::fill(vx.begin(), vx.end(), 0.f);
        std::fill(vy.begin(), vy.end(), 0.f);
        std::fill(vz.begin(), vz.end(), 0.f);
        std::fill(pressure.begin(), pressure.end(), 0.f);
        std::fill(divergence.begin(), divergence.end(), 0.f);
        std::fill(fwdScratch.begin(), fwdScratch.end(), 0.f);
        std::fill(backScratch.begin(), backScratch.end(), 0.f);
        std::fill(vx0.begin(), vx0.end(), 0.f);
        std::fill(vy0.begin(), vy0.end(), 0.f);
        std::fill(vz0.begin(), vz0.end(), 0.f);
        std::fill(vxTilde.begin(), vxTilde.end(), 0.f);
        std::fill(vyTilde.begin(), vyTilde.end(), 0.f);
        std::fill(vzTilde.begin(), vzTilde.end(), 0.f);
        std::fill(vxHat.begin(), vxHat.end(), 0.f);
        std::fill(vyHat.begin(), vyHat.end(), 0.f);
        std::fill(vzHat.begin(), vzHat.end(), 0.f);
        std::fill(vxDiv0.begin(), vxDiv0.end(), 0.f);
        std::fill(vyDiv0.begin(), vyDiv0.end(), 0.f);
        std::fill(vzDiv0.begin(), vzDiv0.end(), 0.f);
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

        glTexImage3D(GL_TEXTURE_3D, 0, GL_R16F, n, n, n, 0, GL_RED, GL_FLOAT, nullptr);

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
            uniform vec3 u_LightDir;
            uniform float u_ShadowStr;
            uniform float u_ShadowStep;
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
                    if (d > 0.001) {
                        hitSmoke = true;
                        float shadow = 1.0;
                        if (u_ShadowStr > 0.0) {
                            vec3 spos = ro + rd * t;
                            float sa = 0.0;
                            for (float st = u_ShadowStep; st < u_BoxSize * 2.0; st += u_ShadowStep) {
                                vec3 sp = spos + u_LightDir * st;
                                if (sp.x < 0.0 || sp.x > u_BoxSize || sp.y < 0.0 || sp.y > u_BoxSize || sp.z < 0.0 || sp.z > u_BoxSize) break;
                                sa += texture(u_Volume, sp / u_BoxSize).r * u_ShadowStep * u_AlphaMul * u_ShadowStr;
                                if (sa > 3.0) break;
                            }
                            shadow = exp(-sa);
                        }
                        float alpha = clamp(d * step * u_AlphaMul, 0.0, 1.0);
                        col.rgb += (1.0 - col.a) * smokeCol * alpha * shadow;
                        col.a += (1.0 - col.a) * alpha;
                    }
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

    void upload(const std::vector<float>& density) {
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_3D, volumeTex);
        glTexSubImage3D(GL_TEXTURE_3D, 0, 0, 0, 0, n, n, n, GL_RED, GL_FLOAT, density.data());
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

    void render(int w, int h, const float* mvp, float cx, float cy, float cz, float stepScale, float boxSize, float alphaMul, float lightX, float lightY, float lightZ, float shadowStr, float shadowStep, FILE* logFile) {
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
        glUniform3f(glGetUniformLocation(prog, "u_LightDir"), lightX, lightY, lightZ);
        glUniform1f(glGetUniformLocation(prog, "u_ShadowStr"), shadowStr);
        glUniform1f(glGetUniformLocation(prog, "u_ShadowStep"), shadowStep);
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
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
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

    bool paused = false, useReflection = true, debugPrintOn = true;
    int maxFrames = 600;
    float buoyancy = 2.5f, sourceStrength = 1.0f;
    int projectIterations = sim.n * 2;
    float alphaMul = 30.0f;
    int debugFrame = 0;
    const float maxDt = 0.033f;
    float azimuth = -1.2f, elevation = 0.3f, dist = sim.boxSize * 2.5f, stepScale = 0.5f;
    float prevBoxSize = sim.boxSize;
    float lightX = 0.3f, lightY = 0.4f, lightZ = 1.0f;
    float shadowStr = 0.3f, shadowStep = 0.1f;
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
                float ke = sim.kineticEnergy();
                float maxDiv = sim.computeDivergenceNorm();
                float maxV = sim.maxVelocity();
                logPrint("f %4d | ke=%.6f maxDiv=%.6f maxV=%.3f smoke=%.3f | %s\n",
                    debugFrame, ke, maxDiv, maxV, dmax,
                    useReflection ? "REFLECT" : "NOREFLECT");
            }
            if (maxFrames > 0 && debugFrame >= maxFrames) {
                logPrint("--- AUTO-QUIT after %d frames ---\n", maxFrames);
                glfwSetWindowShouldClose(window, GLFW_TRUE);
            }
        }

        int w, h; glfwGetFramebufferSize(window, &w, &h);
        glBindFramebuffer(GL_FRAMEBUFFER, 0);
        glClearColor(0.03f, 0.04f, 0.07f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        renderer.upload(sim.smoke);
        float mvp[16], cx, cy, cz;
        buildMVP(mvp, cx, cy, cz, azimuth, elevation, dist, float(w) / float(h), sim.boxSize);
        float ld = 1.0f / std::sqrt(lightX * lightX + lightY * lightY + lightZ * lightZ);
        renderer.render(w, h, mvp, cx, cy, cz, stepScale, sim.boxSize, alphaMul,
            lightX * ld, lightY * ld, lightZ * ld, shadowStr, shadowStep, logFile);

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
            sim.resetState();
        }
        ImGui::Checkbox("Reflection", &useReflection);
        ImGui::Checkbox("MC smoke", &sim.macCormackSmoke);
        ImGui::Checkbox("MC vel", &sim.macCormackVel);
        ImGui::Checkbox("Debug print", &debugPrintOn);
        ImGui::Checkbox("GPU SOR", &sim.gpuSolver.enabled);
        ImGui::Checkbox("GPU Advection", &sim.gpuSolver.advectionEnabled);
        ImGui::SliderFloat("Buoyancy", &buoyancy, 0.f, 10.f);
        ImGui::SliderFloat("Source", &sourceStrength, 0.f, 3.f);
        ImGui::SliderFloat("Smoke decay", &sim.smokeDecay, 0.f, 2.f);
        ImGui::SliderFloat("Box size", &sim.boxSize, 0.5f, 8.f);
        ImGui::SliderFloat("Emitt radius", &sim.emitterRadius, 0.01f, 0.2f);
        ImGui::SliderFloat("Stir", &sim.stirStrength, 0.f, 2.f);
        ImGui::SliderFloat("Alpha mul", &alphaMul, 5.f, 100.f);
        ImGui::SliderInt("SOR its", &projectIterations, 10, 500);
        ImGui::SliderFloat("MC CFL", &sim.cflMc, 0.5f, 10.f);
        ImGui::SliderFloat("Step", &stepScale, 0.1f, 3.f);
        ImGui::SliderFloat("Shadow str", &shadowStr, 0.f, 5.f);
        ImGui::SliderFloat("Shadow step", &shadowStep, 0.05f, 1.f);
        ImGui::SliderFloat3("Light dir", &lightX, -1.f, 1.f);
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
    sim.gpuSolver.shutdown();
    renderer.shutdown();
    glfwDestroyWindow(window);
    glfwTerminate();
    if (logFile) std::fclose(logFile);
}