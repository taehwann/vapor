#include "gpu/GpuMacBackend.hpp"
#include "fluid/MacGridState.hpp"
#include "graphics/GlLoader.hpp"
#include "ShaderSources.hpp"

#include <algorithm>
#include <cstdio>
#include <vector>
#include <stdexcept>
#include <string>

static GLuint compileComputeShader(const char* src) {
    GLuint sh = glCreateShader(GL_COMPUTE_SHADER);
    glShaderSource(sh, 1, &src, nullptr);
    glCompileShader(sh);
    int ok = 0; glGetShaderiv(sh, GL_COMPILE_STATUS, &ok);
    if (!ok) { char log[1024]; glGetShaderInfoLog(sh, 1024, nullptr, log); glDeleteShader(sh); throw std::runtime_error(std::string("Compute compile: ") + log); }
    GLuint prog = glCreateProgram();
    glAttachShader(prog, sh);
    glLinkProgram(prog);
    int linkOk = 0; glGetProgramiv(prog, GL_LINK_STATUS, &linkOk);
    if (!linkOk) { char log[1024]; glGetProgramInfoLog(prog, 1024, nullptr, log); glDeleteShader(sh); glDeleteProgram(prog); throw std::runtime_error(std::string("Compute link: ") + log); }
    glDeleteShader(sh);
    return prog;
}

GpuMacBackend::~GpuMacBackend() { shutdown(); }

void GpuMacBackend::validateGrid(const MacGridState& state) const {
    if (!initialized_ || state.resolution() != n ||
        state.velocityX().size() != static_cast<std::size_t>(vxSz) ||
        state.velocityY().size() != static_cast<std::size_t>(vySz) ||
        state.velocityZ().size() != static_cast<std::size_t>(vzSz)) {
        throw std::invalid_argument("GPU backend is uninitialized or has a different grid layout");
    }
}

void GpuMacBackend::init(const MacGridState& state) {
    if (!glfwGetCurrentContext() || !glDispatchCompute)
        throw std::runtime_error("GPU initialization requires a current OpenGL 4.3 context and loaded functions");
    shutdown();
    const int res = state.resolution();
    const int vxSize = static_cast<int>(state.velocityX().size());
    const int vySize = static_cast<int>(state.velocityY().size());
    const int vzSize = static_cast<int>(state.velocityZ().size());
    try {
    n = res; vxSz = vxSize; vySz = vySize; vzSz = vzSize;

    const char* divSrc = ShaderSources::divSrc;

    const char* rbgsSrc = ShaderSources::rbgsSrc;

    const char* correctSrc = ShaderSources::correctSrc;

    progDiv = compileComputeShader(divSrc);
    progRBGS = compileComputeShader(rbgsSrc);
    progCorrect = compileComputeShader(correctSrc);

    const char* advectSlSrc = ShaderSources::advectSlSrc;

    const char* advectMcSrc = ShaderSources::advectMcSrc;

    const char* copySrc = ShaderSources::copySrc;

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

    initialized_ = true;
    } catch (...) { shutdown(); throw; }
}

void GpuMacBackend::project(float dt, int iterations, float omega, float hInv, float h,
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

void GpuMacBackend::uploadVelocity(const float* vxData, const float* vyData, const float* vzData) {
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, ssboVx);
    glBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, (GLsizeiptr)vxSz * sizeof(float), vxData);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, ssboVy);
    glBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, (GLsizeiptr)vySz * sizeof(float), vyData);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, ssboVz);
    glBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, (GLsizeiptr)vzSz * sizeof(float), vzData);
}

void GpuMacBackend::dispatchSL(int mode, int count, float dx, float dt) {
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

void GpuMacBackend::dispatchMC(int mode, float dx, float dt, float boxSize, float cflMc, GLuint srcBo, GLuint outBo) {
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

void GpuMacBackend::copySSBO(GLuint dstBo, GLuint srcBo, int count) {
    glUseProgram(progCopy);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, srcBo);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, dstBo);
    glUniform1i(glGetUniformLocation(progCopy, "u_count"), count);
    glDispatchCompute((count + 63) / 64, 1, 1);
    glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
}

void GpuMacBackend::advectSmokeGPU(float dt, float* smokeData, const float* vxVel, const float* vyVel, const float* vzVel, float boxSize, float cflMc, bool useMC) {
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

void GpuMacBackend::advectVxGPU(float dt, const float* vxSrc, int srcSz, const float* vxVel, const float* vyVel, const float* vzVel, float* vxOut, float boxSize, float cflMc, bool useMC) {
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

void GpuMacBackend::advectVyGPU(float dt, const float* vySrc, int srcSz, const float* vxVel, const float* vyVel, const float* vzVel, float* vyOut, float boxSize, float cflMc, bool useMC) {
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

void GpuMacBackend::advectVzGPU(float dt, const float* vzSrc, int srcSz, const float* vxVel, const float* vyVel, const float* vzVel, float* vzOut, float boxSize, float cflMc, bool useMC) {
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

void GpuMacBackend::shutdown() {
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
    progDiv = progRBGS = progCorrect = progAdvectSL = progAdvectMC = progCopy = 0;
    ssboP = ssboDiv = ssboVx = ssboVy = ssboVz = 0;
    ssboSmoke = ssboFwd = ssboAdvectSrc = ssboAdvectOut = 0;
    initialized_ = false;
}
