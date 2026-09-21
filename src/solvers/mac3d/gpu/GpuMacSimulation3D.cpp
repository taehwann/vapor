#include "solvers/mac3d/gpu/GpuMacSimulation3D.hpp"
#include "ShaderSources.hpp"
#include "graphics/GlLoader.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

namespace {
unsigned compile(const char *source) {
    GLuint shader = glCreateShader(GL_COMPUTE_SHADER), program = 0;
    try {
        glShaderSource(shader, 1, &source, nullptr);
        glCompileShader(shader);
        GLint ok;
        glGetShaderiv(shader, GL_COMPILE_STATUS, &ok);
        if (!ok) {
            char log[2048];
            glGetShaderInfoLog(shader, sizeof(log), nullptr, log);
            throw std::runtime_error(log);
        }
        program = glCreateProgram();
        glAttachShader(program, shader);
        glLinkProgram(program);
        glGetProgramiv(program, GL_LINK_STATUS, &ok);
        if (!ok) {
            char log[2048];
            glGetProgramInfoLog(program, sizeof(log), nullptr, log);
            throw std::runtime_error(log);
        }
    } catch (...) {
        if (program)
            glDeleteProgram(program);
        glDeleteShader(shader);
        throw;
    }
    glDeleteShader(shader);
    return program;
}
void allocate(unsigned &buffer, int count) {
    glGenBuffers(1, &buffer);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, buffer);
    glBufferData(GL_SHADER_STORAGE_BUFFER, GLsizeiptr(count) * sizeof(float), nullptr, GL_DYNAMIC_DRAW);
}
} // namespace
GpuMacSimulation3D::GpuMacSimulation3D(MacAlgorithm algorithm)
    : GpuMacSimulation3D(32, 4.5f, MacGridParameters3D(algorithm)) {}

GpuMacSimulation3D::GpuMacSimulation3D(int n, float box, MacGridParameters3D p)
    : parameters_(p), boxSize_(box), randomState_(p.randomSeed) {
    try {
        gpu_.init(MacGridState3D(n, box));
        operationsProgram_ = compile(ShaderSources::residentOps);
        for (auto &b : next_)
            allocate(b, gpu_.maxSz);
        allocate(noiseBuffer_, 1);
        if (glGetError() != GL_NO_ERROR)
            throw std::runtime_error("GPU resident allocation failed");
        resetState();
    } catch (...) {
        release();
        throw;
    }
}
GpuMacSimulation3D::~GpuMacSimulation3D() { release(); }
void GpuMacSimulation3D::release() {
    for (auto &b : next_)
        if (b) {
            glDeleteBuffers(1, &b);
            b = 0;
        }
    for (auto &b : reflected_)
        if (b) {
            glDeleteBuffers(1, &b);
            b = 0;
        }
    if (noiseBuffer_) {
        glDeleteBuffers(1, &noiseBuffer_);
        noiseBuffer_ = 0;
    }
    if (operationsProgram_) {
        glDeleteProgram(operationsProgram_);
        operationsProgram_ = 0;
    }
}
std::array<unsigned, 3> GpuMacSimulation3D::velocityBuffers() const {
    return {gpu_.ssboVx, gpu_.ssboVy, gpu_.ssboVz};
}
void GpuMacSimulation3D::resetState() {
    const float zero = 0;
    glMemoryBarrier(GL_BUFFER_UPDATE_BARRIER_BIT);
    for (auto b : {gpu_.ssboVx, gpu_.ssboVy, gpu_.ssboVz, gpu_.ssboSmoke, gpu_.ssboP, gpu_.ssboDiv}) {
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, b);
        glClearBufferData(GL_SHADER_STORAGE_BUFFER, GL_R32F, GL_RED, GL_FLOAT, &zero);
    }
    randomState_ = parameters_.randomSeed;
    lastProjection_ = {};
}
void GpuMacSimulation3D::setBoxSize(float size) {
    if (!std::isfinite(size) || size <= 0)
        throw std::invalid_argument("Invalid GPU box size");
    if (size != boxSize_) {
        boxSize_ = size;
        resetState();
    }
}
RenderData GpuMacSimulation3D::renderData() const {
    return {{}, gpu_.n, gpu_.n, gpu_.n, boxSize_, gpu_.ssboSmoke};
}
void GpuMacSimulation3D::upload(const MacGridState3D &s) {
    gpu_.validateGrid(s);
    if (s.boxSize() != boxSize_)
        throw std::invalid_argument("GPU snapshot box mismatch");
    glMemoryBarrier(GL_BUFFER_UPDATE_BARRIER_BIT);
    const auto buffers = velocityBuffers();
    const std::array<std::span<const float>, 3> fields = {s.velocityX(), s.velocityY(), s.velocityZ()};
    for (int i = 0; i < 3; ++i) {
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, buffers[i]);
        glBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, fields[i].size_bytes(), fields[i].data());
        uploadedBytes_ += fields[i].size_bytes();
    }
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, gpu_.ssboSmoke);
    glBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, s.density().size_bytes(), s.density().data());
    uploadedBytes_ += s.density().size_bytes();
}
void GpuMacSimulation3D::download(MacGridState3D &s) const {
    gpu_.validateGrid(s);
    if (s.boxSize() != boxSize_)
        throw std::invalid_argument("GPU snapshot box mismatch");
    glMemoryBarrier(GL_BUFFER_UPDATE_BARRIER_BIT);
    const auto buffers = velocityBuffers();
    const std::array<std::span<float>, 3> fields = {s.velocityX(), s.velocityY(), s.velocityZ()};
    for (int i = 0; i < 3; ++i) {
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, buffers[i]);
        glGetBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, fields[i].size_bytes(), fields[i].data());
        downloadedBytes_ += fields[i].size_bytes();
    }
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, gpu_.ssboSmoke);
    glGetBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, s.density().size_bytes(), s.density().data());
    downloadedBytes_ += s.density().size_bytes();
}
void GpuMacSimulation3D::operations(int mode) {
    glUseProgram(operationsProgram_);
    glUniform1i(glGetUniformLocation(operationsProgram_, "u_n"), gpu_.n);
    glUniform1i(glGetUniformLocation(operationsProgram_, "u_mode"), mode);
    auto buffers = velocityBuffers();
    for (int i = 0; i < 3; ++i)
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, i, buffers[i]);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 3, gpu_.ssboSmoke);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 4, noiseBuffer_);
    if (mode == 3)
        for (int i = 0; i < 3; ++i)
            glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 5 + i, reflected_[i]);
    glDispatchCompute((gpu_.n + 8) / 8, (gpu_.n + 8) / 8, (gpu_.n + 4) / 4);
    glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
}
void GpuMacSimulation3D::emit(float dt) {
    const auto &p = parameters_;
    if (p.sourceStrength <= 0)
        return;
    const int n = gpu_.n;
    const float h = boxSize_ / n, r = p.emitterRadius * boxSize_;
    const std::array<float, 3> center = {p.emitterCenterX * boxSize_, p.emitterCenterY * boxSize_,
                                         p.emitterCenterZ * boxSize_};
    std::array<int, 3> lo, extent;
    const int radius = int(r / h + 1);
    for (int i = 0; i < 3; ++i) {
        int c = int(center[i] / h);
        lo[i] = std::max(0, c - radius);
        extent[i] = std::max(0, std::min(n, c + radius) - lo[i]);
    }
    const auto count = std::size_t(extent[0]) * extent[1] * extent[2];
    if (!count)
        return;
    noise_.assign(count, 1e20f);
    // Preserve the existing seeded RNG sequence and cell traversal. Only this
    // compact inlet patch crosses to the GPU; physical fields remain resident.
    for (int z = 0; z < extent[2]; ++z)
        for (int y = 0; y < extent[1]; ++y)
            for (int x = 0; x < extent[0]; ++x) {
                float qx = ((x + lo[0] + .5f) * h - center[0]) / r,
                      qy = ((y + lo[1] + .5f) * h - center[1]) / r,
                      qz = ((z + lo[2] + .5f) * h - center[2]) / r;
                if (qx * qx + qy * qy + qz * qz > 1)
                    continue;
                randomState_ = randomState_ * 1664525u + 1013904223u;
                float random = (randomState_ >> 8) * (1.f / 16777216.f);
                noise_[x + extent[0] * (y + extent[1] * z)] = (random - .5f) * p.stirStrength * (dt * 60.f);
            }
    glMemoryBarrier(GL_BUFFER_UPDATE_BARRIER_BIT);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, noiseBuffer_);
    if (count > noiseCapacity_) {
        glBufferData(GL_SHADER_STORAGE_BUFFER, count * sizeof(float), noise_.data(), GL_DYNAMIC_DRAW);
        noiseCapacity_ = count;
    } else
        glBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, count * sizeof(float), noise_.data());
    uploadedBytes_ += count * sizeof(float);
    glUseProgram(operationsProgram_);
    glUniform3i(glGetUniformLocation(operationsProgram_, "u_emitLo"), lo[0], lo[1], lo[2]);
    glUniform3i(glGetUniformLocation(operationsProgram_, "u_emitSize"), extent[0], extent[1], extent[2]);
    glUniform1f(glGetUniformLocation(operationsProgram_, "u_kick"),
                p.emitterKickMultiplier * (p.emitterKickBase + p.emitterKickScale * p.sourceStrength));
    operations(0);
}
void GpuMacSimulation3D::advectVelocity(float dt, const std::array<unsigned, 3> &source) {
    const std::array<int, 3> counts = {gpu_.vxSz, gpu_.vySz, gpu_.vzSz};
    for (int i = 0; i < 3; ++i) {
        gpu_.dispatchSL(i + 1, counts[i], boxSize_ / gpu_.n, dt, source[i],
                        parameters_.macCormackVel ? gpu_.ssboFwd : next_[i]);
        if (parameters_.macCormackVel)
            gpu_.dispatchMC(i + 1, boxSize_ / gpu_.n, dt, boxSize_, parameters_.cflMc, source[i], next_[i]);
    }
    std::swap(gpu_.ssboVx, next_[0]);
    std::swap(gpu_.ssboVy, next_[1]);
    std::swap(gpu_.ssboVz, next_[2]);
}
void GpuMacSimulation3D::project(float dt, int iterations) {
    operations(2);
    gpu_.projectDevice(dt, iterations, parameters_.sorOmega, gpu_.n / boxSize_, boxSize_ / gpu_.n);
    operations(2);
    lastProjection_ = {};
    lastProjection_.linearSolve.iterations = iterations;
}
void GpuMacSimulation3D::advance(float dt) {
    const auto &p = parameters_;
    if (!std::isfinite(p.emitterKickMultiplier) || p.emitterKickMultiplier < 0.f)
        throw std::invalid_argument("Invalid emitter kick multiplier");
    if (!std::isfinite(dt) || dt <= 0 || !std::isfinite(p.sourceStrength) || p.sourceStrength < 0 ||
        !std::isfinite(p.emitterRadius) || p.emitterRadius <= 0 || p.emitterRadius > 1 ||
        !std::isfinite(p.buoyancy) || !std::isfinite(p.smokeDecay) || p.smokeDecay < 0 ||
        !std::isfinite(p.sorOmega) || p.sorOmega <= 0 || p.sorOmega >= 2 || p.projectIterations < 0 ||
        !std::isfinite(p.stirStrength) || !std::isfinite(p.cflMc) || p.cflMc <= 0)
        throw std::invalid_argument("Invalid resident GPU step parameters");
    emit(dt);
    glUseProgram(operationsProgram_);
    glUniform1f(glGetUniformLocation(operationsProgram_, "u_force"), p.buoyancy * dt);
    glUniform1f(glGetUniformLocation(operationsProgram_, "u_split"), p.buoyancySplit);
    operations(1);
    advectVelocity(p.reflection ? .5f * dt : dt, velocityBuffers());
    if (p.reflection) {
        const auto buffers = velocityBuffers();
        const std::array<int, 3> counts = {gpu_.vxSz, gpu_.vySz, gpu_.vzSz};
        for (int i = 0; i < 3; ++i) {
            if (!reflected_[i])
                allocate(reflected_[i], counts[i]);
            gpu_.copySSBO(reflected_[i], buffers[i], counts[i]);
        }
        const int iters = std::max(1, p.projectIterations / 2);
        project(.5f * dt, iters);
        operations(3);
        advectVelocity(.5f * dt, reflected_);
        project(.5f * dt, iters);
    } else
        project(dt, p.projectIterations);
    gpu_.dispatchSL(0, gpu_.n * gpu_.n * gpu_.n, boxSize_ / gpu_.n, dt, gpu_.ssboSmoke,
                    p.macCormackSmoke ? gpu_.ssboFwd : gpu_.ssboAdvectOut);
    if (p.macCormackSmoke)
        gpu_.dispatchMC(0, boxSize_ / gpu_.n, dt, boxSize_, p.cflMc, gpu_.ssboSmoke, gpu_.ssboAdvectOut);
    std::swap(gpu_.ssboSmoke, gpu_.ssboAdvectOut);
    glUseProgram(operationsProgram_);
    glUniform1f(glGetUniformLocation(operationsProgram_, "u_decay"),
                std::clamp(1.f - p.smokeDecay * dt, 0.f, 1.f));
    operations(4);
}
