#pragma once
#include "renderer/RenderData.hpp"
#include "solvers/mac3d/MacGridFluidSolver3D.hpp"
#include "solvers/mac3d/gpu/GpuMacBackend.hpp"
#include <array>

// Fully GPU-resident coordinator. CPU snapshots are explicit, never implicit
// in advance() or renderData(). Requires a current OpenGL 4.3 context.
class GpuMacSimulation3D final {
  public:
    explicit GpuMacSimulation3D(MacAlgorithm algorithm);
    GpuMacSimulation3D(int resolution, float boxSize, MacGridParameters3D parameters);
    ~GpuMacSimulation3D();
    GpuMacSimulation3D(const GpuMacSimulation3D &) = delete;
    GpuMacSimulation3D &operator=(const GpuMacSimulation3D &) = delete;
    void advance(float dt);
    void resetState();
    void setBoxSize(float size);
    float boxSize() const { return boxSize_; }
    MacGridParameters3D &parameters() { return parameters_; }
    const ProjectionResult &lastProjection() const { return lastProjection_; }
    RenderData renderData() const;
    void upload(const MacGridState3D &state);
    void download(MacGridState3D &state) const;
    std::uint64_t uploadedBytes() const { return uploadedBytes_; }
    std::uint64_t downloadedBytes() const { return downloadedBytes_; }

  private:
    void release();
    void operations(int mode);
    void emit(float dt);
    void advectVelocity(float dt, const std::array<unsigned int, 3> &source);
    void project(float dt, int iterations);
    std::array<unsigned int, 3> velocityBuffers() const;
    GpuMacBackend gpu_;
    MacGridParameters3D parameters_;
    float boxSize_;
    std::uint32_t randomState_;
    unsigned int operationsProgram_ = 0, noiseBuffer_ = 0;
    std::array<unsigned int, 3> next_{}, reflected_{};
    std::vector<float> noise_;
    std::size_t noiseCapacity_ = 0;
    ProjectionResult lastProjection_;
    std::uint64_t uploadedBytes_ = 0;
    mutable std::uint64_t downloadedBytes_ = 0;
};
