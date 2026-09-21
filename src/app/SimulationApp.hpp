#pragma once
#include "common/MacAlgorithm.hpp"
#include <filesystem>

enum class SolverChoice {
    Mac2DSimple,
    Mac2DReflectionSL,
    Mac2DReflectionMC,
    Mac3DSimpleSLCPU,
    Mac3DSimpleSLGPU,
    Mac3DSimpleMCCPU,
    Mac3DSimpleMCGPU,
    Mac3DReflectionSLCPU,
    Mac3DReflectionSLGPU,
    Mac3DReflectionMCCPU,
    Mac3DReflectionMCGPU,
    Simplicial2DTeapot,
    BunnyMesh
};

// Application timing/assets only. Physics defaults belong to concrete solvers.
struct RunOptions {
    std::filesystem::path assetDirectory = ".";
    float fixedDt = 1.f / 60.f;
    int maxFrames = 0;
    bool hidden = false; // Integration tests only.
};

class SimulationApp {
  public:
    int run(const std::filesystem::path& assetDirectory);
    int runSolver(SolverChoice choice, const RunOptions& options = {});

  private:
    int run2D(SolverChoice choice, const RunOptions& options);
    int run3D(MacAlgorithm algorithm, bool useGpu, const RunOptions& options);
    int runSimplicial3D(const RunOptions& options);
    int runMesh3D(const RunOptions& options);
};
