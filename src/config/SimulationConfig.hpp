#pragma once

#include "domain/IDomain.hpp"

#include <filesystem>
#include <iosfwd>

enum class DomainType { Mac, Collocated };
enum class AdvectionType { SemiLagrangian, MacCormack };
enum class ComputeBackend { Cpu, Gpu };
enum class FluidSolverType { Simple, Reflection };
enum class ProjectionSolverType { CpuSor, GpuSor };

struct SimulationConfig {
    DomainType domain = DomainType::Mac;
    AdvectionType advection = AdvectionType::SemiLagrangian;
    Dimension dimension = Dimension::D2;
    FluidSolverType fluidSolver = FluidSolverType::Simple;
    ProjectionSolverType projectionSolver = ProjectionSolverType::CpuSor;
    ComputeBackend advectionBackend = ComputeBackend::Cpu;

    int resolution = 64;
    float boxSize = 4.f;
    float emitterRadius = .08f;
    float sourceStrength = 4.f;
    float buoyancy = 3.f;
    float smokeDecay = .08f;
    int sorIterations = 600;
    float sorOmega = 1.9f;
    int maxFrames = 0; // Zero means run until the window closes.
};

// Every supported key must occur exactly once. Unknown, duplicate, malformed,
// and out-of-range values throw. Module availability is checked by the runtime registry.
[[nodiscard]] SimulationConfig parseSimulationConfig(std::istream& input);
[[nodiscard]] SimulationConfig loadSimulationConfig(const std::filesystem::path& path);
void validateSimulationConfig(const SimulationConfig& config);
