#include "config/SimulationConfig.hpp"

#include <cmath>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace {
void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}
template<class F> void rejects(F&& function) {
    try { function(); }
    catch (const std::invalid_argument&) { return; }
    throw std::runtime_error("Invalid configuration was accepted");
}
std::string valid2D() {
    return R"cfg(
# Complete 2D composition, with whitespace and inline comments.
domain = mac
advection=semi_lagrangian
dimension = 2
fluid_solver = simple
projection_solver = cpu_sor
advection_backend = cpu
resolution = 48
box_size = 3.5
emitter_radius = 0.1
source_strength = 2.25
buoyancy = -1.5 # downward is a valid physical choice
smoke_decay = 0.125
sor_iterations = 700
sor_omega = 1.85
max_frames = 17
)cfg";
}
std::string replace(std::string text, const std::string& from, const std::string& to) {
    const auto position = text.find(from);
    if (position == std::string::npos) throw std::runtime_error("Broken test replacement");
    text.replace(position, from.size(), to);
    return text;
}
SimulationConfig parse(const std::string& text) {
    std::istringstream input(text);
    return parseSimulationConfig(input);
}
}

int main() {
    try {
        const auto two = parse(valid2D());
        require(two.domain == DomainType::Mac && two.dimension == Dimension::D2, "2D domain parse");
        require(two.advection == AdvectionType::SemiLagrangian && two.advectionBackend == ComputeBackend::Cpu,
                "2D advection parse");
        require(two.fluidSolver == FluidSolverType::Simple && two.projectionSolver == ProjectionSolverType::CpuSor,
                "2D solver parse");
        require(two.resolution == 48 && std::abs(two.boxSize - 3.5f) < 1e-6f &&
                std::abs(two.emitterRadius - .1f) < 1e-6f && two.sorIterations == 700 && two.maxFrames == 17,
                "Numerical value parse");

        auto threeText = replace(valid2D(), "dimension = 2", "dimension = 3");
        threeText = replace(threeText, "advection=semi_lagrangian", "advection=maccormack");
        threeText = replace(threeText, "fluid_solver = simple", "fluid_solver = reflection");
        threeText = replace(threeText, "projection_solver = cpu_sor", "projection_solver = gpu_sor");
        threeText = replace(threeText, "advection_backend = cpu", "advection_backend = gpu");
        const auto three = parse(threeText);
        require(three.dimension == Dimension::D3 && three.advection == AdvectionType::MacCormack &&
                three.fluidSolver == FluidSolverType::Reflection &&
                three.projectionSolver == ProjectionSolverType::GpuSor &&
                three.advectionBackend == ComputeBackend::Gpu, "3D composition parse");

        rejects([&] { parse(valid2D() + "domain=mac\n"); });
        rejects([&] { parse(replace(valid2D(), "domain = mac", "domain = staggered")); });
        rejects([&] { parse(replace(valid2D(), "dimension = 2", "dimension = 4")); });
        rejects([&] { parse(replace(valid2D(), "resolution = 48\n", "")); });
        rejects([&] { parse(replace(valid2D(), "resolution = 48", "resolution = 0")); });
        rejects([&] { parse(replace(valid2D(), "sor_omega = 1.85", "sor_omega = 2")); });
        rejects([&] { parse(replace(valid2D(), "max_frames = 17", "surprise = 17")); });
        require(parse(replace(valid2D(), "domain = mac", "domain = collocated")).domain ==
                    DomainType::Collocated, "Recognized domain belongs to the parser");
        require(parse(replace(valid2D(), "advection=semi_lagrangian", "advection=maccormack")).advection ==
                    AdvectionType::MacCormack, "Recognized advection belongs to the parser");
        require(parse(replace(valid2D(), "fluid_solver = simple", "fluid_solver = reflection")).fluidSolver ==
                    FluidSolverType::Reflection, "Recognized fluid algorithm belongs to the parser");
        require(parse(replace(valid2D(), "projection_solver = cpu_sor", "projection_solver = gpu_sor")).projectionSolver ==
                    ProjectionSolverType::GpuSor, "Recognized projection belongs to the parser");
        require(parse(replace(valid2D(), "advection_backend = cpu", "advection_backend = gpu")).advectionBackend ==
                    ComputeBackend::Gpu, "Recognized backend belongs to the parser");
        std::cout << "PASS configuration syntax, enums, and numeric validation\n";
    } catch (const std::exception& error) {
        std::cerr << "FAIL configuration: " << error.what() << '\n';
        return 1;
    }
}
