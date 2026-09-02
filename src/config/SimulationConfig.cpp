#include "config/SimulationConfig.hpp"

#include <array>
#include <charconv>
#include <cmath>
#include <fstream>
#include <set>
#include <stdexcept>
#include <string>
#include <string_view>

using namespace std::string_view_literals;

namespace {
std::string_view trim(std::string_view value) {
    const auto first = value.find_first_not_of(" \t\r\n");
    if (first == std::string_view::npos) return {};
    return value.substr(first, value.find_last_not_of(" \t\r\n") - first + 1);
}
template<class Number>
Number number(std::string_view text, std::string_view key, int line) {
    Number result{};
    const auto parsed = std::from_chars(text.data(), text.data() + text.size(), result);
    if (parsed.ec != std::errc{} || parsed.ptr != text.data() + text.size())
        throw std::invalid_argument("Invalid value for '" + std::string(key) + "' on line " + std::to_string(line));
    return result;
}
template<class Enum, std::size_t N>
Enum choice(std::string_view text, std::string_view key, int line,
            const std::array<std::pair<std::string_view, Enum>, N>& values) {
    for (const auto& [name, value] : values) if (text == name) return value;
    throw std::invalid_argument("Unknown value '" + std::string(text) + "' for '" +
                                std::string(key) + "' on line " + std::to_string(line));
}
}

SimulationConfig parseSimulationConfig(std::istream& input) {
    SimulationConfig config;
    std::set<std::string> seen;
    std::string storage;
    int line = 0;
    while (std::getline(input, storage)) {
        ++line;
        auto text = trim(std::string_view(storage).substr(0, storage.find('#')));
        if (text.empty()) continue;
        const auto equal = text.find('=');
        if (equal == std::string_view::npos || text.find('=', equal + 1) != std::string_view::npos)
            throw std::invalid_argument("Expected key=value on line " + std::to_string(line));
        const auto key = trim(text.substr(0, equal));
        const auto value = trim(text.substr(equal + 1));
        if (key.empty() || value.empty()) throw std::invalid_argument("Empty key or value on line " + std::to_string(line));
        if (!seen.insert(std::string(key)).second)
            throw std::invalid_argument("Duplicate key '" + std::string(key) + "' on line " + std::to_string(line));

        if (key == "domain") config.domain = choice(value, key, line,
            std::array{std::pair{"mac"sv, DomainType::Mac}, std::pair{"collocated"sv, DomainType::Collocated}});
        else if (key == "advection") config.advection = choice(value, key, line,
            std::array{std::pair{"semi_lagrangian"sv, AdvectionType::SemiLagrangian}, std::pair{"maccormack"sv, AdvectionType::MacCormack}});
        else if (key == "dimension") {
            const int dimension = number<int>(value, key, line);
            if (dimension != 2 && dimension != 3) throw std::invalid_argument("dimension must be 2 or 3");
            config.dimension = dimension == 2 ? Dimension::D2 : Dimension::D3;
        } else if (key == "fluid_solver") config.fluidSolver = choice(value, key, line,
            std::array{std::pair{"simple"sv, FluidSolverType::Simple}, std::pair{"reflection"sv, FluidSolverType::Reflection}});
        else if (key == "projection_solver") config.projectionSolver = choice(value, key, line,
            std::array{std::pair{"cpu_sor"sv, ProjectionSolverType::CpuSor}, std::pair{"gpu_sor"sv, ProjectionSolverType::GpuSor}});
        else if (key == "advection_backend") config.advectionBackend = choice(value, key, line,
            std::array{std::pair{"cpu"sv, ComputeBackend::Cpu}, std::pair{"gpu"sv, ComputeBackend::Gpu}});
        else if (key == "resolution") config.resolution = number<int>(value, key, line);
        else if (key == "box_size") config.boxSize = number<float>(value, key, line);
        else if (key == "emitter_radius") config.emitterRadius = number<float>(value, key, line);
        else if (key == "source_strength") config.sourceStrength = number<float>(value, key, line);
        else if (key == "buoyancy") config.buoyancy = number<float>(value, key, line);
        else if (key == "smoke_decay") config.smokeDecay = number<float>(value, key, line);
        else if (key == "sor_iterations") config.sorIterations = number<int>(value, key, line);
        else if (key == "sor_omega") config.sorOmega = number<float>(value, key, line);
        else if (key == "max_frames") config.maxFrames = number<int>(value, key, line);
        else throw std::invalid_argument("Unknown configuration key '" + std::string(key) + "' on line " + std::to_string(line));
    }
    static constexpr std::array required{
        "domain"sv, "advection"sv, "dimension"sv, "fluid_solver"sv,
        "projection_solver"sv, "advection_backend"sv, "resolution"sv,
        "box_size"sv, "emitter_radius"sv, "source_strength"sv, "buoyancy"sv, "smoke_decay"sv,
        "sor_iterations"sv, "sor_omega"sv, "max_frames"sv};
    for (const auto key : required)
        if (!seen.contains(std::string(key))) throw std::invalid_argument("Missing configuration key '" + std::string(key) + "'");
    validateSimulationConfig(config);
    return config;
}

SimulationConfig loadSimulationConfig(const std::filesystem::path& path) {
    std::ifstream file(path);
    if (!file) throw std::runtime_error("Could not open simulation configuration: " + path.string());
    return parseSimulationConfig(file);
}

void validateSimulationConfig(const SimulationConfig& c) {
    if (c.dimension != Dimension::D2 && c.dimension != Dimension::D3) throw std::invalid_argument("dimension must be 2 or 3");
    if (c.resolution <= 0 || c.resolution > 512) throw std::invalid_argument("resolution must be between 1 and 512");
    if (!std::isfinite(c.boxSize) || c.boxSize <= 0.f) throw std::invalid_argument("box_size must be finite and positive");
    if (!std::isfinite(c.emitterRadius) || c.emitterRadius <= 0.f || c.emitterRadius > 1.f) throw std::invalid_argument("emitter_radius must be in (0, 1]");
    if (!std::isfinite(c.sourceStrength) || c.sourceStrength < 0.f) throw std::invalid_argument("source_strength must be finite and nonnegative");
    if (!std::isfinite(c.buoyancy)) throw std::invalid_argument("buoyancy must be finite");
    if (!std::isfinite(c.smokeDecay) || c.smokeDecay < 0.f) throw std::invalid_argument("smoke_decay must be finite and nonnegative");
    if (c.sorIterations < 0 || c.sorIterations > 100000) throw std::invalid_argument("sor_iterations must be between 0 and 100000");
    if (!std::isfinite(c.sorOmega) || c.sorOmega <= 0.f || c.sorOmega >= 2.f) throw std::invalid_argument("sor_omega must be between 0 and 2");
    if (c.maxFrames < 0) throw std::invalid_argument("max_frames must be nonnegative");
}
