#include "composition/SimulationRegistry.hpp"

#include "composition/MacSimulationFactory.hpp"

#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>

namespace {
std::string_view name(DomainType value) {
    switch (value) {
    case DomainType::Mac: return "mac";
    case DomainType::Collocated: return "collocated";
    }
    return "unknown";
}

std::string_view name(FluidSolverType value) {
    switch (value) {
    case FluidSolverType::Simple: return "simple";
    case FluidSolverType::Reflection: return "reflection";
    }
    return "unknown";
}
}

std::size_t SimulationRegistry::KeyHash::operator()(const SimulationFactoryKey& key) const noexcept {
    return (static_cast<std::size_t>(key.dimension) << 16) ^
           (static_cast<std::size_t>(key.domain) << 8) ^
           static_cast<std::size_t>(key.fluidSolver);
}

SimulationRegistry::SimulationRegistry() {
    registerBuiltinMacSimulationFactories(*this);
}

SimulationRegistry::~SimulationRegistry() = default;

void SimulationRegistry::registerFactory(SimulationFactoryKey key, Factory factory) {
    if (!factory) throw std::invalid_argument("Simulation factory cannot be empty");
    factories_.insert_or_assign(key, std::move(factory));
}

SimulationInstance SimulationRegistry::create(const SimulationConfig& config) const {
    validateSimulationConfig(config);
    const SimulationFactoryKey key{config.dimension, config.domain, config.fluidSolver};
    const auto factory = factories_.find(key);
    if (factory == factories_.end()) {
        const int dimension = config.dimension == Dimension::D2 ? 2 : 3;
        throw std::invalid_argument(
            "No simulation factory registered for dimension=" + std::to_string(dimension) +
            ", domain=" + std::string(name(config.domain)) +
            ", fluid_solver=" + std::string(name(config.fluidSolver)));
    }
    return factory->second(config);
}
