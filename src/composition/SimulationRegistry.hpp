#pragma once

#include "composition/SimulationInstance.hpp"

#include <functional>
#include <unordered_map>

struct SimulationFactoryKey {
    Dimension dimension;
    DomainType domain;
    FluidSolverType fluidSolver;
    bool operator==(const SimulationFactoryKey&) const noexcept = default;
};

// Application-facing composition root. Its values are factories returning an
// IFluidSolver plus the matching presentation adapter; concrete domain families
// remain below this boundary.
class SimulationRegistry {
public:
    using Factory = std::function<SimulationInstance(const SimulationConfig&)>;

    SimulationRegistry();
    ~SimulationRegistry();
    SimulationRegistry(const SimulationRegistry&) = delete;
    SimulationRegistry& operator=(const SimulationRegistry&) = delete;

    void registerFactory(SimulationFactoryKey key, Factory factory);
    [[nodiscard]] SimulationInstance create(const SimulationConfig& config) const;

private:
    struct KeyHash {
        std::size_t operator()(const SimulationFactoryKey& key) const noexcept;
    };

    std::unordered_map<SimulationFactoryKey, Factory, KeyHash> factories_;
};
