#pragma once

#include <memory>

class MacSimulationRegistry;
class SimulationRegistry;

// Registers MAC as one solver family through SimulationRegistry's public API.
// Other domain families use the same application-facing registration point.
void registerMacSimulationFactories(
    SimulationRegistry& registry,
    std::shared_ptr<MacSimulationRegistry> modules,
    bool include2D,
    bool include3D);

void registerBuiltinMacSimulationFactories(SimulationRegistry& registry);
