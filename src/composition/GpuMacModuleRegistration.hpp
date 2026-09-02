#pragma once

#include <memory>

class GpuMacBackend;
class MacSimulationRegistry;

// Adds GPU implementations and their shared OpenGL resource owner to a typed
// MAC registry.
void registerGpuMacModules(MacSimulationRegistry& registry,
                           std::shared_ptr<GpuMacBackend> backend);
