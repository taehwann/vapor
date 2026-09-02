#pragma once

class SimulationRegistry;

// Adds GPU-capable solver factories without exposing a domain-specific backend
// to the application.
void registerGpuSimulationModules(SimulationRegistry& registry);
