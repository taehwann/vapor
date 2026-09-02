#include "composition/GpuMacModuleRegistration.hpp"

#include "composition/MacSimulationFactory.hpp"
#include "composition/MacSimulationRegistry.hpp"
#include "composition/SimulationRegistry.hpp"
#include "gpu/GpuAdvection.hpp"
#include "gpu/GpuPressureSolver.hpp"

#include <cstdio>
#include <memory>
#include <stdexcept>
#include <utility>

void registerGpuMacModules(MacSimulationRegistry& registry,
                           std::shared_ptr<GpuMacBackend> backend) {
    if (!backend) throw std::invalid_argument("GPU module registration requires a backend");
    registry.registerAdvection3D(
        {DomainType::Mac, AdvectionType::SemiLagrangian, ComputeBackend::Gpu},
        [backend] { return std::make_unique<GpuAdvection>(*backend, false); });
    registry.registerAdvection3D(
        {DomainType::Mac, AdvectionType::MacCormack, ComputeBackend::Gpu},
        [backend] { return std::make_unique<GpuAdvection>(*backend, true); });
    registry.registerPressure3D(
        {DomainType::Mac, ProjectionSolverType::GpuSor},
        [backend] { return std::make_unique<GpuPressureSolver>(*backend); });
    registry.registerPrepare3D([backend](const SimulationConfig&, const MacGridState& state) {
        try {
            backend->init(state);
            return true;
        } catch (const std::exception& error) {
            std::fprintf(stderr, "Configured GPU unavailable; using CPU: %s\n", error.what());
            return false;
        }
    });
}

void registerGpuSimulationModules(SimulationRegistry& registry) {
    auto modules = std::make_shared<MacSimulationRegistry>();
    registerGpuMacModules(*modules, std::make_shared<GpuMacBackend>());
    registerMacSimulationFactories(registry, std::move(modules), false, true);
}
