#include "composition/MacSimulationRegistry.hpp"

#include "fluid/advection/MacCormack.hpp"
#include "fluid/advection/SemiLagrangian.hpp"
#include "fluid/advection/SemiLagrangian2D.hpp"
#include "fluid/projection/CpuPressureSolver.hpp"
#include "fluid/projection/CpuPressureSolver2D.hpp"

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

std::string_view name(AdvectionType value) {
    switch (value) {
    case AdvectionType::SemiLagrangian: return "semi_lagrangian";
    case AdvectionType::MacCormack: return "maccormack";
    }
    return "unknown";
}

std::string_view name(ComputeBackend value) {
    switch (value) {
    case ComputeBackend::Cpu: return "cpu";
    case ComputeBackend::Gpu: return "gpu";
    }
    return "unknown";
}

std::string_view name(ProjectionSolverType value) {
    switch (value) {
    case ProjectionSolverType::CpuSor: return "cpu_sor";
    case ProjectionSolverType::GpuSor: return "gpu_sor";
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

[[noreturn]] void missingDomain(int dimension, DomainType domain) {
    throw std::invalid_argument("No " + std::to_string(dimension) +
                                "D domain factory registered for domain=" + std::string(name(domain)));
}

[[noreturn]] void missingAdvection(int dimension, const MacAdvectionKey& key) {
    throw std::invalid_argument("No " + std::to_string(dimension) +
        "D advection factory registered for domain=" + std::string(name(key.domain)) +
        ", advection=" + std::string(name(key.method)) +
        ", backend=" + std::string(name(key.backend)));
}

[[noreturn]] void missingPressure(int dimension, const MacProjectionKey& key) {
    throw std::invalid_argument("No " + std::to_string(dimension) +
        "D pressure factory registered for domain=" + std::string(name(key.domain)) +
        ", projection_solver=" + std::string(name(key.method)));
}

[[noreturn]] void missingFluid(int dimension, const MacFluidKey& key) {
    throw std::invalid_argument("No " + std::to_string(dimension) +
        "D fluid algorithm factory registered for domain=" + std::string(name(key.domain)) +
        ", fluid_solver=" + std::string(name(key.method)));
}

void configure(MacGridFluidSolver2D& solver, const SimulationConfig& config) {
    auto& parameters = solver.parameters();
    parameters.emitterRadius = config.emitterRadius;
    parameters.sourceStrength = config.sourceStrength;
    parameters.buoyancy = config.buoyancy;
    parameters.smokeDecay = config.smokeDecay;
    parameters.projectIterations = config.sorIterations;
    parameters.sorOmega = config.sorOmega;
}

void configure(MacGridFluidSolver& solver, const SimulationConfig& config) {
    auto& parameters = solver.parameters();
    parameters.macCormackSmoke = parameters.macCormackVel =
        config.advection == AdvectionType::MacCormack;
    parameters.emitterRadius = config.emitterRadius;
    parameters.sourceStrength = config.sourceStrength;
    parameters.buoyancy = config.buoyancy;
    parameters.smokeDecay = config.smokeDecay;
    parameters.projectIterations = config.sorIterations;
    parameters.sorOmega = config.sorOmega;
}

template<class Factory>
void requireFactory(const Factory& factory, const char* kind) {
    if (!factory) throw std::invalid_argument(std::string(kind) + " factory cannot be empty");
}
}

std::size_t MacSimulationRegistry::AdvectionKeyHash::operator()(const MacAdvectionKey& key) const noexcept {
    return (static_cast<std::size_t>(key.domain) << 16) ^
           (static_cast<std::size_t>(key.method) << 8) ^
           static_cast<std::size_t>(key.backend);
}

std::size_t MacSimulationRegistry::ProjectionKeyHash::operator()(const MacProjectionKey& key) const noexcept {
    return (static_cast<std::size_t>(key.domain) << 8) ^ static_cast<std::size_t>(key.method);
}

std::size_t MacSimulationRegistry::FluidKeyHash::operator()(const MacFluidKey& key) const noexcept {
    return (static_cast<std::size_t>(key.domain) << 8) ^ static_cast<std::size_t>(key.method);
}

MacSimulationRegistry::MacSimulationRegistry() {
    registerDomain2D(DomainType::Mac, [](const SimulationConfig& config) {
        return MacGridDomain2D(config.resolution, config.boxSize);
    });
    registerAdvection2D({DomainType::Mac, AdvectionType::SemiLagrangian, ComputeBackend::Cpu}, [] {
        return std::make_unique<SemiLagrangian2D>();
    });
    registerPressure2D({DomainType::Mac, ProjectionSolverType::CpuSor}, [] {
        return std::make_unique<CpuPressureSolver2D>();
    });
    registerFluid2D({DomainType::Mac, FluidSolverType::Simple},
        [](MacGridState2D state, const SimulationConfig&) {
            return std::make_unique<MacGridFluidSolver2D>(std::move(state));
        });

    registerDomain3D(DomainType::Mac, [](const SimulationConfig& config) {
        return MacGridDomain(config.resolution, config.boxSize);
    });
    registerAdvection3D({DomainType::Mac, AdvectionType::SemiLagrangian, ComputeBackend::Cpu}, [] {
        return std::make_unique<SemiLagrangian>();
    });
    registerAdvection3D({DomainType::Mac, AdvectionType::MacCormack, ComputeBackend::Cpu}, [] {
        return std::make_unique<MacCormack>();
    });
    registerPressure3D({DomainType::Mac, ProjectionSolverType::CpuSor}, [] {
        return std::make_unique<CpuPressureSolver>(6.f);
    });
    registerFluid3D({DomainType::Mac, FluidSolverType::Simple},
        [](MacGridState state, const SimulationConfig&) {
            auto solver = std::make_unique<MacGridFluidSolver>(std::move(state));
            solver->parameters().reflection = false;
            return solver;
        });
    registerFluid3D({DomainType::Mac, FluidSolverType::Reflection},
        [](MacGridState state, const SimulationConfig&) {
            auto solver = std::make_unique<MacGridFluidSolver>(std::move(state));
            solver->parameters().reflection = true;
            return solver;
        });
}

void MacSimulationRegistry::registerDomain2D(DomainType type, Domain2DFactory factory) {
    requireFactory(factory, "2D domain");
    domains2D_.insert_or_assign(type, std::move(factory));
}

void MacSimulationRegistry::registerAdvection2D(MacAdvectionKey key, Advection2DFactory factory) {
    requireFactory(factory, "2D advection");
    advections2D_.insert_or_assign(key, std::move(factory));
}

void MacSimulationRegistry::registerPressure2D(MacProjectionKey key, Pressure2DFactory factory) {
    requireFactory(factory, "2D pressure");
    pressures2D_.insert_or_assign(key, std::move(factory));
}

void MacSimulationRegistry::registerFluid2D(MacFluidKey key, Fluid2DFactory factory) {
    requireFactory(factory, "2D fluid");
    fluids2D_.insert_or_assign(key, std::move(factory));
}

void MacSimulationRegistry::registerDomain3D(DomainType type, Domain3DFactory factory) {
    requireFactory(factory, "3D domain");
    domains3D_.insert_or_assign(type, std::move(factory));
}

void MacSimulationRegistry::registerAdvection3D(MacAdvectionKey key, Advection3DFactory factory) {
    requireFactory(factory, "3D advection");
    advections3D_.insert_or_assign(key, std::move(factory));
}

void MacSimulationRegistry::registerPressure3D(MacProjectionKey key, Pressure3DFactory factory) {
    requireFactory(factory, "3D pressure");
    pressures3D_.insert_or_assign(key, std::move(factory));
}

void MacSimulationRegistry::registerFluid3D(MacFluidKey key, Fluid3DFactory factory) {
    requireFactory(factory, "3D fluid");
    fluids3D_.insert_or_assign(key, std::move(factory));
}

void MacSimulationRegistry::registerPrepare3D(Prepare3D prepare) {
    requireFactory(prepare, "3D preparation");
    prepare3D_ = std::move(prepare);
}

MacSimulation2D MacSimulationRegistry::create2D(const SimulationConfig& config) const {
    validateSimulationConfig(config);
    if (config.dimension != Dimension::D2)
        throw std::invalid_argument("A 2D registry cannot compose a non-2D configuration");

    const auto domainFactory = domains2D_.find(config.domain);
    if (domainFactory == domains2D_.end()) missingDomain(2, config.domain);

    const MacAdvectionKey advectionKey{config.domain, config.advection, config.advectionBackend};
    const auto advectionFactory = advections2D_.find(advectionKey);
    if (advectionFactory == advections2D_.end()) missingAdvection(2, advectionKey);

    const MacProjectionKey projectionKey{config.domain, config.projectionSolver};
    const auto pressureFactory = pressures2D_.find(projectionKey);
    if (pressureFactory == pressures2D_.end()) missingPressure(2, projectionKey);

    const MacFluidKey fluidKey{config.domain, config.fluidSolver};
    const auto fluidFactory = fluids2D_.find(fluidKey);
    if (fluidFactory == fluids2D_.end()) missingFluid(2, fluidKey);

    auto domain = domainFactory->second(config);
    auto advection = advectionFactory->second();
    auto pressure = pressureFactory->second();
    auto solver = fluidFactory->second(MacGridState2D(std::move(domain)), config);
    solver->setAdvection(*advection);
    solver->setPressureSolver(*pressure);
    configure(*solver, config);
    return {std::move(advection), std::move(pressure), std::move(solver)};
}

MacSimulation3D MacSimulationRegistry::create3D(const SimulationConfig& config) const {
    validateSimulationConfig(config);
    if (config.dimension != Dimension::D3)
        throw std::invalid_argument("A 3D registry cannot compose a non-3D configuration");

    const auto domainFactory = domains3D_.find(config.domain);
    if (domainFactory == domains3D_.end()) missingDomain(3, config.domain);

    const MacAdvectionKey advectionKey{config.domain, config.advection, config.advectionBackend};
    const auto advectionFactory = advections3D_.find(advectionKey);
    if (advectionFactory == advections3D_.end()) missingAdvection(3, advectionKey);

    const MacProjectionKey projectionKey{config.domain, config.projectionSolver};
    const auto pressureFactory = pressures3D_.find(projectionKey);
    if (pressureFactory == pressures3D_.end()) missingPressure(3, projectionKey);

    const MacFluidKey fluidKey{config.domain, config.fluidSolver};
    const auto fluidFactory = fluids3D_.find(fluidKey);
    if (fluidFactory == fluids3D_.end()) missingFluid(3, fluidKey);

    auto domain = domainFactory->second(config);
    auto advection = advectionFactory->second();
    auto pressure = pressureFactory->second();
    auto solver = fluidFactory->second(MacGridState(std::move(domain)), config);
    solver->setVelocityAdvection(*advection);
    solver->setDensityAdvection(*advection);
    solver->setPressureSolver(*pressure);
    configure(*solver, config);
    const bool needsPreparation = config.projectionSolver == ProjectionSolverType::GpuSor ||
                                  config.advectionBackend == ComputeBackend::Gpu;
    if (needsPreparation && prepare3D_ && !prepare3D_(config, solver->state())) {
        SimulationConfig fallback = config;
        if (fallback.projectionSolver == ProjectionSolverType::GpuSor)
            fallback.projectionSolver = ProjectionSolverType::CpuSor;
        if (fallback.advectionBackend == ComputeBackend::Gpu)
            fallback.advectionBackend = ComputeBackend::Cpu;
        return create3D(fallback);
    }
    return {config, std::move(advection), std::move(pressure), std::move(solver)};
}
