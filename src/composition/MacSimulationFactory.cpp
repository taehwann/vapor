#include "composition/MacSimulationFactory.hpp"

#include "composition/MacSimulationRegistry.hpp"
#include "composition/SimulationInstance.hpp"
#include "composition/SimulationRegistry.hpp"
#include "presentation/MacGridPresentationAdapter.hpp"
#include "presentation/MacGridPresentationAdapter2D.hpp"

#include <memory>
#include <stdexcept>
#include <utility>

namespace {
class OwnedMacSolver2D final : public IFluidSolver {
public:
    explicit OwnedMacSolver2D(MacSimulation2D composition)
        : composition_(std::move(composition)) {}

    void advance(float dt) override { implementation().advance(dt); }
    void resetState() override { implementation().resetState(); }
    [[nodiscard]] const IDomain& domain() const noexcept override { return implementation().domain(); }
    [[nodiscard]] FluidDiagnostics diagnostics() const override { return implementation().diagnostics(); }
    [[nodiscard]] FluidSolverControls controls() const noexcept override { return implementation().controls(); }
    void setControls(const FluidSolverControls& controls) override { implementation().setControls(controls); }
    [[nodiscard]] PressureSolveDiagnostics pressureDiagnostics() const noexcept override {
        return implementation().pressureDiagnostics();
    }
    [[nodiscard]] float boxSize() const noexcept override { return implementation().boxSize(); }
    void setBoxSize(float size) override { implementation().setBoxSize(size); }
    [[nodiscard]] const MacGridState2D& state() const noexcept { return implementation().state(); }

private:
    [[nodiscard]] MacGridFluidSolver2D& implementation() noexcept { return *composition_.solver; }
    [[nodiscard]] const MacGridFluidSolver2D& implementation() const noexcept {
        return *composition_.solver;
    }
    MacSimulation2D composition_;
};

class OwnedMacSolver3D final : public IFluidSolver {
public:
    explicit OwnedMacSolver3D(MacSimulation3D composition)
        : composition_(std::move(composition)) {}

    void advance(float dt) override { implementation().advance(dt); }
    void resetState() override { implementation().resetState(); }
    [[nodiscard]] const IDomain& domain() const noexcept override { return implementation().domain(); }
    [[nodiscard]] FluidDiagnostics diagnostics() const override { return implementation().diagnostics(); }
    [[nodiscard]] FluidSolverControls controls() const noexcept override { return implementation().controls(); }
    void setControls(const FluidSolverControls& controls) override { implementation().setControls(controls); }
    [[nodiscard]] PressureSolveDiagnostics pressureDiagnostics() const noexcept override {
        return implementation().pressureDiagnostics();
    }
    [[nodiscard]] float boxSize() const noexcept override { return implementation().boxSize(); }
    void setBoxSize(float size) override { implementation().setBoxSize(size); }
    [[nodiscard]] const MacGridState& state() const noexcept { return implementation().state(); }
    [[nodiscard]] const SimulationConfig& configuration() const noexcept {
        return composition_.configuration;
    }

private:
    [[nodiscard]] MacGridFluidSolver& implementation() noexcept { return *composition_.solver; }
    [[nodiscard]] const MacGridFluidSolver& implementation() const noexcept {
        return *composition_.solver;
    }
    MacSimulation3D composition_;
};

SimulationInstance createMac2D(
    const std::shared_ptr<MacSimulationRegistry>& modules,
    const SimulationConfig& config) {
    auto solver = std::make_unique<OwnedMacSolver2D>(modules->create2D(config));
    auto presentation = std::make_unique<MacGridPresentationAdapter2D>(solver->state());
    return {config, std::move(solver), std::move(presentation)};
}

SimulationInstance createMac3D(
    const std::shared_ptr<MacSimulationRegistry>& modules,
    const SimulationConfig& config) {
    auto solver = std::make_unique<OwnedMacSolver3D>(modules->create3D(config));
    const SimulationConfig effectiveConfig = solver->configuration();
    auto presentation = std::make_unique<MacGridPresentationAdapter>(solver->state());
    return {effectiveConfig, std::move(solver), std::move(presentation)};
}
}

void registerMacSimulationFactories(
    SimulationRegistry& registry,
    std::shared_ptr<MacSimulationRegistry> modules,
    bool include2D,
    bool include3D) {
    if (!modules) throw std::invalid_argument("MAC factory registration requires module registries");
    if (include2D) {
        registry.registerFactory(
            {Dimension::D2, DomainType::Mac, FluidSolverType::Simple},
            [modules](const SimulationConfig& config) {
                return createMac2D(modules, config);
            });
    }
    if (include3D) {
        const auto factory = [modules](const SimulationConfig& config) {
            return createMac3D(modules, config);
        };
        registry.registerFactory(
            {Dimension::D3, DomainType::Mac, FluidSolverType::Simple}, factory);
        registry.registerFactory(
            {Dimension::D3, DomainType::Mac, FluidSolverType::Reflection}, factory);
    }
}

void registerBuiltinMacSimulationFactories(SimulationRegistry& registry) {
    registerMacSimulationFactories(
        registry, std::make_shared<MacSimulationRegistry>(), true, true);
}
