#include "composition/MacSimulationRegistry.hpp"
#include "composition/SimulationRegistry.hpp"
#include "fluid/advection/SemiLagrangian2D.hpp"

#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <variant>

namespace {
void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

template<class Function>
void rejectsWith(Function&& function, std::string_view expected) {
    try {
        function();
    } catch (const std::invalid_argument& error) {
        require(std::string_view(error.what()).find(expected) != std::string_view::npos,
                "Registry error did not identify the missing module");
        return;
    }
    throw std::runtime_error("Missing registry entry was accepted");
}

class CountingAdvection2D final : public IAdvection<MacAdvectionOperation2D> {
public:
    explicit CountingAdvection2D(std::shared_ptr<int> calls) : calls_(std::move(calls)) {}
    [[nodiscard]] std::string_view name() const noexcept override { return "registered test advection"; }
    void advect(const MacAdvectionOperation2D& operation) const override {
        ++*calls_;
        implementation_.advect(operation);
    }
private:
    std::shared_ptr<int> calls_;
    SemiLagrangian2D implementation_;
};

class TestCollocatedDomain final : public IDomain {
public:
    explicit TestCollocatedDomain(float size) : size_(size) {}
    [[nodiscard]] Dimension dimension() const noexcept override { return Dimension::D2; }
    [[nodiscard]] DomainBounds bounds() const noexcept override {
        return {{0.f, 0.f, 0.f}, {size_, size_, 0.f}};
    }
    void setSize(float size) { size_ = size; }
private:
    float size_;
};

class TestCollocatedSolver final : public IFluidSolver {
public:
    explicit TestCollocatedSolver(float size) : domain_(size) {}
    void advance(float) override { ++steps_; }
    void resetState() override { steps_ = 0; }
    [[nodiscard]] const IDomain& domain() const noexcept override { return domain_; }
    [[nodiscard]] FluidDiagnostics diagnostics() const override {
        return {static_cast<float>(steps_), 0.f, 0.f};
    }
    [[nodiscard]] FluidSolverControls controls() const noexcept override { return controls_; }
    void setControls(const FluidSolverControls& controls) override { controls_ = controls; }
    [[nodiscard]] PressureSolveDiagnostics pressureDiagnostics() const noexcept override { return {}; }
    [[nodiscard]] float boxSize() const noexcept override { return domain_.bounds().max.x; }
    void setBoxSize(float size) override { domain_.setSize(size); }
private:
    TestCollocatedDomain domain_;
    FluidSolverControls controls_;
    int steps_ = 0;
};

class TestCollocatedPresentation final : public IPresentationSource {
public:
    [[nodiscard]] Dimension dimension() const noexcept override { return Dimension::D2; }
    [[nodiscard]] PresentationData presentationData() const noexcept override {
        return ImageRenderData{};
    }
};
}

int main() {
    try {
        MacSimulationRegistry registry;
        SimulationConfig two;
        two.dimension = Dimension::D2;
        two.resolution = 12;
        two.boxSize = 2.5f;
        two.sourceStrength = 0.f;
        two.sorIterations = 20;
        auto modules2D = registry.create2D(two);
        require(modules2D.solver->domain().dimension() == Dimension::D2, "2D domain selection");
        require(modules2D.solver->state().resolution() == 12, "2D domain configuration");
        require(modules2D.solver->state().boxSize() == 2.5f, "2D box configuration");
        require(modules2D.advection->name() == "2D semi-Lagrangian", "2D advection selection");
        require(modules2D.solver->parameters().projectIterations == 20, "2D parameter assignment");

        auto calls = std::make_shared<int>(0);
        registry.registerAdvection2D(
            {DomainType::Mac, AdvectionType::SemiLagrangian, ComputeBackend::Cpu},
            [calls] { return std::make_unique<CountingAdvection2D>(calls); });
        auto overridden = registry.create2D(two);
        overridden.solver->advance(.01f);
        require(*calls == 3, "Registered 2D advection was not injected");

        auto unavailable = two;
        unavailable.advection = AdvectionType::MacCormack;
        rejectsWith([&] { (void)registry.create2D(unavailable); }, "advection=maccormack");
        unavailable = two;
        unavailable.domain = DomainType::Collocated;
        rejectsWith([&] { (void)registry.create2D(unavailable); }, "domain=collocated");
        unavailable = two;
        unavailable.projectionSolver = ProjectionSolverType::GpuSor;
        rejectsWith([&] { (void)registry.create2D(unavailable); }, "projection_solver=gpu_sor");

        SimulationConfig three = two;
        three.dimension = Dimension::D3;
        three.advection = AdvectionType::MacCormack;
        three.fluidSolver = FluidSolverType::Reflection;
        auto modules3D = registry.create3D(three);
        require(modules3D.solver->domain().dimension() == Dimension::D3, "3D domain selection");
        require(modules3D.advection->name() == "MacCormack", "3D advection selection");
        require(modules3D.solver->parameters().reflection, "3D fluid algorithm selection");
        require(modules3D.solver->parameters().projectIterations == 20, "3D parameter assignment");

        three.fluidSolver = FluidSolverType::Simple;
        auto simple = registry.create3D(three);
        require(!simple.solver->parameters().reflection, "Simple fluid algorithm selection");

        three.advectionBackend = ComputeBackend::Gpu;
        rejectsWith([&] { (void)registry.create3D(three); }, "backend=gpu");

        SimulationRegistry simulations;
        auto runtime2D = simulations.create(two);
        IFluidSolver& abstract2D = runtime2D.solver();
        require(abstract2D.domain().dimension() == Dimension::D2,
                "Application registry did not return a 2D IFluidSolver");
        require(std::holds_alternative<ImageRenderData>(
                    runtime2D.presentation().presentationData()),
                "2D presentation was not type-erased with its solver");
        auto commonControls = abstract2D.controls();
        commonControls.sourceStrength = 1.25f;
        abstract2D.setControls(commonControls);
        require(abstract2D.controls().sourceStrength == 1.25f,
                "IFluidSolver controls did not reach the concrete solver");

        three.advectionBackend = ComputeBackend::Cpu;
        three.fluidSolver = FluidSolverType::Reflection;
        auto runtime3D = simulations.create(three);
        require(runtime3D.solver().domain().dimension() == Dimension::D3,
                "Application registry did not return a 3D IFluidSolver");
        require(std::holds_alternative<RenderData>(
                    runtime3D.presentation().presentationData()),
                "3D presentation was not type-erased with its solver");

        unavailable = two;
        unavailable.domain = DomainType::Collocated;
        rejectsWith([&] { (void)simulations.create(unavailable); },
                    "No simulation factory registered");
        simulations.registerFactory(
            {Dimension::D2, DomainType::Collocated, FluidSolverType::Simple},
            [](const SimulationConfig& selected) {
                return SimulationInstance(
                    selected,
                    std::make_unique<TestCollocatedSolver>(selected.boxSize),
                    std::make_unique<TestCollocatedPresentation>());
            });
        auto collocated = simulations.create(unavailable);
        collocated.solver().advance(.01f);
        require(collocated.solver().diagnostics().kineticEnergy == 1.f,
                "Registered non-MAC IFluidSolver factory was not used");
        unavailable = two;
        unavailable.fluidSolver = FluidSolverType::Reflection;
        rejectsWith([&] { (void)simulations.create(unavailable); },
                    "fluid_solver=reflection");

        std::cout << "PASS IFluidSolver factory selection, module lookup, composition, and injection\n";
    } catch (const std::exception& error) {
        std::cerr << "FAIL composition: " << error.what() << '\n';
        return 1;
    }
}
