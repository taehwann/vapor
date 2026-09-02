#include "fluid/MacGridFluidSolver.hpp"
#include "presentation/MacGridPresentationAdapter.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {
void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}
void near(float actual, float expected) {
    require(std::isfinite(actual) && std::abs(actual - expected) < 2e-5f, "Interface changed simulation results");
}

class RecordingPressure final : public IPressureSolver<MacGridState> {
public:
    explicit RecordingPressure(const MacGridState& expectedState) : expected_(&expectedState) {}
    ProjectionResult project(MacGridState& state, const ProjectionOptions& options) override {
        require(&state == expected_, "Pressure must receive solver-owned state");
        steps.push_back(options.dt);
        iterations.push_back(options.linearSolve.maxIterations);
        return cpu_.project(state, options);
    }
    void reset() noexcept override { ++resets; cpu_.reset(); }
    std::vector<float> steps;
    std::vector<int> iterations;
    int resets = 0;
private:
    const MacGridState* expected_;
    CpuPressureSolver cpu_{6.f};
};

class RecordingAdvection final : public IAdvection<MacAdvectionOperation> {
public:
    std::string_view name() const noexcept override { return "Recording transport"; }
    void advect(const MacAdvectionOperation& operation) const override {
        if (operation.axis) {
            require(operation.source.data() != operation.output.data(), "Velocity source/output alias");
            ++faces;
            faceTime += operation.context.dt;
        } else {
            ++scalars;
            scalarTime += operation.context.dt;
        }
        method_.advect(operation);
    }
    mutable int faces = 0, scalars = 0;
    mutable float faceTime = 0.f, scalarTime = 0.f;
private:
    SemiLagrangian method_;
};

void interfaceWiring(bool reflect) {
    MacGridFluidSolver solver(8), reference(8);
    for (auto* item : {&solver, &reference}) {
        auto& options = item->parameters();
        options.stirStrength = 0.f;
        options.emitterRadius = .18f;
        options.emitterCenterZ = .3f;
        options.reflection = reflect;
        options.projectIterations = 16;
        options.macCormackVel = options.macCormackSmoke = false;
    }
    RecordingPressure pressure(solver.state());
    RecordingAdvection velocity, density;
    solver.setPressureSolver(pressure);
    solver.setVelocityAdvection(velocity);
    solver.setDensityAdvection(density);

    IFluidSolver& fluid = solver;
    MacGridPresentationAdapter presentation(solver.state());
    const IDomain& domain = fluid.domain();
    require(domain.dimension() == Dimension::D3, "Domain dimension");
    near(domain.bounds().min.x, 0.f);
    near(domain.bounds().max.z, 4.5f);
    constexpr float dt = 1.f / 60.f;
    for (int frame = 0; frame < 3; ++frame) {
        fluid.advance(dt);
        reference.advance(dt);
    }
    const int stages = reflect ? 2 : 1;
    require(pressure.steps.size() == 3 * stages, "Pressure stage count");
    require(velocity.faces == 9 * stages && velocity.scalars == 0, "Selected velocity backend");
    require(density.scalars == 3 && density.faces == 0, "Selected density backend");
    near(velocity.faceTime, 9 * dt);
    near(density.scalarTime, 3 * dt);
    for (float value : pressure.steps) near(value, dt / stages);
    for (int value : pressure.iterations) require(value == 16 / stages, "Pressure work budget");
    require(solver.lastProjection().residualAvailable, "CPU residual must be marked available");
    const auto output = presentation.renderData();
    const auto expected = reference.state().density();
    require(output.density.data() == solver.state().density().data(), "RenderData must borrow physical density");
    require(output.width == 8 && output.height == 8 && output.depth == 8, "Volume dimensions must match MAC state");
    require(output.density.size() == expected.size(), "Volume density extent must match MAC state");
    near(output.boxSize, domain.bounds().max.x);
    for (std::size_t i = 0; i < output.density.size(); ++i) near(output.density[i], expected[i]);
    near(fluid.diagnostics().kineticEnergy, reference.kineticEnergy());
    near(fluid.diagnostics().maxDivergence, reference.computeDivergenceNorm());

    // Switching implementations must preserve fields and stop calling the old backend.
    const auto previousCalls = pressure.steps.size();
    solver.useDefaultPressureSolver();
    solver.useDefaultAdvection();
    fluid.advance(dt);
    reference.advance(dt);
    require(pressure.steps.size() == previousCalls, "Detached pressure still called");
    near(fluid.diagnostics().kineticEnergy, reference.kineticEnergy());

    solver.setPressureSolver(pressure);
    fluid.resetState();
    require(pressure.resets == 1, "Reset must reach active pressure backend");
    near(fluid.diagnostics().kineticEnergy, 0.f);
    for (float value : presentation.renderData().density) near(value, 0.f);
    solver.setBoxSize(3.f);
    near(domain.bounds().max.x, 3.f);
    near(presentation.renderData().boxSize, 3.f);
}
}

int main() {
    try {
        interfaceWiring(false);
        interfaceWiring(true);
        std::cout << "PASS domain, fluid, pressure, advection, and presentation wiring\n";
    } catch (const std::exception& error) {
        std::cerr << "FAIL interfaces: " << error.what() << '\n';
        return 1;
    }
    return 0;
}
