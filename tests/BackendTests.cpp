#include "solvers/mac3d/MacGridFluidSolver3D.hpp"
#include "renderer/MacRenderData.hpp"

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

class RecordingPressure final : public Projection3D {
public:
    explicit RecordingPressure(const MacGridState3D& expectedState) : expected_(&expectedState) {}
    ProjectionResult project(MacGridState3D& state, const ProjectionOptions& options) override {
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
    const MacGridState3D* expected_;
    CpuPressureSolver3D cpu_{6.f};
};

struct TransportCalls {
    int faces = 0, scalars = 0, destroyed = 0;
    float faceTime = 0.f, scalarTime = 0.f;
    AdvectionMethod velocity = AdvectionMethod::SemiLagrangian;
    AdvectionMethod density = AdvectionMethod::SemiLagrangian;
};
class RecordingAdvection final : public Advection3D {
public:
    explicit RecordingAdvection(std::shared_ptr<TransportCalls> calls) : calls_(std::move(calls)) {}
    ~RecordingAdvection() override { ++calls_->destroyed; }
    void advectScalar(const AdvectionContext& context, std::span<float> field,
                      AdvectionWorkspace& workspace, AdvectionMethod method) const override {
        ++calls_->scalars; calls_->scalarTime += context.dt; calls_->density = method;
        cpu_.advectScalar(context, field, workspace, method);
    }
    void advectFace(const AdvectionContext& context, FaceAxis axis, std::span<const float> source,
                    std::span<float> output, AdvectionWorkspace& workspace, AdvectionMethod method) const override {
        require(source.data() != output.data(), "Velocity source/output alias");
        ++calls_->faces; calls_->faceTime += context.dt; calls_->velocity = method;
        cpu_.advectFace(context, axis, source, output, workspace, method);
    }
private:
    std::shared_ptr<TransportCalls> calls_;
    CpuAdvection cpu_;
};

void interfaceWiring(bool reflect) {
    MacGridFluidSolver3D solver(8), reference(8);
    for (auto* item : {&solver, &reference}) {
        auto& options = item->parameters();
        options.stirStrength = 0.f;
        options.emitterRadius = .18f;
        options.emitterCenterZ = .3f;
        options.reflection = reflect;
        options.projectIterations = 16;
        options.macCormackVel = options.macCormackSmoke = false;
    }
    auto pressureOwner = std::make_unique<RecordingPressure>(solver.state());
    auto* pressure = pressureOwner.get();
    auto calls = std::make_shared<TransportCalls>();
    solver.setPressureSolver(std::move(pressureOwner));
    solver.setAdvection(std::make_unique<RecordingAdvection>(calls));

    auto& fluid = solver;

    const auto& domain = fluid.domain();
    require(domain.dimension() == Dimension::D3, "Domain dimension");
    near(domain.bounds().min.x, 0.f);
    near(domain.bounds().max.z, 4.5f);
    constexpr float dt = 1.f / 60.f;
    for (int frame = 0; frame < 3; ++frame) {
        fluid.advance(dt);
        reference.advance(dt);
    }
    const int stages = reflect ? 2 : 1;
    require(pressure->steps.size() == 3 * stages, "Pressure stage count");
    require(calls->faces == 9 * stages, "Selected velocity backend");
    require(calls->scalars == 3, "Selected density backend");
    near(calls->faceTime, 9 * dt);
    near(calls->scalarTime, 3 * dt);
    for (float value : pressure->steps) near(value, dt / stages);
    for (int value : pressure->iterations) require(value == 16 / stages, "Pressure work budget");
    require(solver.lastProjection().residualAvailable, "CPU residual must be marked available");
    const auto output = volumeData(solver.state());
    const auto expected = reference.state().density();
    require(output.density.data() == solver.state().density().data(), "RenderData must borrow physical density");
    require(output.width == 8 && output.height == 8 && output.depth == 8, "Volume dimensions must match MAC state");
    require(output.density.size() == expected.size(), "Volume density extent must match MAC state");
    near(output.boxSize, domain.bounds().max.x);
    for (std::size_t i = 0; i < output.density.size(); ++i) near(output.density[i], expected[i]);
    near(fluid.diagnostics().kineticEnergy, reference.kineticEnergy());
    near(fluid.diagnostics().maxDivergence, reference.computeDivergenceNorm());

    // Flags must still select methods after installing a backend.
    solver.parameters().macCormackVel = reference.parameters().macCormackVel = true;
    fluid.advance(dt); reference.advance(dt);
    require(calls->velocity == AdvectionMethod::MacCormack &&
            calls->density == AdvectionMethod::SemiLagrangian, "Per-field method settings ignored");
    near(fluid.diagnostics().kineticEnergy, reference.kineticEnergy());
    fluid.resetState();
    require(pressure->resets == 1, "Reset must reach owned pressure backend");
    for (float value : volumeData(solver.state()).density) near(value, 0.f);
    solver.useDefaultPressureSolver(); // Destroys the old backend.
    solver.useDefaultAdvection();
    require(calls->destroyed == 1, "Detached advection must be destroyed by the solver");
    const auto previousCalls = calls->faces;
    fluid.advance(dt);
    require(calls->faces == previousCalls, "Detached backend still called");
    bool rejected = false;
    try { solver.setAdvection(nullptr); } catch (const std::invalid_argument&) { rejected = true; }
    require(rejected, "Null backend accepted");
    solver.setBoxSize(3.f);
    near(domain.bounds().max.x, 3.f);
    near(volumeData(solver.state()).boxSize, 3.f);
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
