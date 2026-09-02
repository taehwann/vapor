#include "fluid/MacGridOperators.hpp"
#include "fluid/MacGridFluidSolver.hpp"
#include "gpu/GpuPressureSolver.hpp"
#include "gpu/GpuAdvection.hpp"
#include "graphics/GlLoader.hpp"
#include "presentation/MacGridPresentationAdapter.hpp"
#include "renderer/OrbitCamera.hpp"
#include "renderer/Renderer.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {
void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}
void constant(std::span<const float> field, float expected) {
    for (float value : field)
        require(std::isfinite(value) && std::abs(value-expected)<1e-5f, "GPU constant-field advection failed");
}
template<class F> void rejects(F function) {
    try { function(); }
    catch (const std::invalid_argument&) { return; }
    throw std::runtime_error("GPU interface accepted invalid input");
}

void advectionParity(GpuMacBackend& backend) {
    MacGridState state(8, 4.f);
    std::fill(state.velocityX().begin(), state.velocityX().end(), .2f);
    std::fill(state.velocityY().begin(), state.velocityY().end(), -.1f);
    std::fill(state.velocityZ().begin(), state.velocityZ().end(), .3f);
    const AdvectionContext context{state, velocityView(state), .1f, 3.f};
    const auto count = state.velocityX().size();
    std::vector<float> source(count), cpuOut(count), gpuOut(count);
    for (std::size_t i = 0; i < count; ++i) source[i] = .2f + .1f * std::sin(float(i) * .13f);
    SemiLagrangian sl;
    MacCormack mc;
    AdvectionWorkspace cpuScratch, gpuScratch;
    for (bool correction : {false, true}) {
        const MacGridAdvection& cpu = correction ? static_cast<const MacGridAdvection&>(mc) : sl;
        GpuAdvection gpu(backend, correction);
        for (auto axis : {FaceAxis::X, FaceAxis::Y, FaceAxis::Z}) {
            cpu.advectFace(context, axis, source, cpuOut, cpuScratch);
            gpu.advectFace(context, axis, source, gpuOut, gpuScratch);
            for (std::size_t i = 0; i < count; ++i) {
                require(std::isfinite(gpuOut[i]) && std::abs(cpuOut[i] - gpuOut[i]) < 2e-5f,
                        "CPU/GPU face advection differs after interface dispatch");
            }
        }
        rejects([&] { gpu.advectFace(context, FaceAxis::X, source, source, gpuScratch); });
        rejects([&] { gpu.advectFace(context, FaceAxis::X, source, state.velocityX(), gpuScratch); });
        rejects([&] { gpu.advectScalar(context, std::span<float>(gpuOut).first(1), gpuScratch); });
    }
    // Scalar MacCormack is intentionally not compared here: the preserved CPU
    // limiter updates its source in place, whereas the GPU reads a saved source.
    std::vector<float> cpuDensity(source.begin(), source.begin() + state.cellCount());
    auto gpuDensity = cpuDensity;
    GpuAdvection gpuSl(backend, false);
    sl.advectScalar(context, cpuDensity, cpuScratch);
    gpuSl.advectScalar(context, gpuDensity, gpuScratch);
    for (std::size_t i = 0; i < cpuDensity.size(); ++i)
        require(std::abs(cpuDensity[i] - gpuDensity[i]) < 2e-5f, "CPU/GPU scalar transport differs");
}
}

int main() {
    if (!glfwInit()) { std::cout << "SKIP: GLFW unavailable\n"; return 77; }
    glfwWindowHint(GLFW_VISIBLE, GLFW_FALSE);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    GLFWwindow* window=glfwCreateWindow(64,64,"Vapor hidden graphics check",nullptr,nullptr);
    if (!window) { glfwTerminate(); std::cout << "SKIP: OpenGL 4.3 unavailable\n"; return 77; }
    glfwMakeContextCurrent(window);
    if (!loadOpenGLFunctions()) { glfwDestroyWindow(window); glfwTerminate(); return 1; }
    int status=0;
    GpuMacBackend gpu;
    Renderer renderer;
    try {
        // Loading GL entry points here must also initialize the functions used
        // by the separate GPU and renderer translation units.
        MacGridState state(8,4.f);
        const auto count=static_cast<int>(state.velocityX().size());
        gpu.init(state);
        require(gpu.initialized(),"GPU initialization failed");
        GpuPressureSolver gpuPressure(gpu);
        IPressureSolver<MacGridState>& pressure = gpuPressure;
        MacGridState wrongSize(4);
        rejects([&] { pressure.project(wrongSize, {}); });
        state.velocityX()[state.idX(4,4,4)]=1.f;
        const float before=MacGridOperators::maxDivergence(state);
        ProjectionOptions projectionOptions;
        projectionOptions.relaxation = 1.5f;
        projectionOptions.linearSolve.maxIterations = 200;
        projectionOptions.linearSolve.fixedIterations = true;
        const auto projection = pressure.project(state, projectionOptions);
        require(!projection.residualAvailable, "GPU must not report an unmeasured residual");
        require(MacGridOperators::maxDivergence(state)<before*.01f,"GPU projection did not reduce divergence");

        state.reset();
        std::fill(state.density().begin(),state.density().end(),.5f);
        std::vector<float> source(count,.125f),output(count);
        AdvectionWorkspace workspace;
        for (bool correction : {false,true}) {
            GpuAdvection gpuMethod(gpu, correction);
            const IAdvection<MacAdvectionOperation>& method = gpuMethod;
            const AdvectionContext context{state, velocityView(state), .02f, 3.f};
            method.advect({context, state.density(), state.density(), workspace, std::nullopt});
            constant(state.density(),.5f);
            for (auto axis : {FaceAxis::X, FaceAxis::Y, FaceAxis::Z}) {
                method.advect({context, source, output, workspace, axis});
                constant(output,.125f);
            }
        }
        require(glGetError()==GL_NO_ERROR,"OpenGL error in compute backend");
        advectionParity(gpu);

        // Exercise the coordinator through all CPU/GPU pressure/advection
        // combinations, including switching backends on a live state.
        GpuAdvection gpuSl(gpu, false), gpuMc(gpu, true);
        MacGridFluidSolver simulation(8);
        simulation.parameters().stirStrength = 0.f;
        simulation.parameters().emitterRadius = .18f;
        simulation.parameters().emitterCenterZ = .3f;
        simulation.parameters().projectIterations = 16;
        IFluidSolver& fluid = simulation;
        for (bool reflect : {false, true}) for (bool correction : {false, true}) {
            simulation.parameters().reflection = reflect;
            simulation.parameters().macCormackVel = correction;
            simulation.parameters().macCormackSmoke = correction;
            for (bool gpuP : {false, true}) for (bool gpuA : {false, true}) {
                if (gpuP) simulation.setPressureSolver(gpuPressure);
                else simulation.useDefaultPressureSolver();
                if (gpuA) {
                    const auto& method = correction ? gpuMc : gpuSl;
                    simulation.setVelocityAdvection(method);
                    simulation.setDensityAdvection(method);
                } else simulation.useDefaultAdvection();
                fluid.advance(1.f / 60.f);
                const auto stats = fluid.diagnostics();
                require(std::isfinite(stats.kineticEnergy) && std::isfinite(stats.maxDivergence),
                        "Backend switch produced invalid state");
                require(simulation.lastProjection().residualAvailable != gpuP,
                        "Pressure diagnostics belong to the wrong backend");
            }
        }
        require(glGetError()==GL_NO_ERROR,"OpenGL error after backend switching");

        MacGridPresentationAdapter presentation(state);
        const RenderData data = presentation.renderData();
        renderer.init(data);
        renderer.updateDomain(data);
        glClearColor(0.f,0.f,0.f,1.f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        float mvp[16],cx,cy,cz;
        buildMVP(mvp,cx,cy,cz,-1.2f,.3f,10.f,1.f,state.boxSize());
        renderer.render(data,64,64,mvp,cx,cy,cz,.5f,30.f,0.f,0.f,1.f,0.f,.1f,nullptr);
        std::vector<unsigned char> pixels(64*64*3);
        glReadBuffer(GL_BACK);
        glReadPixels(0,0,64,64,GL_RGB,GL_UNSIGNED_BYTE,pixels.data());
        require(glGetError()==GL_NO_ERROR,"OpenGL error in renderer");
        require(std::any_of(pixels.begin(),pixels.end(),[](unsigned char v){return v>32;}),
                "Renderer produced no visible volume samples");
        std::cout << "PASS hidden OpenGL compute and volume-rendering smoke test\n";
    } catch (const std::exception& error) {
        std::cerr << "FAIL graphics smoke test: " << error.what() << '\n'; status=1;
    }
    renderer.shutdown();
    gpu.shutdown();
    require(!gpu.initialized(), "GPU shutdown must clear readiness");
    gpu.shutdown(); // Explicit shutdown followed by destruction is safe.
    glfwDestroyWindow(window);
    glfwTerminate();
    return status;
}
