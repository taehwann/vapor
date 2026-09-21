#include "solvers/mac3d/MacGridOperators3D.hpp"
#include "solvers/mac3d/MacGridFluidSolver3D.hpp"
#include "solvers/mac3d/gpu/GpuPressureSolver.hpp"
#include "solvers/mac3d/gpu/GpuAdvection.hpp"
#include "solvers/mac3d/gpu/GpuMacSimulation3D.hpp"
#include "graphics/GlLoader.hpp"
#include "renderer/MacRenderData.hpp"
#include "renderer/OrbitCamera.hpp"
#include "renderer/VolumeRenderer.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {
void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}
void constant(std::span<const float> field, float expected);
void residentParity() {
    for (bool reflection : {false,true}) for (bool mcVelocity : {false,true}) for (bool mcSmoke : {false,true}) {
        MacGridFluidSolver3D reference(MacGridState3D(9,1.f));
        auto& p=reference.parameters();
        p.reflection=reflection; p.macCormackVel=mcVelocity; p.macCormackSmoke=mcSmoke;
        p.emitterRadius=.22f; p.emitterCenterZ=.3f; p.stirStrength=.5f; p.projectIterations=20;
        auto backend=std::make_shared<GpuMacBackend>(); backend->init(reference.state());
        reference.setAdvection(std::make_unique<GpuAdvection>(backend));
        reference.setPressureSolver(std::make_unique<GpuPressureSolver>(backend));
        GpuMacSimulation3D resident(9,1.f,p);
        for(int f=0;f<8;++f){reference.advance(1.f/60.f);resident.advance(1.f/60.f);}
        require(resident.downloadedBytes()==0,"Resident advance downloaded fields");
        require(resident.uploadedBytes()<8*reference.state().density().size_bytes(),"Resident advance uploaded full grids");
        MacGridState3D snapshot(9,1.f);resident.download(snapshot);
        double maximum=0;
        for(auto fields:{std::pair{reference.state().density(),std::span<const float>(snapshot.density())},
                         std::pair{reference.state().velocityX(),std::span<const float>(snapshot.velocityX())},
                         std::pair{reference.state().velocityY(),std::span<const float>(snapshot.velocityY())},
                         std::pair{reference.state().velocityZ(),std::span<const float>(snapshot.velocityZ())}})
            for(size_t i=0;i<fields.first.size();++i){require(std::isfinite(fields.second[i]),"Nonfinite resident field");maximum=std::max(maximum,double(std::abs(fields.first[i]-fields.second[i])));}
        std::cout<<"resident reflection="<<reflection<<" mcV="<<mcVelocity<<" mcS="<<mcSmoke<<" max error="<<maximum<<'\n';
        require(maximum<2e-4,"Resident/staged simulation mismatch");
        resident.resetState();for(int f=0;f<8;++f)resident.advance(1.f/60.f);
        MacGridState3D replay(9,1.f);resident.download(replay);
        require(std::equal(snapshot.density().begin(),snapshot.density().end(),replay.density().begin()),"Resident reset replay differs");
        resident.resetState();
        resident.parameters().emitterKickMultiplier=0;
        resident.parameters().buoyancy=0;
        resident.parameters().stirStrength=0;
        resident.advance(1.f/60); resident.download(replay);
        constant(replay.velocityX(),0.f); constant(replay.velocityY(),0.f); constant(replay.velocityZ(),0.f);
        require(*std::max_element(replay.density().begin(),replay.density().end())>0,"Disabling GPU kick disabled smoke");
        resident.setBoxSize(2.f);MacGridState3D resized(9,2.f);resident.download(resized);
        constant(resized.density(),0.f);
        resident.parameters().sourceStrength=0;resident.advance(1.f/60.f);resident.download(resized);
        constant(resized.velocityX(),0.f);constant(resized.velocityY(),0.f);constant(resized.velocityZ(),0.f);
    }
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

void advectionParity(const std::shared_ptr<GpuMacBackend>& backend) {
    MacGridState3D state(8, 4.f);
    std::fill(state.velocityX().begin(), state.velocityX().end(), .2f);
    std::fill(state.velocityY().begin(), state.velocityY().end(), -.1f);
    std::fill(state.velocityZ().begin(), state.velocityZ().end(), .3f);
    AdvectionContext context{state, velocityView(state), .1f, 3.f};
    const auto count = state.velocityX().size();
    std::vector<float> source(count), cpuOut(count), gpuOut(count);
    for (std::size_t i = 0; i < count; ++i) source[i] = float((i * 17) % 13) / 12.f;
    CpuAdvection sl;
    AdvectionWorkspace cpuScratch, gpuScratch;
    for (float dt : {.1f, 2.f}) for (float cfl : {3.f, .001f})
    for (bool correction : {false, true}) {
        context.dt = dt;
        context.cflLimit = cfl;
        const auto& cpu = sl;
        const auto selected = correction ? AdvectionMethod::MacCormack : AdvectionMethod::SemiLagrangian;
        GpuAdvection gpu(backend);
        for (auto axis : {FaceAxis::X, FaceAxis::Y, FaceAxis::Z}) {
            cpu.advectFace(context, axis, source, cpuOut, cpuScratch, selected);
            gpu.advectFace(context, axis, source, gpuOut, gpuScratch, selected);
            for (std::size_t i = 0; i < count; ++i) {
                require(std::isfinite(gpuOut[i]) && std::abs(cpuOut[i] - gpuOut[i]) < 2e-5f,
                        "CPU/GPU face advection differs after interface dispatch");
            }
        }
        rejects([&] { gpu.advectFace(context, FaceAxis::X, source, source, gpuScratch, selected); });
        rejects([&] { gpu.advectFace(context, FaceAxis::X, source, state.velocityX(), gpuScratch, selected); });
        rejects([&] { gpu.advectScalar(context, std::span<float>(gpuOut).first(1), gpuScratch, selected); });
    }
    for (float dt : {.1f, 2.f}) for (float cfl : {3.f, .001f})
    for (const auto method : {AdvectionMethod::SemiLagrangian, AdvectionMethod::MacCormack}) {
        context.dt = dt;
        context.cflLimit = cfl;
        std::vector<float> cpuDensity(source.begin(), source.begin() + state.cellCount());
        for (std::size_t i = 0; i < cpuDensity.size(); ++i)
            cpuDensity[i] = float((i * 17) % 13) / 12.f;
        auto gpuDensity = cpuDensity;
        GpuAdvection gpuSl(backend);
        sl.advectScalar(context, cpuDensity, cpuScratch, method);
        gpuSl.advectScalar(context, gpuDensity, gpuScratch, method);
        for (std::size_t i = 0; i < cpuDensity.size(); ++i)
            require(std::abs(cpuDensity[i] - gpuDensity[i]) < 2e-5f, "CPU/GPU scalar transport differs");
    }
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
    auto gpu = std::make_shared<GpuMacBackend>();
    VolumeRenderer renderer;
    try {
        // Loading GL entry points here must also initialize the functions used
        // by the separate GPU and renderer translation units.
        MacGridState3D state(8,4.f);
        const auto count=static_cast<int>(state.velocityX().size());
        gpu->init(state);
        require(gpu->initialized(),"GPU initialization failed");
        GpuPressureSolver gpuPressure(gpu);
        Projection3D& pressure = gpuPressure;
        MacGridState3D wrongSize(4);
        rejects([&] { pressure.project(wrongSize, {}); });
        state.velocityX()[state.idX(4,4,4)]=1.f;
        const float before=MacGridOperators3D::maxDivergence(state);
        ProjectionOptions projectionOptions;
        projectionOptions.relaxation = 1.5f;
        projectionOptions.linearSolve.maxIterations = 200;
        projectionOptions.linearSolve.fixedIterations = true;
        const auto projection = pressure.project(state, projectionOptions);
        require(!projection.residualAvailable, "GPU must not report an unmeasured residual");
        require(MacGridOperators3D::maxDivergence(state)<before*.01f,"GPU projection did not reduce divergence");

        state.reset();
        std::fill(state.density().begin(),state.density().end(),.5f);
        std::vector<float> source(count,.125f),output(count);
        AdvectionWorkspace workspace;
        for (bool correction : {false,true}) {
            GpuAdvection gpuMethod(gpu);
            const auto selected = correction ? AdvectionMethod::MacCormack : AdvectionMethod::SemiLagrangian;
            const Advection3D& method = gpuMethod;
            const AdvectionContext context{state, velocityView(state), .02f, 3.f};
            method.advectScalar(context, state.density(), workspace, selected);
            constant(state.density(),.5f);
            for (auto axis : {FaceAxis::X, FaceAxis::Y, FaceAxis::Z}) {
                method.advectFace(context, axis, source, output, workspace, selected);
                constant(output,.125f);
            }
        }
        require(glGetError()==GL_NO_ERROR,"OpenGL error in compute backend");
        advectionParity(gpu);
        residentParity();

        // Exercise the coordinator through all CPU/GPU pressure/advection
        // combinations, including switching backends on a live state.
        MacGridFluidSolver3D simulation(8);
        simulation.parameters().stirStrength = 0.f;
        simulation.parameters().emitterRadius = .18f;
        simulation.parameters().emitterCenterZ = .3f;
        simulation.parameters().projectIterations = 16;
        auto& fluid = simulation;
        for (bool reflect : {false, true}) for (bool correction : {false, true}) {
            simulation.parameters().reflection = reflect;
            simulation.parameters().macCormackVel = correction;
            simulation.parameters().macCormackSmoke = correction;
            for (bool gpuP : {false, true}) for (bool gpuA : {false, true}) {
                if (gpuP) simulation.setPressureSolver(std::make_unique<GpuPressureSolver>(gpu));
                else simulation.useDefaultPressureSolver();
                if (gpuA) simulation.setAdvection(std::make_unique<GpuAdvection>(gpu));
                else simulation.useDefaultAdvection();
                fluid.advance(1.f / 60.f);
                const auto stats = fluid.diagnostics();
                require(std::isfinite(stats.kineticEnergy) && std::isfinite(stats.maxDivergence),
                        "Backend switch produced invalid state");
                require(simulation.lastProjection().residualAvailable != gpuP,
                        "Pressure diagnostics belong to the wrong backend");
            }
        }
        require(glGetError()==GL_NO_ERROR,"OpenGL error after backend switching");

        std::weak_ptr<GpuMacBackend> lifetime;
        {
            MacGridFluidSolver3D owned(8);
            {
                auto backend = std::make_shared<GpuMacBackend>();
                backend->init(owned.state());
                lifetime = backend;
                owned.setAdvection(std::make_unique<GpuAdvection>(backend));
                owned.setPressureSolver(std::make_unique<GpuPressureSolver>(backend));
            }
            require(!lifetime.expired(), "Solver lost its GPU backend after construction");
            owned.advance(1.f / 60.f);
            owned.useDefaultAdvection();
            require(!lifetime.expired(), "Pressure stage must retain the shared backend");
        }
        require(lifetime.expired(), "GPU backend outlived its owning solver");
        require(glGetError()==GL_NO_ERROR, "GPU resource destruction failed");

        const RenderData data = volumeData(state);
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
                "VolumeRenderer produced no visible volume samples");
        // GPU-buffer presentation must produce the same pixels as CPU upload.
        GpuMacSimulation3D residentImage(state.resolution(),state.boxSize(),MacGridParameters3D{});
        residentImage.upload(state);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        renderer.render(residentImage.renderData(),64,64,mvp,cx,cy,cz,.5f,30.f,0.f,0.f,1.f,0.f,.1f,nullptr);
        std::vector<unsigned char> gpuPixels(pixels.size());
        glReadPixels(0,0,64,64,GL_RGB,GL_UNSIGNED_BYTE,gpuPixels.data());
        require(gpuPixels==pixels,"GPU buffer rendering differs from CPU texture upload");
        require(residentImage.downloadedBytes()==0,"Rendering downloaded GPU fields");
        std::cout << "PASS hidden OpenGL compute and volume-rendering smoke test\n";
    } catch (const std::exception& error) {
        std::cerr << "FAIL graphics smoke test: " << error.what() << '\n'; status=1;
    }
    renderer.shutdown();
    gpu->shutdown();
    require(!gpu->initialized(), "GPU shutdown must clear readiness");
    gpu->shutdown(); // Explicit shutdown followed by destruction is safe.
    glfwDestroyWindow(window);
    glfwTerminate();
    return status;
}
