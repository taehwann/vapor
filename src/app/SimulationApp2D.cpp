#include "app/SimulationApp.hpp"
#include "app/WindowContext.hpp"
#include "graphics/GlLoader.hpp"
#include "renderer/MacRenderData2D.hpp"
#include "renderer/Renderer2D.hpp"
#include "renderer/SimplicialWireData2D.hpp"
#include "renderer/TetrahedralWireRenderer.hpp"
#include "solvers/mac2d/MacGridFluidSolver2D.hpp"
#include "solvers/mac2d/ReflectionMacFluidSolver2D.hpp"
#include "solvers/simplicial2d/SimplicialFluidSolver2D.hpp"
#include "solvers/simplicial2d/SimplicialRasterizer2D.hpp"
#include <algorithm>
#include <cstdio>
#include <imgui.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_opengl3.h>

namespace {
// Common window mechanics. Solver construction and stepping stay in explicit loops.
struct Frame2D {
    GLFWwindow* window;
    const RunOptions& options;
    Renderer2D renderer;
    bool paused = false;
    float exposure = 3.f;
    int frame = 0;

    bool begin() {
        if (glfwWindowShouldClose(window) || (options.maxFrames > 0 && frame >= options.maxFrames))
            return false;
        glfwPollEvents();
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();
        return true;
    }
    float timestep() const {
        return options.fixedDt > 0.f ? options.fixedDt : std::min(ImGui::GetIO().DeltaTime, 1.f / 30.f);
    }
    void draw(const ImageRenderData& image) {
        int width, height;
        glfwGetFramebufferSize(window, &width, &height);
        glBindFramebuffer(GL_FRAMEBUFFER, 0);
        glClearColor(.025f, .025f, .035f, 1.f);
        glClear(GL_COLOR_BUFFER_BIT);
        renderer.render(image, width, height, exposure);
        ImGui::SetNextWindowPos({16, 16}, ImGuiCond_FirstUseEver);
        ImGui::SetNextWindowSizeConstraints({300, 100}, {600, std::max(100.f, ImGui::GetIO().DisplaySize.y - 32)});
        ImGui::Begin("2D Controls", nullptr, ImGuiWindowFlags_AlwaysAutoResize);
    }
    bool controls(float& source, float& buoyancy, float& decay) {
        if (ImGui::Button(paused ? "Resume" : "Pause"))
            paused = !paused;
        ImGui::SameLine();
        const bool reset = ImGui::Button("Reset");
        ImGui::SliderFloat("Source", &source, 0.f, 10.f);
        ImGui::SliderFloat("Buoyancy", &buoyancy, 0.f, 10.f);
        ImGui::SliderFloat("Smoke decay", &decay, 0.f, 2.f);
        ImGui::SliderFloat("Brightness", &exposure, .1f, 10.f);
        return reset;
    }
    void finish(const FluidDiagnostics& stats) {
        ImGui::Text("Energy %.4f | max velocity %.3f", stats.kineticEnergy, stats.maxVelocity);
        ImGui::Text("Max divergence %.3e", stats.maxDivergence);
        ImGui::Text("FPS %.1f", ImGui::GetIO().Framerate);
        ImGui::End();
        ImGui::Render();
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
        glfwSwapBuffers(window);
        ++frame;
    }
};
void pressureReport(const LinearSolveResult& solve) {
    ImGui::Text("Residual %.3e (%d iterations)", solve.finalResidual, solve.iterations);
    ImGui::TextUnformatted(solve.converged ? "Solve converged" : "Iteration budget reached");
}

int runMac2D(const RunOptions& options, GLFWwindow* window) {
    MacGridFluidSolver2D solver;
    Frame2D frame{window, options};
    frame.renderer.init(imageData(solver.state()));
    auto& p = solver.parameters();
    float boxSize = solver.boxSize();
    while (frame.begin()) {
        if (!frame.paused)
            solver.advance(frame.timestep());
        auto image = imageData(solver.state());
        image.historicalPalette = true;
        frame.draw(image);
        if (frame.controls(p.sourceStrength, p.buoyancy, p.smokeDecay))
            solver.resetState();
        if (ImGui::CollapsingHeader("Advanced")) {
            ImGui::SliderFloat("Emitter radius", &p.emitterRadius, .02f, .2f);
            ImGui::SliderFloat("Kick multiplier", &p.emitterKickMultiplier, 0.f, 3.f);
            ImGui::SliderInt("SOR iterations", &p.projectIterations, 10, 1000);
            ImGui::SliderFloat("SOR omega", &p.sorOmega, 1.f, 1.95f);
            pressureReport(solver.lastProjection().linearSolve);
            ImGui::SliderFloat("Box size", &boxSize, .5f, 8.f);
        }
        if (boxSize != solver.boxSize())
            solver.setBoxSize(boxSize);
        frame.finish(solver.diagnostics());
    }
    const auto stats = solver.diagnostics();
    std::printf("2D Mac: %d frames, energy=%.6f maxVelocity=%.6f maxDivergence=%.6g\n", frame.frame,
                stats.kineticEnergy, stats.maxVelocity, stats.maxDivergence);
    return 0;
}

int runReflection2D(const RunOptions& options, GLFWwindow* window, AdvectionMethod method) {
    ReflectionMacFluidSolver2D solver(method);
    Frame2D frame{window, options};
    frame.renderer.init(ImageRenderData{solver.smoke, solver.n, solver.n, 1.f, 1.f, true});
    auto& p = solver.parameters();
    while (frame.begin()) {
        if (!frame.paused)
            solver.advance(frame.timestep());
        const ImageRenderData image{solver.smoke, solver.n, solver.n, 1.f, 1.f, true};
        frame.draw(image);
        if (frame.controls(p.sourceStrength, p.buoyancy, p.smokeDecay))
            solver.resetState();
        if (ImGui::CollapsingHeader("Advanced")) {
            ImGui::SliderFloat("Emitter radius", &p.emitterRadius, .02f, .2f);
            ImGui::TextUnformatted("Historical reflection; closed free-slip walls");
            ImGui::SliderFloat("Kick multiplier", &p.emitterKickMultiplier, 0.f, 3.f);
            ImGui::SliderInt("SOR iterations", &p.projectIterations, 10, 200);
            ImGui::SliderFloat("SOR omega", &p.sorOmega, 1.f, 1.95f);
            ImGui::SliderFloat("Stir", &p.stirStrength, 0.f, 2.f);
            ImGui::TextUnformatted("Fixed-budget pressure; inspect divergence");
        }
        frame.finish(solver.diagnostics());
    }
    const auto stats = solver.diagnostics();
    std::printf("2D Reflection: %d frames, energy=%.6f maxVelocity=%.6f maxDivergence=%.6g\n", frame.frame,
                stats.kineticEnergy, stats.maxVelocity, stats.maxDivergence);
    return 0;
}

int runSimplicial2D(const RunOptions& options, GLFWwindow* window) {
    SimplicialFluidSolver2D solver{SimplicialFluidState2D{SimplicialMesh2D::teapot()}};
    SimplicialRasterizer2D rasterizer(solver.state());
    Frame2D frame{window, options};
    frame.renderer.init(rasterizer.renderData());
    auto& p = solver.parameters();
    float boxSize = solver.boxSize();
    TetrahedralWireRenderer wires, outline;
    auto updateWires = [&] {
        wires.init(simplicialWires2D(solver.domain()));
        outline.init(simplicialWires2D(solver.domain(),true));
    };
    updateWires();
    int visualization = 2;
    float simulationSpeed = 1;
    double simulatedTime = 0;
    bool primal = true, dual = false, showOutline = true, capture = false;
    std::string captureStatus;
    while (frame.begin()) {
        if (!frame.paused) {
            float dt = frame.timestep() * simulationSpeed;
            solver.advance(dt);
            simulatedTime += dt;
        }
        const auto image = rasterizer.renderData();
        frame.draw(image);
        int w,h;glfwGetFramebufferSize(window,&w,&h);
        if(visualization==1) {glClearColor(.025f,.025f,.035f,1);glClear(GL_COLOR_BUFFER_BIT);}
        float matrix[16];simplicialMvp2D(matrix,w,h,boxSize);
        if(visualization!=0)wires.render(matrix,w,h,primal,dual,false,boxSize+1,false);
        if(showOutline)outline.render(matrix,w,h,true,false,false,boxSize+1,false);
        if(capture) {
            try {saveWireframeImage("simplicial2d-render.ppm",w,h);captureStatus="Saved simplicial2d-render.ppm";}
            catch(const std::exception& e){captureStatus=e.what();}
            capture=false;
        }
        ImGui::TextUnformatted("Teapot | inviscid, slip walls");
        ImGui::Text("%zu vertices | %zu triangles",solver.domain().vertexCount(),solver.domain().triangleCount());
        ImGui::Combo("Visualization",&visualization,"Smoke\0Wireframe\0Smoke + wireframe\0");
        if(visualization!=0) {ImGui::Checkbox("Primal (white)",&primal);ImGui::Checkbox("Dual (green)",&dual);}
        ImGui::Checkbox("Domain outline", &showOutline);
        if(ImGui::Button("Save render (.ppm)"))capture=true;
        if(!captureStatus.empty())ImGui::TextUnformatted(captureStatus.c_str());
        if (frame.controls(p.sourceStrength, p.buoyancy, p.smokeDecay)) {
            solver.resetState(); simulatedTime = 0;
        }
        ImGui::SliderFloat("Simulation speed", &simulationSpeed, .1f, 3.f);
        ImGui::SliderFloat("Emitter stir", &p.stirStrength, 0.f, 1.f);
        ImGui::Text("Simulated time: %.2f s", simulatedTime);
        if (ImGui::Button("Persistent smoke preset")) { p.sourceStrength = 1; p.smokeDecay = 0; }
        if (ImGui::CollapsingHeader("Advanced")) {
            ImGui::SliderFloat("Emitter radius", &p.emitterRadius, .02f, .3f);
            ImGui::SliderFloat("Emitter horizontal", &p.emitterCenterS, 0.f, 1.f);
            ImGui::SliderFloat("Emitter vertical", &p.emitterCenterT, 0.f, 1.f);
            ImGui::SliderInt("CG iterations", &p.cgIterations, 10, 1000);
            ImGui::SliderFloat("Advection CFL", &p.advectionCfl, .1f, 1.f);
            ImGui::InputDouble("CG absolute tolerance", &p.cgAbsoluteTolerance, 0, 0, "%.2e");
            ImGui::InputDouble("CG relative tolerance", &p.cgRelativeTolerance, 0, 0, "%.2e");
            p.cgAbsoluteTolerance = std::clamp(p.cgAbsoluteTolerance, 1e-12, 1e-3);
            p.cgRelativeTolerance = std::clamp(p.cgRelativeTolerance, 1e-12, 1e-3);
            pressureReport(solver.lastRecovery());
            ImGui::Text("Substeps %d", solver.lastSubstepCount());
            ImGui::SliderFloat("Box size", &boxSize, .5f, 8.f);
        }
        if (boxSize != solver.boxSize()) {
            solver.setBoxSize(boxSize);
            simulatedTime = 0;
            updateWires();
        }
        frame.finish(solver.diagnostics());
    }
    const auto stats = solver.diagnostics();
    std::printf("2D Simplicial: %d frames, energy=%.6f maxVelocity=%.6f maxDivergence=%.6g\n", frame.frame,
                stats.kineticEnergy, stats.maxVelocity, stats.maxDivergence);
    return 0;
}

} // namespace
int SimulationApp::run2D(SolverChoice choice, const RunOptions& options) {
    WindowContext context(960, 960, "vapor - 2D", options.hidden);
    switch (choice) {
    case SolverChoice::Mac2DSimple:
        return runMac2D(options, context.window());
    case SolverChoice::Mac2DReflectionSL:
        return runReflection2D(options, context.window(), AdvectionMethod::SemiLagrangian);
    case SolverChoice::Mac2DReflectionMC:
        return runReflection2D(options, context.window(), AdvectionMethod::MacCormack);
    case SolverChoice::Simplicial2DTeapot:
        return runSimplicial2D(options, context.window());
    default:
        throw std::invalid_argument("Expected a 2D solver");
    }
}
