#include "app/SimulationApp.hpp"
#include "composition/SimulationRegistry.hpp"
#include "renderer/Renderer2D.hpp"
#include "graphics/GlLoader.hpp"
#include <algorithm>
#include <cstdio>
#include <exception>
#include <variant>
#include <imgui.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_opengl3.h>

int SimulationApp::run(const SimulationConfig& config) {
    validateSimulationConfig(config);
    return config.dimension == Dimension::D2 ? run2D(config) : run3D(config);
}

int SimulationApp::run2D(const SimulationConfig& config) {
    SimulationRegistry registry;
    auto simulation = registry.create(config);
    IFluidSolver& fluid = simulation.solver();

    if (!glfwInit()) return 1;
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    GLFWwindow* window = glfwCreateWindow(960, 960, "vapor - 2D", nullptr, nullptr);
    if (!window) { glfwTerminate(); return 1; }
    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);
    if (!loadOpenGLFunctions()) { glfwDestroyWindow(window); glfwTerminate(); return 1; }
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGui::StyleColorsDark();
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init("#version 330");
    int status = 0;
    Renderer2D renderer;
    try {
        renderer.init(std::get<ImageRenderData>(simulation.presentation().presentationData()));
        auto controls = fluid.controls();
        bool paused = false;
        float exposure = 3.f, boxSize = fluid.boxSize();
        int frame = 0;
        while (!glfwWindowShouldClose(window) && (config.maxFrames == 0 || frame < config.maxFrames)) {
            glfwPollEvents();
            ImGui_ImplOpenGL3_NewFrame();
            ImGui_ImplGlfw_NewFrame();
            ImGui::NewFrame();
            if (!paused) fluid.advance(std::min(ImGui::GetIO().DeltaTime, 1.f / 30.f));
            int w = 0, h = 0;
            glfwGetFramebufferSize(window, &w, &h);
            glBindFramebuffer(GL_FRAMEBUFFER, 0);
            glClearColor(.025f, .025f, .035f, 1.f);
            glClear(GL_COLOR_BUFFER_BIT);
            renderer.render(std::get<ImageRenderData>(simulation.presentation().presentationData()),
                            w, h, exposure);

            ImGui::SetNextWindowPos({16, 16}, ImGuiCond_FirstUseEver);
            ImGui::Begin("2D Controls", nullptr, ImGuiWindowFlags_AlwaysAutoResize);
            ImGui::TextUnformatted("Configured 2D fluid solver");
            if (ImGui::Button(paused ? "Resume" : "Pause")) paused = !paused;
            ImGui::SameLine();
            if (ImGui::Button("Reset")) fluid.resetState();
            ImGui::SliderFloat("Source", &controls.sourceStrength, 0.f, 10.f);
            ImGui::SliderFloat("Buoyancy", &controls.buoyancy, 0.f, 10.f);
            ImGui::SliderFloat("Smoke decay", &controls.smokeDecay, 0.f, 2.f);
            ImGui::SliderFloat("Emitter radius", &controls.emitterRadius, .02f, .2f);
            ImGui::SliderFloat("Box size", &boxSize, .5f, 8.f);
            ImGui::SliderFloat("Exposure", &exposure, .1f, 10.f);
            ImGui::SliderInt("SOR iterations", &controls.pressureIterations, 10, 1000);
            ImGui::SliderFloat("SOR omega", &controls.pressureRelaxation, 1.f, 1.95f);
            const auto stats = fluid.diagnostics();
            const auto projection = fluid.pressureDiagnostics();
            ImGui::Text("Energy %.4f | max velocity %.3f", stats.kineticEnergy, stats.maxVelocity);
            ImGui::Text("Max divergence %.3e", stats.maxDivergence);
            ImGui::Text("Pressure residual %.3e (%d iterations)", projection.finalResidual, projection.iterations);
            if (projection.residualAvailable)
                ImGui::TextUnformatted(projection.converged ? "Pressure converged" : "Pressure iteration budget reached");
            ImGui::Text("FPS %.1f", ImGui::GetIO().Framerate);
            ImGui::End();
            fluid.setControls(controls);
            if (boxSize != fluid.boxSize()) fluid.setBoxSize(boxSize);
            ImGui::Render();
            ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
            glfwSwapBuffers(window);
            ++frame;
        }
        const auto stats = fluid.diagnostics();
        const auto solve = fluid.pressureDiagnostics();
        std::printf("2D: %d frames, energy=%.6f maxVelocity=%.6f maxDivergence=%.6g residual=%.6g iterations=%d converged=%s\n",
                    frame, stats.kineticEnergy, stats.maxVelocity, stats.maxDivergence,
                    solve.finalResidual, solve.iterations, solve.converged ? "yes" : "no");
    } catch (const std::exception& error) {
        std::fprintf(stderr, "2D simulation: %s\n", error.what());
        status = 1;
    }
    renderer.shutdown();
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
    glfwDestroyWindow(window);
    glfwTerminate();
    return status;
}
