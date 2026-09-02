#include "app/SimulationApp.hpp"
#include "composition/GpuSimulationRegistration.hpp"
#include "composition/SimulationRegistry.hpp"
#include "graphics/GlLoader.hpp"
#include "renderer/OrbitCamera.hpp"
#include "renderer/Renderer.hpp"

#include <algorithm>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <memory>
#include <variant>
#include <vector>

#include <imgui.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_opengl3.h>

int SimulationApp::run3D(const SimulationConfig& config) {
    if (!glfwInit()) return 1;
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    GLFWwindow* window = glfwCreateWindow(960, 960, "vapor", nullptr, nullptr);
    if (!window) { glfwTerminate(); return 1; }
    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);
    if (!loadOpenGLFunctions()) { glfwDestroyWindow(window); glfwTerminate(); return 1; }

    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGui::StyleColorsDark();
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init("#version 330");
    auto registry = std::make_unique<SimulationRegistry>();
    registerGpuSimulationModules(*registry);
    auto simulation = std::make_unique<SimulationInstance>(registry->create(config));
    IFluidSolver& fluidSolver = simulation->solver();
    const SimulationConfig& effectiveConfig = simulation->configuration();
    auto controls = fluidSolver.controls();
    const bool useGpuPressure = effectiveConfig.projectionSolver == ProjectionSolverType::GpuSor;
    const bool useGpuAdvection = effectiveConfig.advectionBackend == ComputeBackend::Gpu;
    Renderer renderer;

    FILE* logFile = std::fopen("vapor_console.log", "w");
    auto logPrint = [&](const char* fmt, ...) {
        va_list args1, args2;
        va_start(args1, fmt); va_start(args2, fmt);
        std::vprintf(fmt, args1);
        if (logFile) { std::vfprintf(logFile, fmt, args2); std::fflush(logFile); }
        va_end(args2); va_end(args1);
        };

    renderer.init(std::get<RenderData>(simulation->presentation().presentationData()));

    bool paused = false, debugPrintOn = true;
    const bool useReflection = effectiveConfig.fluidSolver == FluidSolverType::Reflection;
    int maxFrames = config.maxFrames;
    float alphaMul = 30.0f;
    int debugFrame = 0, simulationFrame = 0;
    const float maxDt = 0.033f;
    float azimuth = -1.2f, elevation = 0.3f, dist = fluidSolver.boxSize() * 2.5f, stepScale = 0.5f;
    float boxSize = fluidSolver.boxSize();
    float lightX = 0.3f, lightY = 0.4f, lightZ = 1.0f;
    float shadowStr = 0.3f, shadowStep = 0.1f;
    bool dragging = false;
    double lastX = 0, lastY = 0;

    while (!glfwWindowShouldClose(window)) {
        glfwPollEvents();
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        if (!paused) {
            float dt = std::min(ImGui::GetIO().DeltaTime, maxDt);
            fluidSolver.advance(dt);
            ++simulationFrame;
            if (debugPrintOn) {
                ++debugFrame;
                const auto density =
                    std::get<RenderData>(simulation->presentation().presentationData()).density;
                float dmax = *std::max_element(density.begin(), density.end());
                const auto diagnostics = fluidSolver.diagnostics();
                float ke = diagnostics.kineticEnergy;
                float maxDiv = diagnostics.maxDivergence;
                float maxV = diagnostics.maxVelocity;
                logPrint("f %4d | ke=%.6f maxDiv=%.6f maxV=%.3f smoke=%.3f | %s\n",
                    debugFrame, ke, maxDiv, maxV, dmax,
                    useReflection ? "REFLECT" : "NOREFLECT");
            }
            if (maxFrames > 0 && simulationFrame >= maxFrames) {
                logPrint("--- AUTO-QUIT after %d frames ---\n", maxFrames);
                glfwSetWindowShouldClose(window, GLFW_TRUE);
            }
        }

        int w, h; glfwGetFramebufferSize(window, &w, &h);
        glBindFramebuffer(GL_FRAMEBUFFER, 0);
        glClearColor(0.03f, 0.04f, 0.07f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        const RenderData renderData =
            std::get<RenderData>(simulation->presentation().presentationData());
        float mvp[16], cx, cy, cz;
        buildMVP(mvp, cx, cy, cz, azimuth, elevation, dist, float(w) / float(h), renderData.boxSize);
        float ld = 1.0f / std::sqrt(lightX * lightX + lightY * lightY + lightZ * lightZ);
        renderer.render(renderData, w, h, mvp, cx, cy, cz, stepScale, alphaMul,
            lightX * ld, lightY * ld, lightZ * ld, shadowStr, shadowStep, logFile);

        if (debugPrintOn && debugFrame % 5 == 0) {
            int ascW = 80, ascH = 24;
            std::vector<unsigned char> full(w * h * 3);
            glReadBuffer(GL_BACK);
            glPixelStorei(GL_PACK_ALIGNMENT, 1);
            glReadPixels(0, 0, w, h, GL_RGB, GL_UNSIGNED_BYTE, full.data());
            glPixelStorei(GL_PACK_ALIGNMENT, 4);
            logPrint("  framebuffer %dx%d:\n", w, h);
            for (int ry = ascH - 1; ry >= 0; --ry) {
                logPrint("  ");
                for (int rx = 0; rx < ascW; ++rx) {
                    int sx = rx * w / ascW;
                    int sy = ry * h / ascH;
                    int i = (sy * w + sx) * 3;
                    float lum = 0.2126f * full[i] + 0.7152f * full[i + 1] + 0.0722f * full[i + 2];
                    const char* ramp = " .:-=+*#%@";
                    int idx = int(lum / 255.0f * 9.0f);
                    idx = idx < 0 ? 0 : idx > 9 ? 9 : idx;
                    logPrint("%c", ramp[idx]);
                }
                logPrint("\n");
            }
        }

        ImGui::SetNextWindowPos({ 16, 16 }, ImGuiCond_FirstUseEver);
        ImGui::Begin("Controls", nullptr, ImGuiWindowFlags_AlwaysAutoResize);
        if (ImGui::Button(paused ? "Resume" : "Pause")) paused = !paused;
        ImGui::SameLine();
        if (ImGui::Button("Reset")) {
            fluidSolver.resetState();
        }
        ImGui::Text("Fluid algorithm: %s", useReflection ? "reflection" : "simple");
        ImGui::Text("Advection: %s",
                    effectiveConfig.advection == AdvectionType::MacCormack
                        ? "MacCormack" : "Semi-Lagrangian");
        ImGui::Text("Advection backend: %s", useGpuAdvection ? "GPU" : "CPU");
        ImGui::Text("Projection: %s", useGpuPressure ? "GPU SOR" : "CPU SOR");
        ImGui::Checkbox("Debug print", &debugPrintOn);
        ImGui::SliderFloat("Buoyancy", &controls.buoyancy, 0.f, 10.f);
        ImGui::SliderFloat("Source", &controls.sourceStrength, 0.f, 3.f);
        ImGui::SliderFloat("Smoke decay", &controls.smokeDecay, 0.f, 2.f);
        ImGui::SliderFloat("Box size", &boxSize, 0.5f, 8.f);
        ImGui::SliderFloat("Emitt radius", &controls.emitterRadius, 0.01f, 0.2f);
        ImGui::SliderFloat("Stir", &controls.stirStrength, 0.f, 2.f);
        ImGui::SliderFloat("Alpha mul", &alphaMul, 5.f, 100.f);
        ImGui::SliderInt("SOR its", &controls.pressureIterations, 10, 500);
        const auto pressure = fluidSolver.pressureDiagnostics();
        if (pressure.residualAvailable) {
            const auto& solve = pressure;
            ImGui::Text("CPU pressure residual: %.3e (%d its)", solve.finalResidual, solve.iterations);
        }
        ImGui::SliderFloat("MC CFL", &controls.advectionCfl, 0.5f, 10.f);
        ImGui::SliderFloat("Step", &stepScale, 0.1f, 3.f);
        ImGui::SliderFloat("Shadow str", &shadowStr, 0.f, 5.f);
        ImGui::SliderFloat("Shadow step", &shadowStep, 0.05f, 1.f);
        ImGui::SliderFloat3("Light dir", &lightX, -1.f, 1.f);
        ImGui::Text("Drag mouse to orbit, scroll to zoom");
        ImGui::Text("FPS: %.1f", ImGui::GetIO().Framerate);
        ImGui::End();

        fluidSolver.setControls(controls);
        if (boxSize != fluidSolver.boxSize()) {
            fluidSolver.setBoxSize(boxSize);
            renderer.updateDomain(
                std::get<RenderData>(simulation->presentation().presentationData()));
            dist = boxSize * 2.5f;
        }

        ImGui::Render();
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

        if (!ImGui::GetIO().WantCaptureMouse && ImGui::GetIO().MouseWheel != 0.f) {
            dist = std::clamp(dist - ImGui::GetIO().MouseWheel * dist * 0.1f,
                              fluidSolver.boxSize() * 0.3f, fluidSolver.boxSize() * 8.f);
        }

        if (ImGui::IsMouseDown(0) && !ImGui::GetIO().WantCaptureMouse) {
            double mx, my; glfwGetCursorPos(window, &mx, &my);
            if (!dragging) { dragging = true; lastX = mx; lastY = my; }
            else { azimuth += float((mx - lastX) * 0.005f); elevation = std::clamp(elevation + float((lastY - my) * 0.005f), -1.5f, 1.5f); lastX = mx; lastY = my; }
        }
        else dragging = false;

        glfwSwapBuffers(window);
    }
    renderer.shutdown();
    simulation.reset();
    registry.reset(); // Releases GPU resources while the OpenGL context is current.
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
    glfwDestroyWindow(window);
    glfwTerminate();
    if (logFile) std::fclose(logFile);
    return 0;
}
