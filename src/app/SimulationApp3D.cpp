#include "app/SimulationApp.hpp"
#include "app/WindowContext.hpp"
#include "graphics/GlLoader.hpp"
#include "renderer/MacRenderData.hpp"
#include "renderer/OrbitCamera.hpp"
#include "renderer/VolumeRenderer.hpp"
#include "solvers/mac3d/MacGridFluidSolver3D.hpp"
#include "solvers/mac3d/MacGridOperators3D.hpp"
#include "solvers/mac3d/gpu/GpuMacSimulation3D.hpp"

#include <algorithm>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <memory>

#include <vector>

#include <imgui.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_opengl3.h>

int SimulationApp::run3D(MacAlgorithm algorithm, bool useGpu, const RunOptions& options) {
    WindowContext context(960, 960, "vapor", options.hidden);
    GLFWwindow* window = context.window();
    std::unique_ptr<MacGridFluidSolver3D> simulation;
    std::unique_ptr<GpuMacSimulation3D> resident;
    if (useGpu) {
        try {
            resident = std::make_unique<GpuMacSimulation3D>(algorithm);
        } catch (const std::exception& error) {
            std::fprintf(stderr, "GPU initialization failed; using CPU: %s\n", error.what());
        }
    }
    if (!resident)
        simulation = std::make_unique<MacGridFluidSolver3D>(algorithm);
    auto& controls = resident ? resident->parameters() : simulation->parameters();
    auto renderDataForFrame = [&] {
        return resident ? resident->renderData() : volumeData(simulation->state());
    };
    auto reset = [&] {
        if (resident)
            resident->resetState();
        else
            simulation->resetState();
    };
    const bool useGpuPressure = bool(resident);
    const bool useGpuAdvection = bool(resident);
    VolumeRenderer renderer;

    std::unique_ptr<FILE, decltype(&std::fclose)> log(std::fopen("vapor_console.log", "w"), &std::fclose);
    FILE* logFile = log.get();
    auto logPrint = [&](const char* fmt, ...) {
        va_list args1, args2;
        va_start(args1, fmt);
        va_start(args2, fmt);
        std::vprintf(fmt, args1);
        if (logFile) {
            std::vfprintf(logFile, fmt, args2);
            std::fflush(logFile);
        }
        va_end(args2);
        va_end(args1);
    };

    logPrint("3D MAC: resolution=%d box=%.6g seed=%u fixed_dt=%.9g pressure=%s advection=%s\n",
             renderDataForFrame().width, renderDataForFrame().boxSize, controls.randomSeed, options.fixedDt,
             useGpuPressure ? "GPU" : "CPU", useGpuAdvection ? "GPU" : "CPU");
    renderer.init(renderDataForFrame());

    bool paused = false, debugPrintOn = false;
    const bool useReflection = controls.reflection;
    int maxFrames = options.maxFrames;
    float alphaMul = 30.0f;
    int simulationFrame = 0;
    const float maxDt = 0.033f;
    float azimuth = -1.2f, elevation = 0.3f, dist = renderDataForFrame().boxSize * 2.5f, stepScale = 0.5f;
    float boxSize = renderDataForFrame().boxSize;
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
            float dt = options.fixedDt > 0.f ? options.fixedDt : std::min(ImGui::GetIO().DeltaTime, maxDt);
            if (resident)
                resident->advance(dt);
            else
                simulation->advance(dt);
            ++simulationFrame;
            if (debugPrintOn && simulationFrame % 60 == 0) {
                // Explicit opt-in snapshot. Normal frames never download GPU fields.
                std::unique_ptr<MacGridState3D> snapshot;
                if (resident) {
                    snapshot =
                        std::make_unique<MacGridState3D>(renderDataForFrame().width, resident->boxSize());
                    resident->download(*snapshot);
                }
                const auto& state = snapshot ? *snapshot : simulation->state();
                float dmax = *std::max_element(state.density().begin(), state.density().end());
                const float maxDiv = MacGridOperators3D::maxDivergence(state);
                logPrint("f %4d | maxDiv=%.6f smoke=%.3f | %s\n", simulationFrame, maxDiv, dmax,
                         useReflection ? "REFLECT" : "NOREFLECT");
            }
            if (maxFrames > 0 && simulationFrame >= maxFrames) {
                logPrint("--- AUTO-QUIT after %d frames ---\n", maxFrames);
                glfwSetWindowShouldClose(window, GLFW_TRUE);
            }
        }

        int w, h;
        glfwGetFramebufferSize(window, &w, &h);
        glBindFramebuffer(GL_FRAMEBUFFER, 0);
        glClearColor(0.03f, 0.04f, 0.07f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        const RenderData renderData = renderDataForFrame();
        float mvp[16], cx, cy, cz;
        buildMVP(mvp, cx, cy, cz, azimuth, elevation, dist, float(w) / float(h), renderData.boxSize);
        float ld = 1.0f / std::sqrt(lightX * lightX + lightY * lightY + lightZ * lightZ);
        renderer.render(renderData, w, h, mvp, cx, cy, cz, stepScale, alphaMul, lightX * ld, lightY * ld,
                        lightZ * ld, shadowStr, shadowStep, logFile);

        ImGui::SetNextWindowPos({16, 16}, ImGuiCond_FirstUseEver);
        ImGui::Begin("Controls", nullptr, ImGuiWindowFlags_AlwaysAutoResize);
        if (ImGui::Button(paused ? "Resume" : "Pause"))
            paused = !paused;
        ImGui::SameLine();
        if (ImGui::Button("Reset")) {
            reset();
        }
        ImGui::SliderFloat("Buoyancy", &controls.buoyancy, 0.f, 10.f);

        ImGui::SliderFloat("Source", &controls.sourceStrength, 0.f, 3.f);
        ImGui::SliderFloat("Smoke decay", &controls.smokeDecay, 0.f, 2.f);
        ImGui::SliderFloat("Opacity", &alphaMul, 5.f, 100.f);
        if (ImGui::CollapsingHeader("Advanced")) {
            ImGui::Text("Fluid algorithm: %s", useReflection ? "reflection" : "simple");
            ImGui::Text("Advection: %s", controls.macCormackSmoke ? "MacCormack" : "Semi-Lagrangian");
            ImGui::Text("Advection backend: %s", useGpuAdvection ? "GPU" : "CPU");
            ImGui::Text("Projection: %s", useGpuPressure ? "GPU SOR" : "CPU SOR");
            ImGui::TextUnformatted(resident ? "Fields stay on GPU" : "CPU / staged fields");
            ImGui::Checkbox("Debug snapshot (every 60 steps)", &debugPrintOn);
            ImGui::SliderFloat("Kick multiplier", &controls.emitterKickMultiplier, 0.f, 3.f);
            ImGui::SliderFloat("Box size", &boxSize, 0.5f, 8.f);
            ImGui::SliderFloat("Emitt radius", &controls.emitterRadius, 0.01f, 0.2f);
            ImGui::SliderFloat("Stir", &controls.stirStrength, 0.f, 2.f);

            ImGui::SliderInt("SOR its", &controls.projectIterations, 10, 500);
            const auto pressure = resident ? resident->lastProjection() : simulation->lastProjection();
            if (pressure.residualAvailable) {
                const auto& solve = pressure.linearSolve;
                ImGui::Text("CPU pressure residual: %.3e (%d its)", solve.finalResidual, solve.iterations);
            }
            ImGui::SliderFloat("MC CFL", &controls.cflMc, 0.5f, 10.f);
            ImGui::SliderFloat("Step", &stepScale, 0.1f, 3.f);
            ImGui::SliderFloat("Shadow str", &shadowStr, 0.f, 5.f);
            ImGui::SliderFloat("Shadow step", &shadowStep, 0.05f, 1.f);
            ImGui::SliderFloat3("Light dir", &lightX, -1.f, 1.f);
        }
        ImGui::Text("Drag mouse to orbit, scroll to zoom");
        ImGui::Text("FPS: %.1f", ImGui::GetIO().Framerate);
        ImGui::End();

        if (boxSize != renderDataForFrame().boxSize) {
            if (resident)
                resident->setBoxSize(boxSize);
            else
                simulation->setBoxSize(boxSize);
            renderer.updateDomain(renderDataForFrame());
            dist = boxSize * 2.5f;
        }

        ImGui::Render();
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

        if (!ImGui::GetIO().WantCaptureMouse && ImGui::GetIO().MouseWheel != 0.f) {
            dist = std::clamp(dist - ImGui::GetIO().MouseWheel * dist * 0.1f,
                              renderDataForFrame().boxSize * 0.3f, renderDataForFrame().boxSize * 8.f);
        }

        if (ImGui::IsMouseDown(0) && !ImGui::GetIO().WantCaptureMouse) {
            double mx, my;
            glfwGetCursorPos(window, &mx, &my);
            if (!dragging) {
                dragging = true;
                lastX = mx;
                lastY = my;
            } else {
                azimuth += float((mx - lastX) * 0.005f);
                elevation = std::clamp(elevation + float((lastY - my) * 0.005f), -1.5f, 1.5f);
                lastX = mx;
                lastY = my;
            }
        } else
            dragging = false;

        glfwSwapBuffers(window);
    }
    return 0;
}
