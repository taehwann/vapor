#include "SimulationApp.hpp"
#include "WindowContext.hpp"
#include "graphics/GlLoader.hpp"
#include "renderer/OrbitCamera.hpp"
#include "renderer/SimplicialRenderData3D.hpp"
#include "renderer/TetrahedralWireRenderer.hpp"
#include "renderer/SimplicialWireData3D.hpp"
#include "renderer/SilhouetteRenderer.hpp"
#include "renderer/VolumeRenderer.hpp"
#include "solvers/simplicial3d/SimplicialFluidSolver3D.hpp"
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <imgui.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_opengl3.h>
#include <limits>

int SimulationApp::runSimplicial3D(const RunOptions &options) {
  WindowContext context(960, 960, "vapor - 3D simplicial fluids",
                        options.hidden);
  auto *window = context.window();
  SimplicialFluidSolver3D solver(simplicial3d::Mesh::loadDomain(
      options.assetDirectory / "examples/simplicial3d-bunny/bunny-fluid.tet", 2));
  SimplicialRenderData3D raster(64);
  raster.update(solver);
  VolumeRenderer renderer;
  renderer.init(raster.renderData());
  TetrahedralWireRenderer wireRenderer;
  SilhouetteRenderer boundaryRenderer;
  auto updateWires = [&] {
    wireRenderer.init(simplicialWires3D(solver.domain(), solver.dual()));
    boundaryRenderer.init(solver.domain());
  };
  updateWires();
  int visualization = 0;
  bool showPrimal = true, showDual = true, wireDepth = true, saveRender = false;
  float cutaway = 1;
  bool showBoundary = true;
  float simulationSpeed = 1, renderStep = .5f;
  double simulationTime = 0;
  double volumeUpdateMs = 0;
  std::string captureStatus;
  auto &p = solver.parameters();
  bool paused = false, dragging = false;
  int frames = 0;
  float azimuth = -1.2f, elevation = .3f, dist = solver.boxSize() * 2.5f,
        opacity = 30, box = solver.boxSize();
  double lastX = 0, lastY = 0;
  std::printf("3D simplicial: %zu tetrahedra, %zu face fluxes, %zu edge "
              "potentials; circumcentric dual\n",
              solver.domain().tets.size(), solver.domain().faces.size(),
              solver.domain().edges.size());
  while (!glfwWindowShouldClose(window)) {
    glfwPollEvents();
    ImGui_ImplOpenGL3_NewFrame();
    ImGui_ImplGlfw_NewFrame();
    ImGui::NewFrame();
    if (!paused) {
      float dt = simulationSpeed * (options.fixedDt > 0 ? options.fixedDt
                                     : std::min(ImGui::GetIO().DeltaTime, .033f));
      solver.advance(dt);
      simulationTime += dt;
      auto updateStart = std::chrono::steady_clock::now();
      raster.update(solver);
      volumeUpdateMs = std::chrono::duration<double, std::milli>(
                           std::chrono::steady_clock::now() - updateStart)
                           .count();
      if (options.maxFrames > 0 && ++frames >= options.maxFrames)
        glfwSetWindowShouldClose(window, GLFW_TRUE);
    }
    int w, h;
    glfwGetFramebufferSize(window, &w, &h);
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    glClearColor(.03f, .04f, .07f, 1);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    float mvp[16], cx, cy, cz;
    buildMVP(mvp, cx, cy, cz, azimuth, elevation, dist,
             float(w) / std::max(1, h), box);
    if (visualization != 1)
      renderer.render(raster.renderData(), w, h, mvp, cx, cy, cz, renderStep, opacity,
                      .2683f, .3578f, .8944f, .3f, .1f, nullptr);
    if (visualization != 0) {
      glDepthMask(GL_TRUE);
      glClear(GL_DEPTH_BUFFER_BIT);
      wireRenderer.render(mvp, w, h, showPrimal, showDual, false,
                          cutaway >= 1 ? std::numeric_limits<float>::max()
                                       : cutaway * box,
                          wireDepth);
    }
    if (showBoundary) boundaryRenderer.render(mvp, w, h);
    if (saveRender) {
      try {
        saveWireframeImage("simplicial3d-render.ppm", w, h);
        captureStatus = "Saved simplicial3d-render.ppm";
      } catch (const std::exception &e) {
        captureStatus = e.what();
      }
      saveRender = false;
    }
    ImGui::SetNextWindowPos({16, 16}, ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowSizeConstraints({330, 100}, {650, std::max(100.f, ImGui::GetIO().DisplaySize.y - 32)});
    ImGui::Begin("3D simplicial fluids", nullptr,
                 ImGuiWindowFlags_AlwaysAutoResize);
    if (ImGui::Button(paused ? "Resume" : "Pause"))
      paused = !paused;
    ImGui::SameLine();
    if (ImGui::Button("Reset")) {
      solver.resetState();
      simulationTime = 0;
      raster.update(solver);
    }
    ImGui::Text("%s | %zu tetrahedra", "Bunny",
                solver.domain().tets.size());

    ImGui::Combo("Visualization", &visualization,
                 "Smoke + outline\0Mesh inspection\0Smoke + mesh inspection\0");
    ImGui::Checkbox("Bunny silhouette outline", &showBoundary);
    if (visualization != 0) {
      ImGui::Checkbox("Primal (white)", &showPrimal);
      ImGui::Checkbox("Dual (green)", &showDual);
      ImGui::TextUnformatted("Circumcentric dual (solver geometry)");
      ImGui::SliderFloat("Wire X cutaway", &cutaway, 0, 1);
      ImGui::Checkbox("Depth test wires", &wireDepth);
    }
    if (ImGui::Button("Save render (.ppm)"))
      saveRender = true;
    if (!captureStatus.empty())
      ImGui::TextUnformatted(captureStatus.c_str());
    ImGui::SliderFloat("Buoyancy", &p.buoyancy, 0, 10);
    ImGui::SliderFloat("Emitter stir", &p.stirStrength, 0, 10);
    if (p.stirStrength > 0) {
      ImGui::SliderFloat("Stir radius", &p.stirRadius, .04f, .4f);
      ImGui::SliderFloat("Stir frequency (Hz)", &p.stirFrequency, 0, 4);
      ImGui::TextUnformatted("Localized stirring follows the enabled smoke emitter.");
    }
    ImGui::Checkbox("Continuous emitter", &p.emitterEnabled);
    ImGui::SliderFloat("Source concentration", &p.sourceStrength, 0, 5);
    ImGui::SliderFloat("Emitter radius", &p.emitterRadius, .02f, .3f);
    ImGui::SliderFloat("Emitter X", &p.emitterX, 0, 1);
    ImGui::SliderFloat("Emitter Y", &p.emitterY, 0, 1);
    ImGui::SliderFloat("Emitter Z (height)", &p.emitterZ, 0, 1);
    if (ImGui::Button("Persistent smoke preset")) {
      p.emitterEnabled = true; p.sourceStrength = 1; p.smokeDecay = 0;
      p.emitterX = .5f; p.emitterY = .5f; p.emitterZ = .2f; p.emitterRadius = .12f;
    }
    ImGui::SliderFloat("Simulation speed", &simulationSpeed, .1f, 3.f);

    ImGui::SliderFloat("Smoke decay", &p.smokeDecay, 0, 2);
    ImGui::SliderFloat("Opacity", &opacity, 5, 100);
    if (ImGui::CollapsingHeader("Advanced")) {
      const auto &timing = solver.timings();
      double stepMs = timing.emissionMs + timing.advectionMs + timing.forcesMs +
                      timing.recoveryMs + timing.smokeMs;
      ImGui::Text("CPU step: %.1f ms | volume update: %.2f ms", stepMs,
                  volumeUpdateMs);
      ImGui::Text("Substeps: %d", timing.substeps);
      ImGui::SliderFloat("Advection CFL", &p.advectionCfl, .1f, .8f);
      ImGui::SliderFloat("Render step", &renderStep, .25f, 2.f);
      ImGui::InputDouble("CG absolute tolerance", &p.cgAbsoluteTolerance, 0, 0, "%.2e");
      ImGui::InputDouble("CG relative tolerance", &p.cgRelativeTolerance, 0, 0, "%.2e");
      p.cgAbsoluteTolerance = std::clamp(p.cgAbsoluteTolerance, 1e-12, 1e-3);
      p.cgRelativeTolerance = std::clamp(p.cgRelativeTolerance, 1e-12, 1e-3);
      ImGui::SliderFloat("Domain size", &box, .5f, 8);
      ImGui::SliderInt("CG iterations", &p.cgIterations, 100, 3000);
    }
    const auto d = solver.diagnostics();
    ImGui::Text("Divergence: %.2e | Wall flux: %.2e", d.maxDivergence,
                d.maxBoundaryFlux);
    ImGui::Text("Recovery residual: %.2e (%d iterations)",
                solver.lastRecovery().finalResidual,
                solver.lastRecovery().iterations);
    ImGui::Text("Energy: %.4g", d.kineticEnergy);
    ImGui::Text("Simulated time: %.2f s | FPS: %.1f", simulationTime, ImGui::GetIO().Framerate);
    ImGui::Text("Smoke mass: %.4f | Added this step: %.6f", d.smokeMass, solver.timings().emittedMass);
    ImGui::Text("Emitter: %s | source vertices: %d", p.emitterEnabled && p.sourceStrength > 0 ? "on" : "off",
                solver.timings().sourceVertices);
    if (p.emitterEnabled && p.sourceStrength > 0 && solver.timings().sourceVertices == 0 && !paused)
      ImGui::TextWrapped("Emitter misses mesh vertices: move it into the bunny or increase its radius.");
    ImGui::Text("Boundary circulation residual: %.2e",
                solver.lastBoundaryReconstruction().finalResidual);
    ImGui::TextUnformatted("Drag to orbit; scroll to zoom");
    ImGui::End();
    if (box != solver.boxSize()) {
      solver.setBoxSize(box);
      simulationTime = 0;
      raster.update(solver);
      renderer.updateDomain(raster.renderData());
      updateWires();
      dist = box * 2.5f;
    }
    if (!ImGui::GetIO().WantCaptureMouse) {
      dist = std::clamp(dist - ImGui::GetIO().MouseWheel * dist * .1f,
                        box * .3f, box * 8);
      if (ImGui::IsMouseDown(0)) {
        double x, y;
        glfwGetCursorPos(window, &x, &y);
        if (dragging) {
          azimuth += float(x - lastX) * .005f;
          elevation =
              std::clamp(elevation + float(lastY - y) * .005f, -1.5f, 1.5f);
        }
        lastX = x;
        lastY = y;
        dragging = true;
      } else
        dragging = false;
    } else
      dragging = false;
    ImGui::Render();
    ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
    glfwSwapBuffers(window);
  }
  return 0;
}
