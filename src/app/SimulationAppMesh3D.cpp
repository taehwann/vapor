#include "SimulationApp.hpp"
#include "WindowContext.hpp"
#include "graphics/GlLoader.hpp"
#include "renderer/OrbitCamera.hpp"
#include "renderer/TetrahedralWireRenderer.hpp"
#include <imgui.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_opengl3.h>
#include <limits>

int SimulationApp::runMesh3D(const RunOptions& options) {
    const auto data =
        TetrahedralWireData::load(options.assetDirectory / "examples/simplicial3d-bunny/bunny.tet", 2.f);
    WindowContext context(1100, 900, "vapor - Stanford bunny tetrahedra", options.hidden);
    auto* window = context.window();
    TetrahedralWireRenderer renderer;
    renderer.init(data);
    bool primal = true, dual = true, barycentric = false, depth = true, dragging = false, save = false;
    float azimuth = -1.2f, elevation = .3f, dist = 2.f * 2, cut = 1;
    double lastX = 0, lastY = 0;
    int frames = 0;
    std::string status;
    while (!glfwWindowShouldClose(window)) {
        glfwPollEvents();
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();
        int w, h;
        glfwGetFramebufferSize(window, &w, &h);
        glBindFramebuffer(GL_FRAMEBUFFER, 0);
        glDepthMask(GL_TRUE);
        glClearColor(.018f, .024f, .035f, 1);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        float mvp[16], cx, cy, cz;
        buildMVP(mvp, cx, cy, cz, azimuth, elevation, dist, float(w) / std::max(h, 1), 2.f);
        renderer.render(mvp, w, h, primal, dual, barycentric,
                        cut >= 1 ? std::numeric_limits<float>::max() : cut * 2.f, depth);
        if (save) {
            try {
                saveWireframeImage("bunny-wireframe.ppm", w, h);
                status = "Saved bunny-wireframe.ppm";
            } catch (const std::exception& e) {
                status = e.what();
            }
            save = false;
        }
        ImGui::SetNextWindowPos({16, 16}, ImGuiCond_FirstUseEver);
        ImGui::Begin("Stanford bunny mesh", nullptr, ImGuiWindowFlags_AlwaysAutoResize);
        ImGui::Text("%zu vertices | %zu tetrahedra", data.vertexCount, data.tetCount);
        ImGui::Text("%zu primal edges | %zu boundary faces", data.primal.size(), data.boundaryFaces);
        ImGui::Checkbox("Primal (white)", &primal);
        ImGui::Checkbox("Dual (green)", &dual);
        ImGui::Checkbox("Barycentric dual", &barycentric);
        ImGui::TextUnformatted(barycentric ? "Dual: barycenters through face centroids"
                                           : "Dual: circumcenters through face circumcenters");
        if (!barycentric)
            ImGui::TextUnformatted("Circumcenters may extend outside the bunny.");
        ImGui::SliderFloat("X cutaway", &cut, 0, 1, "%.2f");
        ImGui::Checkbox("Depth test wires", &depth);
        if (ImGui::Button("Reset view")) {
            azimuth = -1.2f;
            elevation = .3f;
            dist = 2.f * 2;
            cut = 1;
        }
        if (ImGui::Button("Save render (.ppm)"))
            save = true;
        if (!status.empty())
            ImGui::TextUnformatted(status.c_str());
        ImGui::TextUnformatted("Geometry test; fluid simulation requires box boundaries.");
        ImGui::TextUnformatted("Drag to orbit; scroll to zoom");
        ImGui::End();
        if (!ImGui::GetIO().WantCaptureMouse) {
            dist = std::clamp(dist - ImGui::GetIO().MouseWheel * dist * .1f, 2.f * .8f, 2.f * 8);
            if (ImGui::IsMouseDown(0)) {
                double x, y;
                glfwGetCursorPos(window, &x, &y);
                if (dragging) {
                    azimuth += float(x - lastX) * .005f;
                    elevation = std::clamp(elevation + float(lastY - y) * .005f, -1.5f, 1.5f);
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
        if (options.maxFrames > 0 && ++frames >= options.maxFrames)
            glfwSetWindowShouldClose(window, GLFW_TRUE);
    }
    return 0;
}
