#include "app/SimulationApp.hpp"
#include "app/WindowContext.hpp"
#include "graphics/GlLoader.hpp"
#include <imgui.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_opengl3.h>
#include <optional>
#include <stdexcept>

namespace {
struct SolverItem {
    const char* name;
    SolverChoice choice;
};
constexpr SolverItem solvers[] = {{"Mac2DSimple-SemiLagrangian", SolverChoice::Mac2DSimple},
                                  {"Mac2DReflection-SemiLagrangian", SolverChoice::Mac2DReflectionSL},
                                  {"Mac2DReflection-MacCormack", SolverChoice::Mac2DReflectionMC},
                                  {"Mac3DSimple-SemiLagrangian-CPU", SolverChoice::Mac3DSimpleSLCPU},
                                  {"Mac3DSimple-SemiLagrangian-GPU", SolverChoice::Mac3DSimpleSLGPU},
                                  {"Mac3DSimple-MacCormack-CPU", SolverChoice::Mac3DSimpleMCCPU},
                                  {"Mac3DSimple-MacCormack-GPU", SolverChoice::Mac3DSimpleMCGPU},
                                  {"Mac3DReflection-SemiLagrangian-CPU", SolverChoice::Mac3DReflectionSLCPU},
                                  {"Mac3DReflection-SemiLagrangian-GPU", SolverChoice::Mac3DReflectionSLGPU},
                                  {"Mac3DReflection-MacCormack-CPU", SolverChoice::Mac3DReflectionMCCPU},
                                  {"Mac3DReflection-MacCormack-GPU", SolverChoice::Mac3DReflectionMCGPU},
                                  {"Simplicial2D-Teapot", SolverChoice::Simplicial2DTeapot},
                                  {"Simplicial3D-Bunny", SolverChoice::BunnyMesh}};

std::optional<SolverChoice> selectSolver(int& selected) {
    WindowContext context(720, 480, "Vapor - choose a simulation");
    auto* window = context.window();
    while (!glfwWindowShouldClose(window)) {
        glfwPollEvents();
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();
        ImGui::SetNextWindowPos({0, 0});
        ImGui::SetNextWindowSize(ImGui::GetIO().DisplaySize);
        ImGui::Begin("Choose a simulation", nullptr,
                     ImGuiWindowFlags_NoDecoration | ImGuiWindowFlags_NoMove |
                         ImGuiWindowFlags_NoSavedSettings);
        ImGui::TextUnformatted("Vapor");
        ImGui::Spacing();
        ImGui::TextUnformatted("Choose a solver, then start the simulation.");
        ImGui::SetNextItemWidth(-1);
        if (ImGui::BeginCombo("##solver", solvers[selected].name)) {
            for (int i = 0; i < int(sizeof(solvers) / sizeof(solvers[0])); ++i) {
                if (ImGui::Selectable(solvers[i].name, i == selected))
                    selected = i;
                if (i == selected)
                    ImGui::SetItemDefaultFocus();
            }
            ImGui::EndCombo();
        }
        ImGui::Spacing();
        if (solvers[selected].choice == SolverChoice::BunnyMesh)
            ImGui::TextWrapped("Inviscid smoke inside a tetrahedral bunny domain, with primal and dual wireframe controls.");
        else if (solvers[selected].choice == SolverChoice::Simplicial2DTeapot)
            ImGui::TextWrapped("Inviscid smoke in a triangular teapot domain, with solid outer and handle walls.");
        else
            ImGui::TextWrapped("Uses the selected solver's built-in scene and numerical defaults.");
        ImGui::Spacing();
        const bool start = ImGui::Button("Start simulation", {180, 36});
        ImGui::TextWrapped("Close the simulation window to return here.");
        ImGui::End();
        ImGui::Render();
        int width, height;
        glfwGetFramebufferSize(window, &width, &height);
        glViewport(0, 0, width, height);
        glClearColor(.025f, .03f, .045f, 1.f);
        glClear(GL_COLOR_BUFFER_BIT);
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
        glfwSwapBuffers(window);
        if (start)
            return solvers[selected].choice;
    }
    return {};
}
} // namespace

int SimulationApp::run(const std::filesystem::path& assetDirectory) {
    RunOptions options;
    options.assetDirectory = assetDirectory;
    int selected = 0;
    while (auto choice = selectSolver(selected)) {
        const int status = runSolver(*choice, options);
        if (status != 0)
            return status;
    }
    return 0;
}

int SimulationApp::runSolver(SolverChoice choice, const RunOptions& options) {
    switch (choice) {
    case SolverChoice::Mac2DSimple:
    case SolverChoice::Mac2DReflectionSL:
    case SolverChoice::Mac2DReflectionMC:
        return run2D(choice, options);
    case SolverChoice::Simplicial2DTeapot: {
        auto timing = options;
        timing.fixedDt = 1.f / 30.f;
        return run2D(choice, timing);
    }
    case SolverChoice::BunnyMesh:
        return runSimplicial3D(options);
    case SolverChoice::Mac3DSimpleSLCPU:
        return run3D(MacAlgorithm::SimpleSemiLagrangian, false, options);
    case SolverChoice::Mac3DSimpleSLGPU:
        return run3D(MacAlgorithm::SimpleSemiLagrangian, true, options);
    case SolverChoice::Mac3DSimpleMCCPU:
        return run3D(MacAlgorithm::SimpleMacCormack, false, options);
    case SolverChoice::Mac3DSimpleMCGPU:
        return run3D(MacAlgorithm::SimpleMacCormack, true, options);
    case SolverChoice::Mac3DReflectionSLCPU:
        return run3D(MacAlgorithm::ReflectionSemiLagrangian, false, options);
    case SolverChoice::Mac3DReflectionSLGPU:
        return run3D(MacAlgorithm::ReflectionSemiLagrangian, true, options);
    case SolverChoice::Mac3DReflectionMCCPU:
        return run3D(MacAlgorithm::ReflectionMacCormack, false, options);
    case SolverChoice::Mac3DReflectionMCGPU:
        return run3D(MacAlgorithm::ReflectionMacCormack, true, options);
    }
    throw std::invalid_argument("Unknown solver selection");
}
