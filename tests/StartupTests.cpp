#include "app/SimulationApp.hpp"
#include <GLFW/glfw3.h>
#include <iostream>

int main() {
    if (!glfwInit())
        return 77;
    glfwWindowHint(GLFW_VISIBLE, GLFW_FALSE);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    auto* probe = glfwCreateWindow(32, 32, "probe", nullptr, nullptr);
    if (!probe) {
        glfwTerminate();
        return 77;
    }
    glfwDestroyWindow(probe);
    glfwTerminate();
    try {
        SimulationApp app;
        RunOptions options;
        options.assetDirectory = std::filesystem::current_path();
        options.maxFrames = 1;
        options.hidden = true;
        for (auto choice :
             {SolverChoice::Mac2DSimple, SolverChoice::Mac2DReflectionSL, SolverChoice::Mac2DReflectionMC,
              SolverChoice::Mac3DSimpleSLCPU, SolverChoice::Mac3DSimpleSLGPU, SolverChoice::Mac3DSimpleMCCPU,
              SolverChoice::Mac3DSimpleMCGPU, SolverChoice::Mac3DReflectionSLCPU,
              SolverChoice::Mac3DReflectionSLGPU, SolverChoice::Mac3DReflectionMCCPU,
              SolverChoice::Mac3DReflectionMCGPU, SolverChoice::Simplicial2DTeapot,
              SolverChoice::BunnyMesh}) {
            if (app.runSolver(choice, options) != 0)
                return 1;
        }
        std::cout << "PASS every menu selection constructs, steps, renders, and releases its context\n";
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
