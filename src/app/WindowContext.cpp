#include "app/WindowContext.hpp"
#include "graphics/GlLoader.hpp"
#include <imgui.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_opengl3.h>
#include <stdexcept>

WindowContext::WindowContext(int width, int height, const char* title, bool hidden) {
    try {
        if (!(glfw_ = glfwInit() != 0))
            throw std::runtime_error("GLFW initialization failed");
        glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
        glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
        glfwWindowHint(GLFW_VISIBLE, hidden ? GLFW_FALSE : GLFW_TRUE);
        window_ = glfwCreateWindow(width, height, title, nullptr, nullptr);
        if (!window_)
            throw std::runtime_error("OpenGL 4.3 window creation failed");
        glfwMakeContextCurrent(window_);
        glfwSwapInterval(1);
        if (!loadOpenGLFunctions())
            throw std::runtime_error("OpenGL entry point loading failed");
        IMGUI_CHECKVERSION();
        ImGui::CreateContext();
        imgui_ = true;
        ImGui::StyleColorsDark();
        if (!(platform_ = ImGui_ImplGlfw_InitForOpenGL(window_, true)))
            throw std::runtime_error("ImGui GLFW initialization failed");
        if (!(renderer_ = ImGui_ImplOpenGL3_Init("#version 330")))
            throw std::runtime_error("ImGui OpenGL initialization failed");
    } catch (...) {
        shutdown();
        throw;
    }
}
WindowContext::~WindowContext() { shutdown(); }
void WindowContext::shutdown() noexcept {
    if (window_)
        glfwMakeContextCurrent(window_);
    if (renderer_)
        ImGui_ImplOpenGL3_Shutdown();
    if (platform_)
        ImGui_ImplGlfw_Shutdown();
    if (imgui_)
        ImGui::DestroyContext();
    if (window_)
        glfwDestroyWindow(window_);
    if (glfw_)
        glfwTerminate();
}
