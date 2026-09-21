#pragma once
struct GLFWwindow;

// Declare before renderers and GPU solvers so the context is destroyed last.
class WindowContext {
  public:
    WindowContext(int width, int height, const char* title, bool hidden = false);
    ~WindowContext();
    WindowContext(const WindowContext&) = delete;
    WindowContext& operator=(const WindowContext&) = delete;
    GLFWwindow* window() const noexcept { return window_; }

  private:
    void shutdown() noexcept;
    GLFWwindow* window_ = nullptr;
    bool glfw_ = false, imgui_ = false, platform_ = false, renderer_ = false;
};
