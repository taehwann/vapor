#include "solvers/mac2d/MacGridFluidSolver2D.hpp"
#include "renderer/MacRenderData2D.hpp"
#include "renderer/Renderer2D.hpp"
#include "graphics/GlLoader.hpp"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {
void require(bool condition, const char* message) { if (!condition) throw std::runtime_error(message); }
}
int main(int argc, char** argv) {
    if (!glfwInit()) { std::cout << "SKIP GLFW unavailable\n"; return 77; }
    glfwWindowHint(GLFW_VISIBLE, GLFW_FALSE);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    GLFWwindow* window = glfwCreateWindow(512, 512, "Vapor 2D graphics test", nullptr, nullptr);
    if (!window) { glfwTerminate(); std::cout << "SKIP OpenGL unavailable\n"; return 77; }
    glfwMakeContextCurrent(window);
    if (!loadOpenGLFunctions()) { glfwDestroyWindow(window); glfwTerminate(); return 1; }
    int status = 0;
    Renderer2D renderer;
    try {
        int w, h;
        glfwGetFramebufferSize(window, &w, &h);
        std::vector<unsigned char> pixels(std::size_t(w) * h * 3);
        auto draw = [&](const ImageRenderData& data, float exposure = 3.f) {
            glViewport(0, 0, w, h);
            glClearColor(0.f, 0.f, 0.f, 1.f);
            glClear(GL_COLOR_BUFFER_BIT);
            renderer.render(data, w, h, exposure);
            glReadBuffer(GL_BACK);
            glPixelStorei(GL_PACK_ALIGNMENT, 1);
            glReadPixels(0, 0, w, h, GL_RGB, GL_UNSIGNED_BYTE, pixels.data());
            require(glGetError() == GL_NO_ERROR, "OpenGL error in 2D rendering");
        };
        auto blueAt = [&](int x, int y) { return pixels[3 * (x + w * y) + 2]; };
        MacGridState2D state(8);
        for (int y = 0; y < 4; ++y) for (int x = 0; x < 4; ++x) state.density()[state.idC(x, y)] = 1.f;

        renderer.init(imageData(state));
        draw(imageData(state));
        require(blueAt(w / 4, h / 4) > 200, "2D density is not visible in the lower-left quadrant");
        require(blueAt(w / 4, 3 * h / 4) < 5 && blueAt(3 * w / 4, h / 4) < 5, "2D image orientation or sampling is wrong");

        MacGridState2D resized(16);
        for (int y = 8; y < 16; ++y) for (int x = 8; x < 16; ++x) resized.density()[resized.idC(x, y)] = 1.f;

        draw(imageData(resized));
        require(blueAt(3 * w / 4, 3 * h / 4) > 200 && blueAt(w / 4, h / 4) < 5, "Texture resize retained stale data");
        draw(imageData(resized), 0.f);
        require(*std::max_element(pixels.begin(), pixels.end()) < 5, "Zero exposure should produce a black image");
        auto historical = imageData(resized);
        historical.historicalPalette = true;
        draw(historical);
        require(std::abs(int(blueAt(w/4,h/4))-18)<=1 &&
                std::abs(int(blueAt(3*w/4,3*h/4))-222)<=1, "Historical palette differs from ba20987");

        std::fill(resized.density().begin(), resized.density().end(), 1.f);
        auto wide = imageData(resized);
        wide.worldWidth = 2.f; wide.worldHeight = 1.f;
        draw(wide);
        require(blueAt(w / 2, h / 2) > 200 && blueAt(w / 2, h / 8) < 5, "2D physical aspect ratio is not preserved");
        auto invalid = wide;
        invalid.density = invalid.density.first(1);
        bool rejected = false;
        try { renderer.render(invalid, w, h); } catch (const std::invalid_argument&) { rejected = true; }
        require(rejected, "Invalid image extent was accepted");

        // Real solver -> borrowed adapter -> texture upload -> pixels.
        MacGridFluidSolver2D solver(32);
        auto& fluid = solver;
        require(fluid.domain().dimension() == Dimension::D2, "2D fluid dimension");
        for (int i = 0; i < 120; ++i) fluid.advance(1.f / 60.f);

        renderer.shutdown();
        renderer.shutdown();
        renderer.init(imageData(solver.state()));
        draw(imageData(solver.state()));
        const auto lit = std::count_if(pixels.begin(), pixels.end(), [](unsigned char v) { return v > 32; });
        require(lit > std::ptrdiff_t(pixels.size() / 100), "2D smoke solver rendered no substantial plume");
        if (argc > 1) {
            std::ofstream capture(argv[1], std::ios::binary);
            capture << "P6\n" << w << ' ' << h << "\n255\n";
            for (int y = h - 1; y >= 0; --y)
                capture.write(reinterpret_cast<const char*>(pixels.data() + std::size_t(y) * w * 3), w * 3);
            require(bool(capture), "Could not save 2D graphics capture");
        }
        std::cout << "PASS 2D image orientation, resize, exposure, aspect ratio, and simulated smoke pixels\n";
    } catch (const std::exception& error) {
        std::cerr << "FAIL 2D graphics: " << error.what() << '\n'; status = 1;
    }
    renderer.shutdown();
    glfwDestroyWindow(window);
    glfwTerminate();
    return status;
}
