#include "graphics/GlLoader.hpp"
#include "renderer/OrbitCamera.hpp"
#include "renderer/SimplicialRenderData3D.hpp"
#include "renderer/TetrahedralWireRenderer.hpp"
#include "renderer/SimplicialWireData3D.hpp"
#include "renderer/SilhouetteRenderer.hpp"
#include "renderer/VolumeRenderer.hpp"
#include "solvers/simplicial3d/SimplicialFluidSolver3D.hpp"
#include <iostream>
int main(int argc, char **argv) {
  if (!glfwInit())
    return 77;
  glfwWindowHint(GLFW_VISIBLE, GLFW_FALSE);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
  auto *window =
      glfwCreateWindow(960, 960, "Bunny validation", nullptr, nullptr);
  if (!window) {
    glfwTerminate();
    return 77;
  }
  glfwMakeContextCurrent(window);
  if (!loadOpenGLFunctions())
    return 1;
  int status = 0;
  try {
    SimplicialFluidSolver3D solver(simplicial3d::Mesh::loadDomain(
        "examples/simplicial3d-bunny/bunny-fluid.tet", 2));
    solver.parameters().sourceStrength = 1;
    solver.parameters().cgIterations = 1500;
    int steps = argc > 2 ? std::stoi(argv[2]) : 90;
    for (int i = 0; i < steps; ++i)
      solver.advance(1.f / 30);
    SimplicialRenderData3D raster(80);
    raster.update(solver);
    VolumeRenderer smoke;
    smoke.init(raster.renderData());
    auto data = simplicialWires3D(solver.domain(), solver.dual());
    auto surface = simplicialWires3D(solver.domain(), solver.dual(), true);
    size_t boundaryEdges = 0, wallEdges = 0;
    for (const auto& e : solver.domain().edges) boundaryEdges += e.boundary;
    for (const auto& e : solver.dual().edges) wallEdges += e.wall;
    if (surface.primal.size() != boundaryEdges || surface.circumcentric.size() != wallEdges ||
        surface.primal.size() >= data.primal.size() || surface.circumcentric.size() >= data.circumcentric.size())
      throw std::runtime_error("Boundary view contains interior edges");
    TetrahedralWireRenderer wires, outline;
    wires.init(data);
    outline.init(surface);
    SilhouetteRenderer silhouette;
    silhouette.init(solver.domain());
    std::vector<unsigned char> smokePixels;
    float mvp[16], cx, cy, cz;
    buildMVP(mvp, cx, cy, cz, -1.2f, .2f, 4.2f, 1, 2);
    for (int mode = 0; mode < 8; ++mode) {
      glDepthMask(GL_TRUE);
      glClearColor(.03f, .04f, .07f, 1);
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
      if (mode == 0 || mode == 7) {
        smoke.render(raster.renderData(), 960, 960, mvp, cx, cy, cz, .5f, 35,
                     .2683f, .3578f, .8944f, .3f, .1f, nullptr);
      } else if (mode == 6) silhouette.render(mvp, 960, 960);
      else if (mode >= 4)
        outline.render(mvp, 960, 960, mode == 4, mode == 5, false, 3, true);
      else
        wires.render(mvp, 960, 960, mode != 2, mode != 1, false, 3, true);
      if (mode == 7) silhouette.render(mvp, 960, 960);
      std::vector<unsigned char> pixels(960 * 960 * 3);
      glReadPixels(0, 0, 960, 960, GL_RGB, GL_UNSIGNED_BYTE, pixels.data());
      if (mode == 0) smokePixels = pixels;
      if (mode == 7) {
        size_t changed = 0;
        for (size_t i = 0; i < pixels.size(); i += 3) {
          if (pixels[i] != smokePixels[i] || pixels[i+1] != smokePixels[i+1] || pixels[i+2] != smokePixels[i+2]) {
            ++changed;
            if (pixels[i] != 255 || pixels[i+1] != 255 || pixels[i+2] != 255)
              throw std::runtime_error("Outline changed smoke pixels instead of only drawing its contour");
          }
        }
        if (changed < 100 || changed > 12000) throw std::runtime_error("Outline missing or filling surface");
      }
      size_t white = 0, green = 0, lit = 0;
      for (size_t i = 0; i < pixels.size(); i += 3) {
        if (pixels[i] > 230 && pixels[i + 1] > 230 && pixels[i + 2] > 230)
          ++white;
        if (pixels[i] < 80 && pixels[i + 1] > 200 && pixels[i + 2] < 100)
          ++green;
        if (pixels[i] > 40)
          ++lit;
      }
      if (glGetError() != GL_NO_ERROR || (mode == 0 && lit < 200) ||
          (mode == 1 && (white < 1000 || green)) ||
          (mode == 2 && (green < 1000 || white)) ||
          (mode == 3 && (white < 500 || green < 500)) ||
          (mode == 4 && (white < 500 || green)) ||
          (mode == 5 && (green < 500 || white)) ||
          (mode == 6 && (white < 100 || white > 12000 || green)))
        throw std::runtime_error("Bunny render check failed");
      if (mode == 0)
        silhouette.render(mvp, 960, 960);
      if (argc > 1)
        saveWireframeImage(std::string(argv[1]) + "-" + std::to_string(mode) +
                               ".ppm",
                           960, 960);
      std::cout << "PASS bunny render " << mode << ": lit=" << lit
                << " white=" << white << " green=" << green << '\n';
    }
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    status = 1;
  }
  glfwDestroyWindow(window);
  glfwTerminate();
  return status;
}
