#include "renderer/SimplicialRenderData3D.hpp"
#include "solvers/simplicial3d/SimplicialFluidSolver3D.hpp"
#include <chrono>
#include <iostream>
int main() {
  SimplicialFluidSolver3D s(simplicial3d::Mesh::loadDomain(
      "examples/simplicial3d-bunny/bunny-fluid.tet", 2));
  SimplicialRenderData3D raster(64);
  for (int i = 0; i < 40; ++i) {
    s.advance(1.f / 30);
    auto begin = std::chrono::steady_clock::now();
    raster.update(s);
    double rendering = std::chrono::duration<double, std::milli>(
                           std::chrono::steady_clock::now() - begin)
                           .count();
    const auto &t = s.timings();
    if (i % 5 == 0)
      std::cout << i << " sub=" << t.substeps << " emit=" << t.emissionMs
                << " adv=" << t.advectionMs << " force=" << t.forcesMs
                << " recovery=" << t.recoveryMs << " smoke=" << t.smokeMs
                << " raster=" << rendering
                << " cg=" << s.lastRecovery().iterations
                << " boundary=" << s.lastBoundaryReconstruction().iterations
                << std::endl;
  }
}
