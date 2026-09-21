#include "renderer/SimplicialRenderData3D.hpp"
#include "solvers/simplicial3d/SimplicialFluidSolver3D.hpp"
#include <cmath>
#include <iostream>
#include <map>
using namespace simplicial3d;
void require(bool v, const char *why) {
  if (!v)
    throw std::runtime_error(why);
}
int main(int argc, char **argv) {
  try {
    SimplicialFluidSolver3D solver(
        Mesh::loadDomain("examples/simplicial3d-bunny/bunny-fluid.tet", 2));
    auto &p = solver.parameters();
    p.sourceStrength = 1;
    p.cgIterations = 1500;
    if (argc > 2) p.smokeDecay = std::stof(argv[2]);
    const auto &mesh = solver.domain();
    const auto &dual = solver.dual();
    require(int(mesh.vertices.size()) - int(mesh.edges.size()) +
                    int(mesh.faces.size()) - int(mesh.tets.size()) ==
                1,
            "Bunny topology changed");
    std::vector<std::vector<int>> adjacency(dual.vertices.size());
    for (auto e : dual.edges) {
      adjacency[e.a].push_back(e.b);
      adjacency[e.b].push_back(e.a);
      require(length(dual.vertices[e.a].position -
                     dual.vertices[e.b].position) > 1e-10,
              "Zero-length dual edge");
    }
    std::vector<bool> visited(adjacency.size());
    std::vector<int> queue{0};
    visited[0] = true;
    for (size_t i = 0; i < queue.size(); ++i)
      for (int v : adjacency[queue[i]])
        if (!visited[v]) {
          visited[v] = true;
          queue.push_back(v);
        }
    require(queue.size() == adjacency.size(), "Disconnected dual graph");
    for (auto loop : dual.loops) {
      std::map<int, int> balance;
      for (auto s : loop) {
        auto e = dual.edges[s.edge];
        balance[e.a] += s.sign;
        balance[e.b] -= s.sign;
      }
      for (auto [v, b] : balance)
        require(b == 0, "Open dual face");
    }
    for (auto v : dual.vertices)
      require(mesh.locate(v.position) >= 0, "Dual vertex outside bunny");
    for (auto star : {mesh.star0, mesh.star1, mesh.star2, mesh.star3})
      for (double x : star)
        require(x > 0 && std::isfinite(x), "Invalid Hodge metric");
    require(mesh.locate({0, 0, 0}) < 0, "Bunny exterior accepted");
    bool reentry = false;
    for (size_t a = 0; a < mesh.tets.size() && !reentry; a += 97)
      for (size_t b = a + 1; b < mesh.tets.size() && !reentry; b += 113) {
        Point start = mesh.tets[a].center, end = mesh.tets[b].center;
        Point hit = mesh.clipSegment(start, end);
        if (length(hit - end) > 1e-5) {
          require(mesh.locate(hit) >= 0, "First exit moved outside the domain");
          require(mesh.locate(hit + (end - start) * 1e-5) < 0,
                  "Trace stopped before the first exit");
          reentry = true;
        }
      }
    require(reentry, "No concave-domain reentry trace tested");
    int walls = 0;
    for (const auto &face : mesh.faces)
      if (face.boundary()) {
        int t = face.tets[0] >= 0 ? face.tets[0] : face.tets[1];
        Point start = mesh.tets[t].center;
        Point end = face.center + (face.center - start) * 2;
        Point hit = mesh.clipSegment(start, end);
        require(length(hit - face.center) < 1e-6,
                "Trace did not stop at surface triangle");
        ++walls;
      }
    std::cout << "PASS bunny geometry: " << mesh.tets.size() << " tets, "
              << walls << " boundary triangles, one connected dual graph\n";
    int steps = argc > 1 ? std::stoi(argv[1]) : 180;
    double maximumSpeed = 0;
    std::vector<double> initialFlux;
    std::vector<float> initialDensity;
    for (int i = 0; i < steps; ++i) {
      solver.advance(1.f / 30);
      if (i == 0) {
        initialFlux = solver.flux();
        initialDensity = solver.density();
      }
      auto d = solver.diagnostics();
      require(solver.timings().sourceVertices > 0 && solver.timings().emittedMass > 0,
              "Continuous source stopped replenishing smoke");
      maximumSpeed = std::max(maximumSpeed, d.maxSpeed);
      require(solver.lastRecovery().converged &&
                  solver.lastBoundaryReconstruction().converged,
              "Unconverged solve");
      require(std::isfinite(d.kineticEnergy) && d.maxSpeed < 100,
              "Unbounded fluid velocity");
      require(d.maxDivergence < 1e-7 && d.maxBoundaryFlux < 1e-12,
              "Divergence or wall leakage");
      for (float density : solver.density())
        require(std::isfinite(density) && density >= 0 && density <= 1.00001,
                "Invalid smoke density");
      if (i % 30 == 29)
        std::cout << "t=" << (i + 1) / 30. << " speed=" << d.maxSpeed
                  << " E=" << d.kineticEnergy << " mass=" << d.smokeMass
                  << " emitted=" << solver.timings().emittedMass
                  << " div=" << d.maxDivergence << std::endl;
    }
    require(solver.diagnostics().smokeMass > 0, "No smoke emitted");
    for (const auto &face : mesh.faces)
      if (face.boundary()) {
        Point sample =
            (mesh.vertices[face.vertices[0]] + mesh.vertices[face.vertices[1]] +
             mesh.vertices[face.vertices[2]]) *
            (1. / 3);
        require(std::abs(dot(solver.sampleVelocity(sample), face.areaNormal)) /
                        length(face.areaNormal) <
                    1e-7,
                "Nonzero wall normal velocity");
      }
    std::cout << "PASS bunny simulation: " << steps
              << " steps, max speed=" << maximumSpeed << '\n';
    SimplicialRenderData3D raster(16);
    auto checkRaster = [&] {
      raster.update(solver);
      auto image = raster.renderData();
      double size = solver.boxSize();
      for (int z = 0; z < 16; ++z)
        for (int y = 0; y < 16; ++y)
          for (int x = 0; x < 16; ++x)
            require(std::abs(image.density[x + 16 * (y + 16 * z)] -
                             solver.sampleDensity(
                                 {size * (x + .5) / 16, size * (y + .5) / 16,
                                  size * (z + .5) / 16})) < 1e-6,
                    "Cached raster differs from direct sampling");
    };
    checkRaster();
    solver.parameters().emitterEnabled = false;
    solver.advance(1.f / 30);
    require(solver.timings().emittedMass == 0, "Disabled emitter added smoke");
    solver.parameters().emitterEnabled = true;
    solver.parameters().emitterX = 0;
    solver.parameters().emitterY = 0;
    solver.parameters().emitterZ = 0;
    solver.parameters().emitterRadius = .01f;
    solver.advance(1.f / 30);
    require(solver.timings().sourceVertices == 0 && solver.timings().emittedMass == 0,
            "Exterior source emitted into the bunny");
    solver.parameters().emitterX = .5f;
    solver.parameters().emitterY = .5f;
    solver.parameters().emitterZ = .2f;
    solver.parameters().emitterRadius = .12f;
    solver.resetState();
    checkRaster();
    solver.advance(1.f / 30);
    require(initialFlux == solver.flux() && initialDensity == solver.density(),
            "Parallel reset replay changed results");
    checkRaster();
    solver.setBoxSize(3);
    require(!solver.domain().boxDomain && solver.domain().tets.size() == 9921,
            "Resize lost bunny domain");
    require(solver.diagnostics().smokeMass == 0, "Resize did not reset smoke");
    checkRaster();
    solver.advance(1.f / 30);
    require(solver.diagnostics().smokeMass > 0 &&
                solver.diagnostics().maxBoundaryFlux == 0,
            "Resized bunny failed");
    checkRaster();
    solver.parameters().buoyancy = 0;
    solver.parameters().stirStrength = 4;
    solver.resetState();
    solver.advance(1.f/30);
    const auto stirredFlux = solver.flux();
    require(solver.diagnostics().kineticEnergy > 1e-9, "Stirring did not create flow");
    solver.resetState();
    solver.advance(1.f/30);
    require(stirredFlux == solver.flux(), "Stirring reset did not replay deterministically");
    for (int i=0;i<180;++i) {
      solver.advance(1.f/30);
      auto d=solver.diagnostics();
      require(std::isfinite(d.kineticEnergy) && d.maxSpeed < 100 &&
                  d.maxDivergence < 1e-7 && d.maxBoundaryFlux == 0,
              "Stirring broke finite divergence-free wall-constrained flow");
    }
    solver.resetState();
    solver.parameters().emitterEnabled = false;
    solver.advance(1.f/30);
    require(solver.diagnostics().kineticEnergy == 0, "Disabled emitter still stirred");
    solver.parameters().emitterEnabled = true;
    solver.parameters().stirStrength = 0;
    solver.advance(1.f/30);
    require(solver.diagnostics().kineticEnergy == 0, "Zero stirring created flow");
    std::cout << "PASS stirring: drives flow, preserves walls/divergence, resets and disables\n";
  } catch (const std::exception &e) {
    std::cerr << "FAIL " << e.what() << '\n';
    return 1;
  }
}
