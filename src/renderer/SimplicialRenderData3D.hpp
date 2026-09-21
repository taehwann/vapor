#pragma once
#include "RenderData.hpp"
#include "solvers/simplicial3d/SimplicialFluidSolver3D.hpp"
#include <algorithm>
#include <stdexcept>

// Rendering-only resampling. No Cartesian velocity or pressure participates in
// the tetrahedral simulation. Owns the buffer borrowed by RenderData.
class SimplicialRenderData3D {
public:
  explicit SimplicialRenderData3D(int resolution = 48) : n_(resolution) {
    if (n_ < 1 || n_ > 512)
      throw std::invalid_argument("Invalid volume resolution");
    density_.resize(size_t(n_) * n_ * n_);
  }
  void update(const SimplicialFluidSolver3D &solver) {
    const auto &mesh = solver.domain();
    if (mesh.tets.data() != meshKey_ || size_ != solver.boxSize()) {
      size_ = solver.boxSize();
      meshKey_ = mesh.tets.data();
      samples_.clear();
      std::fill(density_.begin(), density_.end(), 0.f);
      for (int z = 0; z < n_; ++z)
        for (int y = 0; y < n_; ++y)
          for (int x = 0; x < n_; ++x) {
            simplicial3d::Point p{size_ * (x + .5) / n_, size_ * (y + .5) / n_,
                                  size_ * (z + .5) / n_};
            int ti = mesh.locate(p);
            if (ti < 0)
              continue;
            auto w = mesh.barycentric(ti, p);
            double total = 0;
            for (auto &value : w) {
              value = std::max(0., value);
              total += value;
            }
            for (auto &value : w)
              value /= total;
            samples_.push_back(
                {x + n_ * (y + n_ * z), mesh.tets[ti].vertices, w});
          }
    }
    const auto &values = solver.vertexDensity();
    for (const auto &sample : samples_) {
      double value = 0;
      for (int k = 0; k < 4; ++k)
        value += sample.weights[k] * values[sample.vertices[k]];
      density_[sample.pixel] = float(value);
    }
  }
  RenderData renderData() const { return {density_, n_, n_, n_, size_}; }

private:
  int n_;
  float size_ = 1;
  std::vector<float> density_;
  struct Sample {
    int pixel;
    std::array<int, 4> vertices;
    std::array<double, 4> weights;
  };
  const simplicial3d::Tet *meshKey_ = nullptr;
  std::vector<Sample> samples_;
};
