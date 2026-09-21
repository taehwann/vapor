#pragma once

#include <cstddef>
#include <vector>

// Symmetric seven-point structured matrix. Negative-direction coefficients live
// in the corresponding neighbor's positive-direction entry. No fluid semantics.
struct StencilRow {
    float diagonal = 0.f;
    float right = 0.f;
    float up = 0.f;
    float front = 0.f;
};

struct StencilLinearSystem {
    int resolution = 0;
    std::vector<StencilRow> matrix;
    std::vector<float> x;
    std::vector<float> b;

    void resize(int newResolution);
    void validate() const;
    [[nodiscard]] std::size_t index(int i, int j, int k) const noexcept {
        return i + std::size_t(resolution) * (j + std::size_t(resolution) * k);
    }
    [[nodiscard]] double residualInfinityNorm() const;
};
