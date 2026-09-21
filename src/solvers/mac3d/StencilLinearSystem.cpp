#include "solvers/mac3d/StencilLinearSystem.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace {
std::size_t checkedCount(int resolution) {
    if (resolution <= 0) throw std::invalid_argument("Linear system resolution must be positive");
    const auto n = static_cast<std::size_t>(resolution);
    const auto limit = static_cast<std::size_t>(std::numeric_limits<int>::max());
    if (n > limit / n || n * n > limit / n) {
        throw std::length_error("Linear system exceeds the supported index range");
    }
    return n * n * n;
}
}

void StencilLinearSystem::resize(int newResolution) {
    const auto count = checkedCount(newResolution);
    matrix.assign(count, {});
    x.assign(count, 0.f);
    b.assign(count, 0.f);
    resolution = newResolution;
}

void StencilLinearSystem::validate() const {
    const auto count = checkedCount(resolution);
    if (matrix.size() != count || x.size() != count || b.size() != count) {
        throw std::invalid_argument("Linear system storage does not match its resolution");
    }
    const int n = resolution;
    for (int k = 0; k < n; ++k) for (int j = 0; j < n; ++j) for (int i = 0; i < n; ++i) {
        const auto idx = index(i, j, k);
        const auto& a = matrix[idx];
        if (!std::isfinite(a.diagonal) || a.diagonal < 0.f ||
            !std::isfinite(a.right) || !std::isfinite(a.up) || !std::isfinite(a.front) ||
            !std::isfinite(x[idx]) || !std::isfinite(b[idx]) ||
            (i == n - 1 && a.right != 0.f) || (j == n - 1 && a.up != 0.f) ||
            (k == n - 1 && a.front != 0.f)) {
            throw std::invalid_argument("Invalid stencil coefficient, boundary entry, or vector value");
        }
        if (a.diagonal == 0.f && (b[idx] != 0.f || a.right != 0.f || a.up != 0.f || a.front != 0.f ||
            (i > 0 && matrix[idx - 1].right != 0.f) ||
            (j > 0 && matrix[idx - n].up != 0.f) ||
            (k > 0 && matrix[idx - n * n].front != 0.f))) {
            throw std::invalid_argument("Only an isolated zero equation may have zero diagonal");
        }
    }
}

double StencilLinearSystem::residualInfinityNorm() const {
    validate();
    double result = 0.0;
    const int n = resolution;
    for (int k = 0; k < n; ++k) for (int j = 0; j < n; ++j) for (int i = 0; i < n; ++i) {
        const auto idx = index(i, j, k);
        const auto& a = matrix[idx];
        double ax = double(a.diagonal) * x[idx];
        if (i > 0) ax += double(matrix[idx - 1].right) * x[idx - 1];
        if (i + 1 < n) ax += double(a.right) * x[idx + 1];
        if (j > 0) ax += double(matrix[idx - n].up) * x[idx - n];
        if (j + 1 < n) ax += double(a.up) * x[idx + n];
        if (k > 0) ax += double(matrix[idx - n * n].front) * x[idx - n * n];
        if (k + 1 < n) ax += double(a.front) * x[idx + n * n];
        result = std::max(result, std::abs(double(b[idx]) - ax));
    }
    return result;
}
