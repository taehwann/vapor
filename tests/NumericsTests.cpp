#include "fluid/MacGridState.hpp"
#include "fluid/MacGridOperators.hpp"
#include "fluid/projection/CpuPressureSolver.hpp"
#include "numerics/SOR.hpp"

#include <algorithm>
#include <cmath>
#include <exception>
#include <iostream>
#include <limits>
#include <numeric>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

namespace {
// These checks remain active in Release builds (unlike assert).
void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

void near(double actual, double expected, double tolerance, const char* message) {
    if (!std::isfinite(actual) || std::abs(actual - expected) > tolerance) {
        throw std::runtime_error(std::string(message) + ": actual=" + std::to_string(actual) +
                                 ", expected=" + std::to_string(expected));
    }
}

template<class Exception = std::invalid_argument, class Function>
void rejects(Function&& function) {
    try { function(); }
    catch (const Exception&) { return; }
    throw std::runtime_error("Expected invalid input to be rejected");
}

void same(std::span<const float> actual, std::span<const float> expected, float tolerance) {
    require(actual.size() == expected.size(), "Field size mismatch");
    for (std::size_t i = 0; i < actual.size(); ++i) {
        near(actual[i], expected[i], tolerance, "Field value mismatch");
    }
}

void sameVelocities(const MacGridState& a, const MacGridState& b, float tolerance) {
    same(a.velocityX(), b.velocityX(), tolerance);
    same(a.velocityY(), b.velocityY(), tolerance);
    same(a.velocityZ(), b.velocityZ(), tolerance);
}

void seed(MacGridState& state) {
    int component = 1;
    for (auto field : {state.velocityX(), state.velocityY(), state.velocityZ()}) {
        for (std::size_t i = 0; i < field.size(); ++i) {
            field[i] = std::sin(float(i + 1) * (0.13f * component));
        }
        ++component;
    }
    for (std::size_t i = 0; i < state.cellCount(); ++i) state.density()[i] = float(i % 11) / 11.f;
}

ProjectionOptions accurateProjection(float dt = 0.02f) {
    ProjectionOptions options;
    options.dt = dt;
    options.linearSolve.maxIterations = 4000;
    options.linearSolve.absoluteTolerance = 1e-5;
    options.linearSolve.relativeTolerance = 1e-6;
    return options;
}

void gridLayout() {
    static_assert(std::is_same_v<decltype(std::declval<const MacGridState&>().density()),
                                 std::span<const float>>);
    static_assert(std::is_same_v<decltype(std::declval<const MacGridState&>().velocityX()),
                                 std::span<const float>>);
    for (int n : {1, 2, 7}) {
        MacGridState state(n, 3.5f);
        require(state.cellCount() == std::size_t(n * n * n), "Cell count");
        near(state.cellSize(), 3.5 / n, 1e-6, "Cell spacing");
        for (int axis = 0; axis < 4; ++axis) {
            const auto size = axis == 0 ? state.cellCount() : std::size_t(n * n * (n + 1));
            const auto field = axis == 0 ? state.density() : axis == 1 ? state.velocityX() :
                               axis == 2 ? state.velocityY() : state.velocityZ();
            require(field.size() == size, "Staggered field size");
            std::vector<int> visits(size, 0);
            for (int z = 0; z < n + (axis == 3); ++z)
                for (int y = 0; y < n + (axis == 2); ++y)
                    for (int x = 0; x < n + (axis == 1); ++x) {
                        const auto i = axis == 0 ? state.idC(x,y,z) : axis == 1 ? state.idX(x,y,z) :
                                       axis == 2 ? state.idY(x,y,z) : state.idZ(x,y,z);
                        require(i < size, "Grid index out of range");
                        ++visits[i];
                    }
            require(std::all_of(visits.begin(), visits.end(), [](int v) { return v == 1; }),
                    "Grid indexing must visit every entry exactly once");
        }
    }
}

void gridValidation() {
    rejects([] { MacGridState state(0); });
    rejects([] { MacGridState state(-2); });
    rejects<std::length_error>([] { MacGridState state(std::numeric_limits<int>::max()); });
    rejects([] { MacGridState state(2, 0.f); });
    rejects([] { MacGridState state(2, std::numeric_limits<float>::infinity()); });
    MacGridState state(3);
    seed(state);
    const auto original = state;
    state.setBoxSize(state.boxSize());
    sameVelocities(state, original, 0.f);
    rejects([&] { state.setBoxSize(-1.f); });
    rejects([&] { state.setBoxSize(std::numeric_limits<float>::quiet_NaN()); });
    sameVelocities(state, original, 0.f);
    state.setBoxSize(2.f);
    near(state.boxSize(), 2.f, 0, "Updated domain");
    for (auto field : {state.density(), state.velocityX(), state.velocityY(), state.velocityZ()}) {
        require(std::all_of(field.begin(), field.end(), [](float v) { return v == 0.f; }),
                "Changing domain must clear state");
    }
    seed(state);
    state.reset();
    sameVelocities(state, MacGridState(3, 2.f), 0.f);
    require(std::all_of(state.density().begin(), state.density().end(), [](float v) { return v == 0.f; }),
            "Reset must clear density");
}

void closedBoundaries() {
    MacGridState state(4);
    std::fill(state.velocityX().begin(), state.velocityX().end(), 1.f);
    std::fill(state.velocityY().begin(), state.velocityY().end(), 2.f);
    std::fill(state.velocityZ().begin(), state.velocityZ().end(), 3.f);
    MacGridOperators::enforceClosedBoundaries(state);
    const int n = state.resolution();
    for (int z = 0; z < n; ++z) for (int y = 0; y < n; ++y) for (int x = 0; x <= n; ++x)
        near(state.velocityX()[state.idX(x,y,z)], x == 0 || x == n ? 0.f : 1.f, 0, "X wall");
    for (int z = 0; z < n; ++z) for (int y = 0; y <= n; ++y) for (int x = 0; x < n; ++x)
        near(state.velocityY()[state.idY(x,y,z)], y == 0 || y == n ? 0.f : 2.f, 0, "Y wall");
    for (int z = 0; z <= n; ++z) for (int y = 0; y < n; ++y) for (int x = 0; x < n; ++x)
        near(state.velocityZ()[state.idZ(x,y,z)], z == 0 || z == n ? 0.f : 3.f, 0, "Z wall");
    std::vector<float> divergence(state.cellCount());
    MacGridOperators::computeDivergence(state, divergence);
    near(std::accumulate(divergence.begin(), divergence.end(), 0.0), 0.0, 1e-5, "Closed flux sum");
    rejects([&] { MacGridOperators::computeDivergence(state, {}); });
    rejects([&] { MacGridOperators::applyPressureGradient(state, {}, 0.1f); });
    rejects([&] { MacGridOperators::applyPressureGradient(state, divergence, 0.f); });
}

StencilLinearSystem manufactured(std::vector<float>& expected) {
    StencilLinearSystem system;
    constexpr int n = 5;
    system.resize(n);
    expected.resize(system.x.size());
    for (int z = 0; z < n; ++z) for (int y = 0; y < n; ++y) for (int x = 0; x < n; ++x)
        expected[system.index(x,y,z)] = float(1 + x + 2*y + 3*z) * 0.125f;
    for (int z = 0; z < n; ++z) for (int y = 0; y < n; ++y) for (int x = 0; x < n; ++x) {
        const auto i = system.index(x,y,z);
        system.matrix[i] = {6.f, x+1 < n ? -1.f : 0.f, y+1 < n ? -1.f : 0.f, z+1 < n ? -1.f : 0.f};
        float rhs = 6.f * expected[i];
        if (x > 0) rhs -= expected[system.index(x-1,y,z)];
        if (x+1 < n) rhs -= expected[system.index(x+1,y,z)];
        if (y > 0) rhs -= expected[system.index(x,y-1,z)];
        if (y+1 < n) rhs -= expected[system.index(x,y+1,z)];
        if (z > 0) rhs -= expected[system.index(x,y,z-1)];
        if (z+1 < n) rhs -= expected[system.index(x,y,z+1)];
        system.b[i] = rhs;
    }
    return system;
}

void sorManufactured() {
    std::vector<float> expected;
    auto system = manufactured(expected);
    SOR solver(1.35f);
    LinearSolveOptions options;
    options.relativeTolerance = 1e-6;
    const auto result = solver.solve(system, options);
    require(result.converged, "Manufactured system must converge");
    require(result.iterations > 0 && result.iterations < options.maxIterations, "Residual stopping");
    same(system.x, expected, 2e-5f);
    near(result.finalResidual, system.residualInfinityNorm(), 0.0, "Reported residual");
}

void sorOptions() {
    std::vector<float> expected;
    auto system = manufactured(expected);
    SOR solver;
    LinearSolveOptions options;
    options.maxIterations = 0;
    const auto result = solver.solve(system, options);
    require(!result.converged && result.iterations == 0, "Zero budget must not report convergence");
    near(result.initialResidual, result.finalResidual, 0.0, "Zero budget residual");
    system.x = expected;
    options.maxIterations = 3;
    const auto alreadySolved = solver.solve(system, options);
    require(alreadySolved.converged && alreadySolved.iterations == 0, "Already solved early exit");
    options.fixedIterations = true;
    const auto fixed = solver.solve(system, options);
    require(fixed.iterations == 3 && fixed.converged, "Fixed iteration budget");
}

void sorValidation() {
    rejects([] { SOR solver(0.f); });
    rejects([] { SOR solver(2.f); });
    rejects([] { SOR solver(std::numeric_limits<float>::quiet_NaN()); });
    rejects([] { SOR solver(1.f, -1.f); });
    StencilLinearSystem system;
    rejects([&] { system.resize(0); });
    rejects<std::length_error>([&] { system.resize(std::numeric_limits<int>::max()); });
    system.resize(1);
    SOR solver;
    LinearSolveOptions options;
    require(solver.solve(system, options).converged, "Isolated zero equation");
    system.b[0] = 1.f;
    rejects([&] { solver.solve(system, options); });
    system.matrix[0].diagonal = 1.f;
    options.residualCheckInterval = 0;
    rejects([&] { solver.solve(system, options); });
    options = {};
    options.maxIterations = -1;
    rejects([&] { solver.solve(system, options); });
    options = {};
    options.relativeTolerance = -1;
    rejects([&] { solver.solve(system, options); });
    options = {};
    system.matrix[0].right = -1.f;
    rejects([&] { solver.solve(system, options); });
    system.matrix[0].right = 0.f;
    system.x[0] = std::numeric_limits<float>::infinity();
    rejects([&] { solver.solve(system, options); });
    system.x.clear();
    rejects([&] { solver.solve(system, options); });
}

void projectionRest() {
    SOR solver;
    CpuPressureSolver projection(solver);
    for (int n : {1, 5, 2}) { // Reuse on differently sized grids exercises workspace rebuilding.
        MacGridState state(n);
        auto options = accurateProjection();
        const auto result = projection.project(state, options);
        require(result.linearSolve.converged && result.linearSolve.iterations == 0, "Rest projection");
        near(result.divergenceAfter, 0, 0, "Rest divergence");
        require(projection.pressure().size() == state.cellCount(), "Pressure workspace size");
        sameVelocities(state, MacGridState(n), 0.f);
        options.dt = 0.f;
        rejects([&] { projection.project(state, options); });
        options.dt = std::numeric_limits<float>::quiet_NaN();
        rejects([&] { projection.project(state, options); });
        options.dt = std::numeric_limits<float>::denorm_min();
        rejects([&] { projection.project(state, options); });
    }
}

void projectionGradient() {
    MacGridState state(8, 4.f);
    const int n = state.resolution();
    const auto options = accurateProjection(0.1f);
    std::vector<float> potential(state.cellCount());
    for (int z = 0; z < n; ++z) for (int y = 0; y < n; ++y) for (int x = 0; x < n; ++x)
        potential[state.idC(x,y,z)] = std::sin(0.3f*x) * std::cos(0.4f*y) + 0.2f*z*z/(n*n);
    // Build the input gradient independently from the production gradient routine.
    const float factor = options.dt / state.cellSize();
    for (int z = 0; z < n; ++z) for (int y = 0; y < n; ++y) for (int x = 1; x < n; ++x)
        state.velocityX()[state.idX(x,y,z)] = factor * (potential[state.idC(x,y,z)] - potential[state.idC(x-1,y,z)]);
    for (int z = 0; z < n; ++z) for (int y = 1; y < n; ++y) for (int x = 0; x < n; ++x)
        state.velocityY()[state.idY(x,y,z)] = factor * (potential[state.idC(x,y,z)] - potential[state.idC(x,y-1,z)]);
    for (int z = 1; z < n; ++z) for (int y = 0; y < n; ++y) for (int x = 0; x < n; ++x)
        state.velocityZ()[state.idZ(x,y,z)] = factor * (potential[state.idC(x,y,z)] - potential[state.idC(x,y,z-1)]);
    SOR solver(1.5f);
    CpuPressureSolver projection(solver);
    const auto result = projection.project(state, options);
    require(result.linearSolve.converged, "Gradient projection must converge");
    sameVelocities(state, MacGridState(n, 4.f), 1e-5f);
    const double mean = std::accumulate(potential.begin(), potential.end(), 0.0) / potential.size();
    for (std::size_t i = 0; i < potential.size(); ++i)
        near(projection.pressure()[i], potential[i] - mean, 1e-4, "Recovered pressure up to constant");
}

void projectionDivergence() {
    SOR solver(1.4f);
    CpuPressureSolver projection(solver);
    for (float box : {4.5f, 1.5f}) {
        MacGridState state(10, box);
        seed(state);
        const auto original = state;
        const auto result = projection.project(state, accurateProjection());
        require(result.linearSolve.converged, "Pressure solve convergence");
        require(result.divergenceBefore > 1.f, "Nontrivial input divergence");
        require(result.divergenceAfter < result.divergenceBefore * 2e-5f, "Projection divergence reduction");
        near(result.divergenceAfter, MacGridOperators::maxDivergence(state), 0, "Reported divergence");
        same(state.density(), original.density(), 0.f);
        const double mean = std::accumulate(projection.pressure().begin(), projection.pressure().end(), 0.0) /
                            projection.pressure().size();
        near(mean, 0, 1e-5, "Pressure zero-mean gauge");
        const auto projected = state;
        MacGridOperators::enforceClosedBoundaries(state);
        sameVelocities(state, projected, 0.f);
        projection.reset();
        require(std::all_of(projection.pressure().begin(), projection.pressure().end(),
                            [](float p) { return p == 0.f; }), "Reset pressure workspace");
        require(std::all_of(projection.inputDivergence().begin(), projection.inputDivergence().end(),
                            [](float d) { return d == 0.f; }), "Reset divergence workspace");
        state.setBoxSize(box + 1.f);
        const auto resetResult = projection.project(state, accurateProjection());
        near(resetResult.divergenceAfter, 0, 0, "Rest after domain change");
    }
}

void projectionSolenoidal() {
    constexpr int n = 4;
    MacGridState state(n, float(n));
    auto psi = [](int x, int y) { return float(x * (n-x) * y * (n-y)); };
    for (int z = 0; z < n; ++z) for (int y = 0; y < n; ++y) for (int x = 0; x <= n; ++x)
        state.velocityX()[state.idX(x,y,z)] = psi(x,y+1) - psi(x,y);
    for (int z = 0; z < n; ++z) for (int y = 0; y <= n; ++y) for (int x = 0; x < n; ++x)
        state.velocityY()[state.idY(x,y,z)] = -(psi(x+1,y) - psi(x,y));
    const auto original = state;
    near(MacGridOperators::maxDivergence(state), 0, 0, "Discrete curl is divergence-free");
    SOR solver;
    CpuPressureSolver projection(solver);
    const auto result = projection.project(state, accurateProjection());
    require(result.linearSolve.converged && result.linearSolve.iterations == 0, "Solenoidal early exit");
    sameVelocities(state, original, 0.f);
}

void projectionInjection() {
    static_assert(!std::is_constructible_v<CpuPressureSolver, SOR&&>);
    struct RecordingSolver final : ILinearSolver {
        mutable int calls = 0;
        LinearSolveResult solve(StencilLinearSystem& system, const LinearSolveOptions& options) const override {
            ++calls;
            const int n = system.resolution;
            near(system.matrix[system.index(0,0,0)].diagonal, 3, 0, "Neumann corner diagonal");
            near(system.matrix[system.index(1,1,1)].diagonal, 6, 0, "Interior diagonal");
            near(system.matrix[system.index(n-1,1,1)].right, 0, 0, "No external coupling");
            near(std::accumulate(system.b.begin(), system.b.end(), 0.0), 0, 1e-4, "Compatible RHS");
            return SOR(1.3f).solve(system, options);
        }
    } solver;
    MacGridState state(4);
    seed(state);
    CpuPressureSolver projection(solver);
    IPressureSolver<MacGridState>& interface = projection;
    require(interface.project(state, accurateProjection()).linearSolve.converged, "Injected solver convergence");
    require(solver.calls == 1, "Projection must call injected linear solver");
}

// Independent reference for the original main.cpp CPU pressure loop. Boundary
// neighbors use lagged p_current ghost values and the relaxation denominator is 6.
void legacyProjection(MacGridState& state, float dt, int iterations, float omega) {
    MacGridOperators::enforceClosedBoundaries(state);
    const int n = state.resolution();
    const float h = state.cellSize(), hInv = 1.f / h;
    std::vector<float> div(state.cellCount()), p(state.cellCount(), 0.f);
    for (int z = 0; z < n; ++z) for (int y = 0; y < n; ++y) for (int x = 0; x < n; ++x)
        div[state.idC(x,y,z)] = (state.velocityX()[state.idX(x+1,y,z)] - state.velocityX()[state.idX(x,y,z)] +
            state.velocityY()[state.idY(x,y+1,z)] - state.velocityY()[state.idY(x,y,z)] +
            state.velocityZ()[state.idZ(x,y,z+1)] - state.velocityZ()[state.idZ(x,y,z)]) * hInv;
    for (int iter = 0; iter < iterations; ++iter) for (int parity = 0; parity < 2; ++parity)
        for (int z = 0; z < n; ++z) for (int y = 0; y < n; ++y) for (int x = (y+z+parity)&1; x < n; x += 2) {
            const auto i = state.idC(x,y,z);
            float sum = 0.f;
            sum += x > 0 ? p[state.idC(x-1,y,z)] : p[i];
            sum += x+1 < n ? p[state.idC(x+1,y,z)] : p[i];
            sum += y > 0 ? p[state.idC(x,y-1,z)] : p[i];
            sum += y+1 < n ? p[state.idC(x,y+1,z)] : p[i];
            sum += z > 0 ? p[state.idC(x,y,z-1)] : p[i];
            sum += z+1 < n ? p[state.idC(x,y,z+1)] : p[i];
            p[i] = (1.f - omega) * p[i] + omega * (sum - div[i] * h*h / dt) / 6.f;
        }
    MacGridOperators::applyPressureGradient(state, p, dt);
    MacGridOperators::enforceClosedBoundaries(state);
}

void projectionLegacy() {
    for (int n : {2, 3, 6, 12}) for (int iterations : {1, 8, 64}) {
        MacGridState state(n);
        seed(state);
        auto reference = state;
        constexpr float dt = 1.f/60.f, omega = 1.95f;
        SOR solver(omega, 6.f);
        CpuPressureSolver projection(solver);
        ProjectionOptions options;
        options.dt = dt;
        options.linearSolve.maxIterations = iterations;
        options.linearSolve.fixedIterations = true;
        const auto result = projection.project(state, options);
        require(result.linearSolve.iterations == iterations, "Preserved app iteration count");
        legacyProjection(reference, dt, iterations, omega);
        sameVelocities(state, reference, 2e-5f);
    }
}

struct Test { const char* name; void (*run)(); };
const Test tests[] = {
    {"grid_layout", gridLayout}, {"grid_validation", gridValidation},
    {"closed_boundaries", closedBoundaries}, {"sor_manufactured", sorManufactured},
    {"sor_options", sorOptions}, {"sor_validation", sorValidation},
    {"projection_rest", projectionRest}, {"projection_gradient", projectionGradient},
    {"projection_divergence", projectionDivergence}, {"projection_solenoidal", projectionSolenoidal},
    {"projection_injection", projectionInjection}, {"projection_legacy", projectionLegacy}
};
}

int main(int argc, char** argv) {
    int passed = 0, failed = 0;
    for (const auto& test : tests) {
        if (argc > 1 && std::string_view(argv[1]) != test.name) continue;
        try {
            test.run();
            std::cout << "PASS " << test.name << '\n';
            ++passed;
        } catch (const std::exception& error) {
            std::cerr << "FAIL " << test.name << ": " << error.what() << '\n';
            ++failed;
        }
    }
    if (passed + failed == 0) {
        std::cerr << "Unknown test name\n";
        return 1;
    }
    return failed == 0 ? 0 : 1;
}
