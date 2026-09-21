#include "solvers/mac2d/MacGridFluidSolver2D.hpp"
#include "solvers/mac2d/ReflectionMacFluidSolver2D.hpp"
#include "solvers/mac2d/MacGridOperators2D.hpp"
#include "renderer/MacRenderData2D.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <numbers>
#include <numeric>
#include <stdexcept>
#include <string>
#include <string_view>

namespace {
void require(bool condition, const char* message) { if (!condition) throw std::runtime_error(message); }
void near(double actual, double expected, double tolerance, const char* message) {
    if (!std::isfinite(actual) || std::abs(actual - expected) > tolerance)
        throw std::runtime_error(std::string(message) + ": " + std::to_string(actual) + " vs " + std::to_string(expected));
}
template<class Exception = std::invalid_argument, class F> void rejects(F&& f) {
    try { f(); } catch (const Exception&) { return; }
    throw std::runtime_error("Expected invalid input to be rejected");
}
void same(std::span<const float> a, std::span<const float> b, float tolerance = 1e-6f) {
    require(a.size() == b.size(), "Field extent mismatch");
    for (std::size_t i = 0; i < a.size(); ++i) near(a[i], b[i], tolerance, "Field mismatch");
}
void walls(const MacGridState2D& s) {
    for (int i = 0; i < s.resolution(); ++i) {
        near(s.velocityX()[s.idX(0, i)], 0, 0, "Left wall");
        near(s.velocityX()[s.idX(s.resolution(), i)], 0, 0, "Right wall");
        near(s.velocityY()[s.idY(i, 0)], 0, 0, "Bottom wall");
        near(s.velocityY()[s.idY(i, s.resolution())], 0, 0, "Top wall");
    }
}
ProjectionOptions accurate(float dt = .02f) {
    ProjectionOptions o;
    o.dt = dt;
    o.relaxation = 1.6f;
    o.linearSolve.maxIterations = 4000;
    o.linearSolve.absoluteTolerance = 3e-6;
    o.linearSolve.relativeTolerance = 1e-6;
    return o;
}
void grid() {
    for (int n : {1, 2, 7}) {
        MacGridState2D s(n, 3.f);
        const auto& domain = s.domain();
        require(domain.dimension() == Dimension::D2, "2D dimension tag");
        near(domain.bounds().max.z, 0, 0, "Planar bounds");
        near(s.cellSize(), 3.0 / n, 1e-6, "2D cell spacing");
        require(s.cellCount() == std::size_t(n * n), "2D has no hidden Z cells");
        for (int field = 0; field < 3; ++field) {
            const int w = n + (field == 1), h = n + (field == 2);
            auto values = field == 0 ? s.density() : field == 1 ? s.velocityX() : s.velocityY();
            require(values.size() == std::size_t(w * h), "2D staggered extent");
            std::vector<int> visited(w * h);
            for (int y = 0; y < h; ++y) for (int x = 0; x < w; ++x) {
                auto i = field == 0 ? s.idC(x, y) : field == 1 ? s.idX(x, y) : s.idY(x, y);
                require(i < visited.size(), "2D index bounds");
                ++visited[i];
            }
            for (int visits : visited) require(visits == 1, "2D index uniqueness");
            std::fill(values.begin(), values.end(), 1.f);
        }
        s.setBoxSize(2.f);
        for (auto f : {s.density(), s.velocityX(), s.velocityY()}) for (float v : f) near(v, 0, 0, "Resize resets fields");
    }
    rejects([] { MacGridState2D s(0); });
    rejects([] { MacGridState2D s(-1); });
    rejects<std::length_error>([] { MacGridState2D s(std::numeric_limits<int>::max()); });
    rejects([] { MacGridState2D s(4, 0); });
    rejects([] { MacGridState2D s(4, std::numeric_limits<float>::infinity()); });
}
void advection() {
    MacGridState2D s(8, 1.f);
    std::fill(s.velocityX().begin(), s.velocityX().end(), .5f);
    std::fill(s.velocityY().begin(), s.velocityY().end(), -.25f);
    SemiLagrangian2D sl;
    const SemiLagrangian2D& transport = sl;
    constexpr float dt = .125f;
    for (auto kind : {MacField2D::Density, MacField2D::VelocityX, MacField2D::VelocityY}) {
        const int w = 8 + (kind == MacField2D::VelocityX), h = 8 + (kind == MacField2D::VelocityY);
        const float ox = kind == MacField2D::VelocityX ? 0.f : .5f;
        const float oy = kind == MacField2D::VelocityY ? 0.f : .5f;
        std::vector<float> source(w * h), output(w * h);
        for (int y = 0; y < h; ++y) for (int x = 0; x < w; ++x)
            source[x + w * y] = 1.f + .2f * (x + ox) * s.cellSize() + .3f * (y + oy) * s.cellSize();
        transport.advect({s.domain(), s.velocityX(), s.velocityY(), source, output, 0.f, kind});
        same(source, output);
        transport.advect({s.domain(), s.velocityX(), s.velocityY(), source, output, dt, kind});
        for (int y = 1; y < h - 1; ++y) for (int x = 1; x < w - 1; ++x) {
            const float expected = 1.f + .2f * ((x + ox) * s.cellSize() - .5f * dt) +
                                   .3f * ((y + oy) * s.cellSize() + .25f * dt);
            near(output[x + w * y], expected, 3e-7, "Bilinear translation in world units");
        }
        const auto [low, high] = std::minmax_element(source.begin(), source.end());
        for (float value : output) require(value >= *low && value <= *high, "Semi-Lagrangian introduced an extremum");
        std::fill(source.begin(), source.end(), .7f);
        transport.advect({s.domain(), s.velocityX(), s.velocityY(), source, output, 100.f, kind});
        same(source, output);
        rejects([&] { transport.advect({s.domain(), s.velocityX(), s.velocityY(), source, source, dt, kind}); });
        rejects([&] { transport.advect({s.domain(), s.velocityX(), s.velocityY(), source, output, -dt, kind}); });
    }
    std::vector<float> out(s.cellCount());
    rejects([&] { transport.advect({s.domain(), s.velocityX().first(1), s.velocityY(), s.density(), out, dt, MacField2D::Density}); });
    rejects([&] { transport.advect({s.domain(), s.velocityX(), s.velocityY(), s.velocityX(), s.velocityX(), dt, MacField2D::VelocityX}); });
    s.velocityX()[0] = std::numeric_limits<float>::quiet_NaN();
    rejects([&] { transport.advect({s.domain(), s.velocityX(), s.velocityY(), s.density(), out, dt, MacField2D::Density}); });
}
void projectionGradient() {
    CpuPressureSolver2D solver;
    CpuPressureSolver2D& projection = solver;
    for (int n : {1, 8, 15}) for (float box : {1.f, 3.f}) for (float dt : {.01f, .1f}) {
        MacGridState2D s(n, box);
        std::vector<float> exact(n * n);
        for (int y = 0; y < n; ++y) for (int x = 0; x < n; ++x)
            exact[s.idC(x, y)] = std::cos(std::numbers::pi_v<float> * (x + .5f) / n) *
                                 std::cos(2.f * std::numbers::pi_v<float> * (y + .5f) / n);
        const double mean = std::accumulate(exact.begin(), exact.end(), 0.0) / exact.size();
        for (auto& p : exact) p -= float(mean);
        for (int y = 0; y < n; ++y) for (int x = 1; x < n; ++x)
            s.velocityX()[s.idX(x, y)] = dt / s.cellSize() * (exact[s.idC(x, y)] - exact[s.idC(x - 1, y)]);
        for (int y = 1; y < n; ++y) for (int x = 0; x < n; ++x)
            s.velocityY()[s.idY(x, y)] = dt / s.cellSize() * (exact[s.idC(x, y)] - exact[s.idC(x, y - 1)]);
        const auto result = projection.project(s, accurate(dt));
        require(result.residualAvailable && result.linearSolve.converged, "Manufactured pressure did not converge");
        require(solver.pressure().size() == exact.size(), "Pressure extent");
        for (std::size_t i = 0; i < exact.size(); ++i) near(solver.pressure()[i], exact[i], 8e-5, "Manufactured pressure recovery");
        require(result.divergenceAfter < 2e-4f, "Projection failed to remove gradient divergence");
        for (auto velocity : {s.velocityX(), s.velocityY()})
            for (float value : velocity) near(value, 0, 3e-5, "Projection must remove a pure gradient");
        walls(s);
    }
}
void projectionSolenoidal() {
    MacGridState2D s(12, 2.f);
    const int n = s.resolution();
    auto psi = [n](int x, int y) {
        if (x == 0 || y == 0 || x == n || y == n) return 0.f;
        return std::sin(std::numbers::pi_v<float> * x / n) * std::sin(std::numbers::pi_v<float> * y / n);
    };
    for (int y = 0; y < n; ++y) for (int x = 0; x <= n; ++x)
        s.velocityX()[s.idX(x, y)] = (psi(x, y + 1) - psi(x, y)) / s.cellSize();
    for (int y = 0; y <= n; ++y) for (int x = 0; x < n; ++x)
        s.velocityY()[s.idY(x, y)] = -(psi(x + 1, y) - psi(x, y)) / s.cellSize();
    const auto original = s;
    CpuPressureSolver2D pressure;
    pressure.project(s, accurate());
    same(s.velocityX(), original.velocityX(), 3e-6f);
    same(s.velocityY(), original.velocityY(), 3e-6f);
    walls(s);
}
void projectionOptions() {
    MacGridState2D s(8);
    CpuPressureSolver2D solver;
    auto options = accurate();
    require(solver.project(s, options).linearSolve.iterations == 0, "Rest projection should exit early");
    options.linearSolve.fixedIterations = true;
    options.linearSolve.maxIterations = 17;
    require(solver.project(s, options).linearSolve.iterations == 17, "Fixed SOR budget");
    options = accurate();
    s.velocityX()[s.idX(4, 3)] = 1.f;
    const auto original = s;
    options.linearSolve.maxIterations = 0;
    require(!solver.project(s, options).linearSolve.converged, "Zero budget must not claim convergence");
    same(s.velocityX(), original.velocityX());
    for (float dt : {0.f, -1.f, std::numeric_limits<float>::infinity()}) {
        options = accurate(dt);
        rejects([&] { solver.project(s, options); });
    }
    options = accurate(); options.relaxation = 2.f;
    rejects([&] { solver.project(s, options); });
    options = accurate(); options.linearSolve.residualCheckInterval = 0;
    rejects([&] { solver.project(s, options); });
    options = accurate(); options.linearSolve.maxIterations = -1;
    rejects([&] { solver.project(s, options); });
    options = accurate(); options.linearSolve.absoluteTolerance = std::numeric_limits<double>::quiet_NaN();
    rejects([&] { solver.project(s, options); });
    solver.reset();
    for (double p : solver.pressure()) near(p, 0, 0, "Pressure reset");
    options = accurate();
    const auto result = solver.project(s, options);
    require(result.divergenceAfter < result.divergenceBefore * 1e-4f, "Localized divergence reduction");
    s.velocityY()[0] = std::numeric_limits<float>::quiet_NaN();
    rejects([&] { solver.project(s, options); });

    // Default 64x64 settings must resolve a large, smooth pressure mode, not
    // only converge on the smaller numerical fixtures above.
    MacGridState2D large(64);
    MacGridParameters2D defaults;
    options = {};
    options.relaxation = defaults.sorOmega;
    options.linearSolve.maxIterations = defaults.projectIterations;
    for (int y = 1; y < 64; ++y) for (int x = 0; x < 64; ++x)
        large.velocityY()[large.idY(x, y)] = 4.f * options.dt;
    const auto largeResult = solver.project(large, options);
    require(largeResult.linearSolve.converged, "Default SOR budget failed on a 64x64 pressure mode");
    require(largeResult.divergenceAfter < 1e-4f, "Default 64x64 projection divergence");
}
double centerY(const MacGridState2D& s) {
    double mass = 0, moment = 0;
    for (int y = 0; y < s.resolution(); ++y) for (int x = 0; x < s.resolution(); ++x) {
        const float d = s.density()[s.idC(x, y)];
        mass += d; moment += d * (y + .5) * s.cellSize();
    }
    require(mass > 0, "Expected emitted smoke");
    return moment / mass;
}
void simulation() {
    // Isolate the source from transport/projection. It must agree with the
    // reference inlet, including weak/disabled sources and saturated density.
    for (float strength : {0.f, .1f, 1.f, 4.f}) {
        ReflectionMacFluidSolver2D reference(64);
        MacGridFluidSolver2D source(64, 1.f);
        source.parameters().sourceStrength = strength;
        source.parameters().buoyancy = source.parameters().smokeDecay = 0.f;
        source.parameters().projectIterations = 0;
        reference.emit(strength, 1e-10f);
        reference.applyBoundary();
        source.advance(1e-10f);
        same(source.state().density(), reference.smoke);
        same(source.state().velocityY(), reference.vy);
        source.advance(1e-10f);
        same(source.state().density(), reference.smoke);
        same(source.state().velocityY(), reference.vy);
    }
    MacGridFluidSolver2D solver(24);
    auto& fluid = solver;
    solver.parameters().sourceStrength = 0.f;
    for (int i = 0; i < 5; ++i) fluid.advance(1.f / 60.f);
    near(fluid.diagnostics().kineticEnergy, 0, 0, "Rest energy");
    solver.parameters().sourceStrength = 4.f;
    for (int i = 0; i < 120; ++i) {
        fluid.advance(1.f / 60.f);
        require(fluid.diagnostics().maxDivergence < .005f, "Smoke step divergence too high");
        walls(solver.state());
        for (float d : solver.state().density()) require(std::isfinite(d) && d >= 0.f && d <= 8.01f, "Smoke density bounds");
    }
    require(fluid.diagnostics().kineticEnergy > .001f, "Buoyancy must produce motion");
    const double before = centerY(solver.state());
    solver.parameters().sourceStrength = 0.f;
    for (int i = 0; i < 30; ++i) fluid.advance(1.f / 60.f);
    require(centerY(solver.state()) > before + .02, "Smoke should rise after emission stops");
    std::cout << "plume energy=" << fluid.diagnostics().kineticEnergy << " divergence=" << fluid.diagnostics().maxDivergence
              << " centerY=" << centerY(solver.state()) << '\n';
    const auto original = solver.state();
    rejects([&] { fluid.advance(0.f); });
    rejects([&] { fluid.advance(-.1f); });
    same(solver.state().density(), original.density());
    fluid.resetState();
    for (float d : solver.state().density()) near(d, 0, 0, "Reset smoke");
    near(fluid.diagnostics().kineticEnergy, 0, 0, "Reset energy");
}
void reflection() {
    // A projected force has an independently known first-order response dt*P(f).
    // This catches missing/doubled buoyancy in the reflected half-step.
    MacGridState2D initial(20, 4.f);
    for (int y=0; y<20; ++y) for (int x=0; x<20; ++x)
        initial.density()[initial.idC(x,y)] = float(std::exp(-.08*((x-8)*(x-8)+(y-7)*(y-7))));
    auto expected = initial;
    constexpr float dt = 1e-4f;
    for (int y=1; y<20; ++y) for (int x=0; x<20; ++x)
        expected.velocityY()[expected.idY(x,y)] = dt * 3.f * .5f *
            (std::pow(initial.density()[initial.idC(x,y-1)], .25f) +
             std::pow(initial.density()[initial.idC(x,y)], .25f));
    CpuPressureSolver2D pressure;
    pressure.project(expected, accurate(dt));
    MacGridFluidSolver2D simpleForced(initial);
    simpleForced.parameters().sourceStrength = simpleForced.parameters().smokeDecay = 0.f;
    simpleForced.advance(dt);
    same(simpleForced.state().velocityX(), expected.velocityX(), 2e-7f);
    same(simpleForced.state().velocityY(), expected.velocityY(), 2e-7f);
    MacGridFluidSolver2D forced(initial);
    forced.parameters().reflection = true;
    forced.parameters().sourceStrength = forced.parameters().smokeDecay = 0.f;
    forced.advance(dt);
    same(forced.state().velocityX(), expected.velocityX(), 2e-7f);
    same(forced.state().velocityY(), expected.velocityY(), 2e-7f);

    // A nonuniform solenoidal vortex exercises transport of distinct reflected
    // components through one immutable midpoint velocity.
    const auto psi = [](int x, int y) {
        return .8 * std::sin(std::numbers::pi * x / 20.) * std::sin(2. * std::numbers::pi * y / 20.);
    };
    for (int y=0; y<20; ++y) for (int x=0; x<=20; ++x)
        initial.velocityX()[initial.idX(x,y)] = float((psi(x,y+1)-psi(x,y)) / initial.cellSize());
    for (int y=0; y<=20; ++y) for (int x=0; x<20; ++x)
        initial.velocityY()[initial.idY(x,y)] = float(-(psi(x+1,y)-psi(x,y)) / initial.cellSize());
    MacGridFluidSolver2D reflected(initial), simple(initial), replay(initial);
    for (auto* solver : {&reflected, &simple, &replay}) {
        solver->parameters().sourceStrength = 0.f;
        solver->parameters().buoyancy = 0.f;
        solver->parameters().reflection = solver != &simple;
    }
    for (int frame=0; frame<40; ++frame) {
        reflected.advance(.02f);
        simple.advance(.02f);
        replay.advance(.02f);
        require(reflected.midpointProjection().linearSolve.converged, "Midpoint pressure failed");
        require(reflected.lastProjection().linearSolve.converged, "Final pressure failed");
        require(reflected.diagnostics().maxDivergence < 1e-4f, "Reflection divergence");
        walls(reflected.state());
    }
    same(reflected.state().velocityX(), replay.state().velocityX(), 0.f);
    same(reflected.state().velocityY(), replay.state().velocityY(), 0.f);
    double difference=0;
    for (std::size_t i=0;i<reflected.state().velocityX().size();++i)
        difference=std::max(difference,double(std::abs(reflected.state().velocityX()[i]-simple.state().velocityX()[i])));
    std::cout << "reflection energy=" << reflected.diagnostics().kineticEnergy
              << " simple energy=" << simple.diagnostics().kineticEnergy << " difference=" << difference << '\n';
    require(difference > 1e-4, "Reflection must differ from simple stepping");
    reflected.resetState();
    require(!reflected.midpointProjection().residualAvailable, "Reset midpoint diagnostics");
    reflected.advance(.02f);
    near(reflected.diagnostics().kineticEnergy, 0., 0., "Reflection reset/rest");

    // Default-resolution plume with forcing and both projection stages active.
    MacGridFluidSolver2D plume(64);
    plume.parameters().reflection = true;
    for (int frame=0;frame<240;++frame) {
        plume.advance(1.f/60.f);
        require(plume.midpointProjection().linearSolve.converged, "Plume midpoint convergence");
        require(plume.lastProjection().linearSolve.converged, "Plume final convergence");
        require(plume.diagnostics().maxDivergence < .001f, "Reflected plume divergence");
        for (float d : plume.state().density()) require(std::isfinite(d) && d>=0.f, "Reflected plume density");
    }
    std::cout << "reflection plume energy=" << plume.diagnostics().kineticEnergy
              << " divergence=" << plume.diagnostics().maxDivergence << '\n';
}

void presentation() {
    MacGridFluidSolver2D solver(16);
    for (int i = 0; i < 10; ++i) solver.advance(.02f);
    const auto image = imageData(solver.state());
    require(image.density.data() == solver.state().density().data() && image.width == 16 && image.height == 16,
            "2D borrowed image");
    solver.resetState();
    for (float d : solver.state().density()) near(d, 0, 0, "Reset density");
    solver.setBoxSize(2.f);
    near(imageData(solver.state()).worldHeight, 2, 0, "View follows domain change");
    solver.advance(.02f);
    require(solver.lastProjection().residualAvailable, "CPU projection diagnostics");
}

}
int main(int argc, char** argv) {
    try {
        const std::string_view name = argc > 1 ? argv[1] : "";
        if (name == "mac2d_grid") grid();
        else if (name == "mac2d_advection") advection();
        else if (name == "mac2d_projection_gradient") projectionGradient();
        else if (name == "mac2d_projection_solenoidal") projectionSolenoidal();
        else if (name == "mac2d_projection_options") projectionOptions();
        else if (name == "mac2d_simulation") simulation();
        else if (name == "mac2d_reflection") reflection();
        else if (name == "mac2d_presentation") presentation();
        else throw std::invalid_argument("Unknown 2D test case");
        std::cout << "PASS " << name << '\n';
    } catch (const std::exception& error) {
        std::cerr << "FAIL 2D: " << error.what() << '\n'; return 1;
    }
}
