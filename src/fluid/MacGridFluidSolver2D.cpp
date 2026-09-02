#include "fluid/MacGridFluidSolver2D.hpp"
#include "fluid/MacGridOperators2D.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

MacGridFluidSolver2D::MacGridFluidSolver2D(int resolution, float boxSize)
    : MacGridFluidSolver2D(MacGridState2D(resolution, boxSize)) {}

MacGridFluidSolver2D::MacGridFluidSolver2D(MacGridState2D state) : state_(std::move(state)),
    uNext_(state_.velocityX().size()), vNext_(state_.velocityY().size()), densityNext_(state_.cellCount()) {}

void MacGridFluidSolver2D::advance(float dt) {
    const auto& p = parameters_;
    if (!std::isfinite(dt) || dt <= 0.f || !std::isfinite(p.sourceStrength) || p.sourceStrength < 0.f ||
        !std::isfinite(p.emitterRadius) || p.emitterRadius <= 0.f || p.emitterRadius > 1.f ||
        !std::isfinite(p.emitterCenterX) || p.emitterCenterX < 0.f || p.emitterCenterX > 1.f ||
        !std::isfinite(p.emitterCenterY) || p.emitterCenterY < 0.f || p.emitterCenterY > 1.f ||
        !std::isfinite(p.buoyancy) || !std::isfinite(p.smokeDecay) || p.smokeDecay < 0.f ||
        !std::isfinite(p.sorOmega) || p.sorOmega <= 0.f || p.sorOmega >= 2.f || p.projectIterations < 0)
        throw std::invalid_argument("Invalid 2D simulation time step or parameters");
    const int n = state_.resolution();
    auto density = state_.density();
    for (int y = 0; y < n; ++y) for (int x = 0; x < n; ++x) {
        const float dx = ((x + .5f) / n - p.emitterCenterX) / p.emitterRadius;
        const float dy = ((y + .5f) / n - p.emitterCenterY) / p.emitterRadius;
        const float weight = std::max(0.f, 1.f - dx * dx - dy * dy);
        density[state_.idC(x, y)] += dt * p.sourceStrength * weight;
    }
    MacGridOperators2D::enforceClosedBoundaries(state_);
    // Both components use the same old velocity; commit only after both passes.
    advection_->advect({state_.domain(), state_.velocityX(), state_.velocityY(), state_.velocityX(), uNext_, dt, MacField2D::VelocityX});
    advection_->advect({state_.domain(), state_.velocityX(), state_.velocityY(), state_.velocityY(), vNext_, dt, MacField2D::VelocityY});
    std::copy(uNext_.begin(), uNext_.end(), state_.velocityX().begin());
    std::copy(vNext_.begin(), vNext_.end(), state_.velocityY().begin());
    for (int y = 1; y < n; ++y) for (int x = 0; x < n; ++x) {
        const float smoke = .5f * (density[state_.idC(x, y - 1)] + density[state_.idC(x, y)]);
        state_.velocityY()[state_.idY(x, y)] += dt * p.buoyancy * smoke;
    }
    ProjectionOptions options;
    options.dt = dt;
    options.relaxation = p.sorOmega;
    options.linearSolve.maxIterations = p.projectIterations;
    lastProjection_ = pressure_->project(state_, options);
    advection_->advect({state_.domain(), state_.velocityX(), state_.velocityY(), density, densityNext_, dt, MacField2D::Density});
    const float decay = std::exp(-dt * p.smokeDecay);
    for (std::size_t i = 0; i < density.size(); ++i) density[i] = decay * densityNext_[i];
}

void MacGridFluidSolver2D::resetState() {
    state_.reset();
    pressure_->reset();
    for (auto* scratch : {&uNext_, &vNext_, &densityNext_}) std::fill(scratch->begin(), scratch->end(), 0.f);
    lastProjection_ = {};
}
void MacGridFluidSolver2D::setBoxSize(float size) {
    const float previous = state_.boxSize();
    state_.setBoxSize(size);
    if (size != previous) resetState();
}
FluidSolverControls MacGridFluidSolver2D::controls() const noexcept {
    return {
        parameters_.sourceStrength,
        parameters_.emitterRadius,
        parameters_.buoyancy,
        parameters_.smokeDecay,
        parameters_.projectIterations,
        parameters_.sorOmega
    };
}
void MacGridFluidSolver2D::setControls(const FluidSolverControls& controls) {
    parameters_.sourceStrength = controls.sourceStrength;
    parameters_.emitterRadius = controls.emitterRadius;
    parameters_.buoyancy = controls.buoyancy;
    parameters_.smokeDecay = controls.smokeDecay;
    parameters_.projectIterations = controls.pressureIterations;
    parameters_.sorOmega = controls.pressureRelaxation;
}
PressureSolveDiagnostics MacGridFluidSolver2D::pressureDiagnostics() const noexcept {
    return {
        lastProjection_.residualAvailable,
        lastProjection_.linearSolve.converged,
        lastProjection_.linearSolve.iterations,
        lastProjection_.linearSolve.finalResidual
    };
}
void MacGridFluidSolver2D::setPressureSolver(IPressureSolver<MacGridState2D>& solver) noexcept {
    pressure_ = &solver;
    lastProjection_ = {};
}
FluidDiagnostics MacGridFluidSolver2D::diagnostics() const {
    FluidDiagnostics result;
    double sum = 0.0;
    // Face quadrature for closed walls, per unit depth and density.
    for (auto field : {state_.velocityX(), state_.velocityY()}) for (float value : field) {
        sum += double(value) * value;
        result.maxVelocity = std::max(result.maxVelocity, std::abs(value));
    }
    result.kineticEnergy = float(.5 * state_.cellSize() * state_.cellSize() * sum);
    result.maxDivergence = MacGridOperators2D::maxDivergence(state_);
    return result;
}
