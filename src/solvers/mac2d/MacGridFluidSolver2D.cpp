#include "solvers/mac2d/MacGridFluidSolver2D.hpp"
#include "solvers/mac2d/MacGridOperators2D.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

MacGridFluidSolver2D::MacGridFluidSolver2D() : MacGridFluidSolver2D(128, 1.f) {
    parameters_.emitterRadius = .04f;
    parameters_.sourceStrength = 1.f;
    parameters_.buoyancy = .8f;
    parameters_.smokeDecay = .3007525f;
    parameters_.projectIterations = 120;
}
MacGridFluidSolver2D::MacGridFluidSolver2D(int resolution, float boxSize)
    : MacGridFluidSolver2D(MacGridState2D(resolution, boxSize)) {}

MacGridFluidSolver2D::MacGridFluidSolver2D(MacGridState2D state) : state_(std::move(state)),
    uNext_(state_.velocityX().size()), vNext_(state_.velocityY().size()), densityNext_(state_.cellCount()),
    uHat_(state_.velocityX().size()), vHat_(state_.velocityY().size()) {}

void MacGridFluidSolver2D::advance(float dt) {
    const auto& p = parameters_;
    if (!std::isfinite(p.emitterKickMultiplier) || p.emitterKickMultiplier < 0.f)
        throw std::invalid_argument("Invalid emitter kick multiplier");
    if (!std::isfinite(dt) || dt <= 0.f || !std::isfinite(p.sourceStrength) || p.sourceStrength < 0.f ||
        !std::isfinite(p.emitterRadius) || p.emitterRadius <= 0.f || p.emitterRadius > 1.f ||
        !std::isfinite(p.emitterCenterX) || p.emitterCenterX < 0.f || p.emitterCenterX > 1.f ||
        !std::isfinite(p.emitterCenterY) || p.emitterCenterY < 0.f || p.emitterCenterY > 1.f ||
        !std::isfinite(p.buoyancy) || !std::isfinite(p.smokeDecay) || p.smokeDecay < 0.f ||
        !std::isfinite(p.sorOmega) || p.sorOmega <= 0.f || p.sorOmega >= 2.f || p.projectIterations < 0)
        throw std::invalid_argument("Invalid 2D simulation time step or parameters");
    MacGridOperators2D::enforceClosedBoundaries(state_);
    midpointProjection_ = {};
    if (p.reflection) {
        advanceReflection(dt);
        return;
    }
    advectVelocity(state_.velocityX(), state_.velocityY(), dt);
    addBuoyancy(dt);
    emit();
    lastProjection_ = project(dt);
    advectDensity(dt);
}

void MacGridFluidSolver2D::emit() {
    const auto& p = parameters_;
    if (p.sourceStrength <= 0.f) return;
    const int n = state_.resolution();
    auto density = state_.density();
    // Match the reflection plume's replenished Gaussian source and inlet kick.
    // Positions/radius are fractions of the box; velocity is in world units.
    const float r2 = p.emitterRadius * p.emitterRadius;
    for (int y = 0; y < n; ++y) for (int x = 0; x < n; ++x) {
        const float dx = (x + .5f) / n - p.emitterCenterX;
        const float dy = (y + .5f) / n - p.emitterCenterY;
        const float q = (dx * dx + dy * dy) / r2;
        if (q > 1.f) continue;
        const float weight = std::exp(-3.5f * q);
        const auto i = state_.idC(x, y);
        density[i] = std::max(density[i], std::min(1.f, p.sourceStrength * 3.f * weight));
        auto& v = state_.velocityY()[state_.idY(x, y)];
        if (p.emitterKickMultiplier > 0.f) v = std::max(v, p.emitterKickMultiplier * (.6f + p.sourceStrength) * weight);
    }
}

void MacGridFluidSolver2D::advectVelocity(std::span<const float> u, std::span<const float> v, float dt) {
    // Both components use the same old velocity; commit only after both passes.
    advection_.advect({state_.domain(), state_.velocityX(), state_.velocityY(), u, uNext_, dt, MacField2D::VelocityX});
    advection_.advect({state_.domain(), state_.velocityX(), state_.velocityY(), v, vNext_, dt, MacField2D::VelocityY});
    std::copy(uNext_.begin(), uNext_.end(), state_.velocityX().begin());
    std::copy(vNext_.begin(), vNext_.end(), state_.velocityY().begin());
}

void MacGridFluidSolver2D::addBuoyancy(float dt) {
    const auto& p = parameters_;
    const int n = state_.resolution();
    const auto density = state_.density();
    for (int y = 1; y < n; ++y) for (int x = 0; x < n; ++x) {
        // Apply the nonlinear force per cell before averaging onto the face,
        // as in ReflectionMacFluidSolver2D's two half-face contributions.
        const float smoke = .5f * (std::pow(density[state_.idC(x, y - 1)], .25f) +
                                  std::pow(density[state_.idC(x, y)], .25f));
        state_.velocityY()[state_.idY(x, y)] += dt * p.buoyancy * smoke;
    }
}

ProjectionResult MacGridFluidSolver2D::project(float dt) {
    const auto& p = parameters_;
    ProjectionOptions options;
    options.dt = dt;
    options.relaxation = p.sorOmega;
    options.linearSolve.maxIterations = p.projectIterations;
    return pressure_.project(state_, options);
}

void MacGridFluidSolver2D::advectDensity(float dt) {
    auto density = state_.density();
    advection_.advect({state_.domain(), state_.velocityX(), state_.velocityY(), density, densityNext_, dt, MacField2D::Density});
    const float decay = std::exp(-dt * parameters_.smokeDecay);
    for (std::size_t i = 0; i < density.size(); ++i) density[i] = decay * densityNext_[i];
}

void MacGridFluidSolver2D::advanceReflection(float dt) {
    const float halfDt = .5f * dt;
    advectDensity(halfDt);
    advectVelocity(state_.velocityX(), state_.velocityY(), halfDt);
    // uNext/vNext retain the advected, pre-force velocity (u_tilde).
    // Reflect the projected half-force update to account for a full force step.
    addBuoyancy(halfDt);
    emit();
    midpointProjection_ = project(halfDt);
    for (std::size_t i = 0; i < uHat_.size(); ++i)
        uHat_[i] = 2.f * state_.velocityX()[i] - uNext_[i];
    for (std::size_t i = 0; i < vHat_.size(); ++i)
        vHat_[i] = 2.f * state_.velocityY()[i] - vNext_[i];
    // The solid-wall reflection must not feed an artificial normal velocity
    // into the second advection pass.
    const int n = state_.resolution();
    for (int i = 0; i < n; ++i) {
        uHat_[state_.idX(0,i)] = uHat_[state_.idX(n,i)] = 0.f;
        vHat_[state_.idY(i,0)] = vHat_[state_.idY(i,n)] = 0.f;
    }
    advectDensity(halfDt);
    // Both reflected components trace through the same projected midpoint.
    advectVelocity(uHat_, vHat_, halfDt);
    lastProjection_ = project(halfDt);
}

void MacGridFluidSolver2D::resetState() {
    state_.reset();
    pressure_.reset();
    for (auto* scratch : {&uNext_, &vNext_, &densityNext_, &uHat_, &vHat_}) std::fill(scratch->begin(), scratch->end(), 0.f);
    lastProjection_ = {};
    midpointProjection_ = {};
}
void MacGridFluidSolver2D::setBoxSize(float size) {
    const float previous = state_.boxSize();
    state_.setBoxSize(size);
    if (size != previous) resetState();
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
