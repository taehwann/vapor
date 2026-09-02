#include "fluid/MacGridFluidSolver.hpp"
#include "fluid/MacGridOperators.hpp"
#include "math/Vec3.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <utility>

static float frand() {
    static uint32_t state = 0x9E3779B9u;
    state = state * 1664525u + 1013904223u;
    return (state >> 8) * (1.0f / 16777216.0f);
}

MacGridFluidSolver::MacGridFluidSolver(int res)
    : MacGridFluidSolver(MacGridState(res)) {}

MacGridFluidSolver::MacGridFluidSolver(MacGridState state)
    : state_(std::move(state)) {
    parameters_.projectIterations = state_.resolution() * 2;
    allocateBuffers();
}

void MacGridFluidSolver::allocateBuffers() {
    const auto NX = state_.velocityX().size();
    const auto NY = state_.velocityY().size();
    const auto NZ = state_.velocityZ().size();

    advectionWorkspace_.resize(state_);

    vx0.assign(NX, 0.f);
    vy0.assign(NY, 0.f);
    vz0.assign(NZ, 0.f);

    vxTilde.assign(NX, 0.f);
    vyTilde.assign(NY, 0.f);
    vzTilde.assign(NZ, 0.f);

    vxHat.assign(NX, 0.f);
    vyHat.assign(NY, 0.f);
    vzHat.assign(NZ, 0.f);

    vxDiv0.assign(NX, 0.f);
    vyDiv0.assign(NY, 0.f);
    vzDiv0.assign(NZ, 0.f);
}

int MacGridFluidSolver::idC(int x, int y, int z) const { return static_cast<int>(state_.idC(x, y, z)); }

int MacGridFluidSolver::idX(int x, int y, int z) const { return static_cast<int>(state_.idX(x, y, z)); }

int MacGridFluidSolver::idY(int x, int y, int z) const { return static_cast<int>(state_.idY(x, y, z)); }

int MacGridFluidSolver::idZ(int x, int y, int z) const { return static_cast<int>(state_.idZ(x, y, z)); }

float MacGridFluidSolver::h() const { return state_.cellSize(); }

void MacGridFluidSolver::advectScalar(std::span<float> field, float dt) {
    densityAdvection().advect({{state_, velocityView(state_), dt, parameters_.cflMc}, field, field, advectionWorkspace_, std::nullopt});
}

void MacGridFluidSolver::advectVx(float dt, std::span<const float> vxSrc,
    std::span<const float> vxVel, std::span<const float> vyVel, std::span<const float> vzVel,
    std::span<float> out) {
    velocityAdvection().advect({{state_, {vxVel, vyVel, vzVel}, dt, parameters_.cflMc},
                                vxSrc, out, advectionWorkspace_, FaceAxis::X});
}

void MacGridFluidSolver::advectVy(float dt, std::span<const float> vySrc,
    std::span<const float> vxVel, std::span<const float> vyVel, std::span<const float> vzVel,
    std::span<float> out) {
    velocityAdvection().advect({{state_, {vxVel, vyVel, vzVel}, dt, parameters_.cflMc},
                                vySrc, out, advectionWorkspace_, FaceAxis::Y});
}

void MacGridFluidSolver::advectVz(float dt, std::span<const float> vzSrc,
    std::span<const float> vxVel, std::span<const float> vyVel, std::span<const float> vzVel,
    std::span<float> out) {
    velocityAdvection().advect({{state_, {vxVel, vyVel, vzVel}, dt, parameters_.cflMc},
                                vzSrc, out, advectionWorkspace_, FaceAxis::Z});
}

void MacGridFluidSolver::enforceVelocityBoundaries() {
    MacGridOperators::enforceClosedBoundaries(state_);
}

void MacGridFluidSolver::project(float dt, int iterations) {
    ProjectionOptions options;
    options.dt = dt;
    options.relaxation = parameters_.sorOmega;
    options.linearSolve.maxIterations = iterations;
    options.linearSolve.fixedIterations = true;
    lastProjection_ = pressureSolver_->project(state_, options);
}

void MacGridFluidSolver::emit(float strength, float dt) {
    if (strength <= 0.f) return;
    const float dx = h();
    const float dtScale = dt * 60.0f;
    const float ecx = parameters_.emitterCenterX * state_.boxSize(), ecy = parameters_.emitterCenterY * state_.boxSize(), ecz = parameters_.emitterCenterZ * state_.boxSize();
    const float erad = parameters_.emitterRadius * state_.boxSize();
    const int emitR = int(erad / dx + 1);
    const int cx = int(ecx / dx), cy = int(ecy / dx), cz = int(ecz / dx);
    for (int z = std::max(0, cz - emitR); z < std::min(state_.resolution(), cz + emitR); ++z)
        for (int y = std::max(0, cy - emitR); y < std::min(state_.resolution(), cy + emitR); ++y)
            for (int x = std::max(0, cx - emitR); x < std::min(state_.resolution(), cx + emitR); ++x) {
                Vec3 p{ (x + .5f) * dx, (y + .5f) * dx, (z + .5f) * dx };
                float qx = (p.x - ecx) / erad, qy = (p.y - ecy) / erad, qz = (p.z - ecz) / erad;
                if (qx * qx + qy * qy + qz * qz > 1.0f) continue;
                const int i = idC(x, y, z);
                state_.density()[i] = 1.0f;
                const float kick = parameters_.emitterKickBase + parameters_.emitterKickScale * strength;
                state_.velocityZ()[idZ(x, y, z)] = std::max(state_.velocityZ()[idZ(x, y, z)], kick);
                state_.velocityZ()[idZ(x, y, z + 1)] = std::max(state_.velocityZ()[idZ(x, y, z + 1)], kick);
                const float j = (frand() - 0.5f) * parameters_.stirStrength * dtScale;
                state_.velocityX()[idX(x, y, z)] += j; state_.velocityX()[idX(x + 1, y, z)] += j;
                state_.velocityY()[idY(x, y, z)] += j; state_.velocityY()[idY(x, y + 1, z)] += j;
            }
}

void MacGridFluidSolver::applyBuoyancy(float dtStep, float buoyancyVal) {
#pragma omp parallel for
    for (int z = 0; z < state_.resolution(); ++z) {
        for (int y = 0; y < state_.resolution(); ++y) {
            for (int x = 0; x < state_.resolution(); ++x) {
                const float f = buoyancyVal * state_.density()[idC(x, y, z)] * dtStep;
                // Skip bottom and top solid wall faces.
                if (z > 0)     state_.velocityZ()[idZ(x, y, z)] += parameters_.buoyancySplit * f;
                if (z < state_.resolution() - 1) state_.velocityZ()[idZ(x, y, z + 1)] += parameters_.buoyancySplit * f;
            }
        }
    }
}

void MacGridFluidSolver::applyDissipation(float dt) {
    const float sFactor = std::clamp(1.0f - parameters_.smokeDecay * dt, 0.0f, 1.0f);

#pragma omp parallel for
    for (int i = 0; i < (int)state_.density().size(); ++i) {
        state_.density()[i] *= sFactor;
    }
}

void MacGridFluidSolver::step(float dt, float buoyancy, float sourceStrength, int projectIterations) {
    const float halfDt = 0.5f * dt;
    const int halfIters = std::max(1, projectIterations / 2);

    emit(sourceStrength, dt);
    applyBuoyancy(dt, buoyancy);
    enforceVelocityBoundaries();

    vx0.assign(state_.velocityX().begin(), state_.velocityX().end()); vy0.assign(state_.velocityY().begin(), state_.velocityY().end()); vz0.assign(state_.velocityZ().begin(), state_.velocityZ().end());

    advectVx(halfDt, vx0, vx0, vy0, vz0, state_.velocityX());
    vxTilde.assign(state_.velocityX().begin(), state_.velocityX().end());
    advectVy(halfDt, vy0, vx0, vy0, vz0, state_.velocityY());
    vyTilde.assign(state_.velocityY().begin(), state_.velocityY().end());
    advectVz(halfDt, vz0, vx0, vy0, vz0, state_.velocityZ());
    vzTilde.assign(state_.velocityZ().begin(), state_.velocityZ().end());

    project(halfDt, halfIters);

    vxDiv0.assign(state_.velocityX().begin(), state_.velocityX().end()); vyDiv0.assign(state_.velocityY().begin(), state_.velocityY().end()); vzDiv0.assign(state_.velocityZ().begin(), state_.velocityZ().end());

#pragma omp parallel for
    for (int i = 0; i < (int)state_.velocityX().size(); ++i) vxHat[i] = 2.f * state_.velocityX()[i] - vxTilde[i];
#pragma omp parallel for
    for (int i = 0; i < (int)state_.velocityY().size(); ++i) vyHat[i] = 2.f * state_.velocityY()[i] - vyTilde[i];
#pragma omp parallel for
    for (int i = 0; i < (int)state_.velocityZ().size(); ++i) vzHat[i] = 2.f * state_.velocityZ()[i] - vzTilde[i];

    advectVx(halfDt, vxHat, vxDiv0, vyDiv0, vzDiv0, state_.velocityX());
    advectVy(halfDt, vyHat, vxDiv0, vyDiv0, vzDiv0, state_.velocityY());
    advectVz(halfDt, vzHat, vxDiv0, vyDiv0, vzDiv0, state_.velocityZ());

    project(halfDt, halfIters);

    advectScalar(state_.density(), dt);
    applyDissipation(dt);
}

void MacGridFluidSolver::stepNoReflect(float dt, float buoyancy, float sourceStrength, int projectIterations) {
    emit(sourceStrength, dt);
    applyBuoyancy(dt, buoyancy);
    enforceVelocityBoundaries();

    vx0.assign(state_.velocityX().begin(), state_.velocityX().end()); vy0.assign(state_.velocityY().begin(), state_.velocityY().end()); vz0.assign(state_.velocityZ().begin(), state_.velocityZ().end());

    advectVx(dt, vx0, vx0, vy0, vz0, state_.velocityX());
    advectVy(dt, vy0, vx0, vy0, vz0, state_.velocityY());
    advectVz(dt, vz0, vx0, vy0, vz0, state_.velocityZ());

    project(dt, projectIterations);

    advectScalar(state_.density(), dt);
    applyDissipation(dt);
}

float MacGridFluidSolver::kineticEnergy() const {
    const float dx = h();
    double ke = 0.0;
    for (int z = 0; z < state_.resolution(); ++z) for (int y = 0; y < state_.resolution(); ++y) for (int x = 0; x <= state_.resolution(); ++x)
        ke += (double)state_.velocityX()[idX(x, y, z)] * state_.velocityX()[idX(x, y, z)];
    for (int z = 0; z < state_.resolution(); ++z) for (int y = 0; y <= state_.resolution(); ++y) for (int x = 0; x < state_.resolution(); ++x)
        ke += (double)state_.velocityY()[idY(x, y, z)] * state_.velocityY()[idY(x, y, z)];
    for (int z = 0; z <= state_.resolution(); ++z) for (int y = 0; y < state_.resolution(); ++y) for (int x = 0; x < state_.resolution(); ++x)
        ke += (double)state_.velocityZ()[idZ(x, y, z)] * state_.velocityZ()[idZ(x, y, z)];
    return float(0.5 * ke * dx * dx * dx);
}

float MacGridFluidSolver::maxVelocity() const {
    float m = 0.f;
    for (float f : state_.velocityX()) m = std::max(m, std::abs(f));
    for (float f : state_.velocityY()) m = std::max(m, std::abs(f));
    for (float f : state_.velocityZ()) m = std::max(m, std::abs(f));
    return m;
}

float MacGridFluidSolver::computeDivergenceNorm() const {
    return MacGridOperators::maxDivergence(state_);
}

void MacGridFluidSolver::resetState() {
    state_.reset();
    defaultPressure_.reset();
    if (pressureSolver_ != &defaultPressure_) pressureSolver_->reset();
    lastProjection_ = {};
    advectionWorkspace_.reset();
    std::fill(vx0.begin(), vx0.end(), 0.f);
    std::fill(vy0.begin(), vy0.end(), 0.f);
    std::fill(vz0.begin(), vz0.end(), 0.f);
    std::fill(vxTilde.begin(), vxTilde.end(), 0.f);
    std::fill(vyTilde.begin(), vyTilde.end(), 0.f);
    std::fill(vzTilde.begin(), vzTilde.end(), 0.f);
    std::fill(vxHat.begin(), vxHat.end(), 0.f);
    std::fill(vyHat.begin(), vyHat.end(), 0.f);
    std::fill(vzHat.begin(), vzHat.end(), 0.f);
    std::fill(vxDiv0.begin(), vxDiv0.end(), 0.f);
    std::fill(vyDiv0.begin(), vyDiv0.end(), 0.f);
    std::fill(vzDiv0.begin(), vzDiv0.end(), 0.f);
}

int MacGridFluidSolver::resolution() const { return state_.resolution(); }

float MacGridFluidSolver::boxSize() const noexcept { return state_.boxSize(); }

const ProjectionResult& MacGridFluidSolver::lastProjection() const { return lastProjection_; }

void MacGridFluidSolver::setBoxSize(float size) {
    if (size != state_.boxSize()) {
        state_.setBoxSize(size);
        resetState();
    }
}

const IAdvection<MacAdvectionOperation>& MacGridFluidSolver::velocityAdvection() const {
    if (velocityAdvection_) return *velocityAdvection_;
    return parameters_.macCormackVel ? static_cast<const IAdvection<MacAdvectionOperation>&>(macCormack)
                                    : static_cast<const IAdvection<MacAdvectionOperation>&>(semiLagrangian);
}

const IAdvection<MacAdvectionOperation>& MacGridFluidSolver::densityAdvection() const {
    if (densityAdvection_) return *densityAdvection_;
    return parameters_.macCormackSmoke ? static_cast<const IAdvection<MacAdvectionOperation>&>(macCormack)
                                      : static_cast<const IAdvection<MacAdvectionOperation>&>(semiLagrangian);
}

void MacGridFluidSolver::setPressureSolver(IPressureSolver<MacGridState>& solver) noexcept {
    pressureSolver_ = &solver;
    lastProjection_ = {};
}

void MacGridFluidSolver::useDefaultPressureSolver() noexcept {
    setPressureSolver(defaultPressure_);
}

void MacGridFluidSolver::setVelocityAdvection(const IAdvection<MacAdvectionOperation>& method) noexcept {
    velocityAdvection_ = &method;
}

void MacGridFluidSolver::setDensityAdvection(const IAdvection<MacAdvectionOperation>& method) noexcept {
    densityAdvection_ = &method;
}

void MacGridFluidSolver::useDefaultAdvection() noexcept {
    velocityAdvection_ = densityAdvection_ = nullptr;
}

void MacGridFluidSolver::advance(float dt) {
    if (parameters_.reflection)
        step(dt, parameters_.buoyancy, parameters_.sourceStrength, parameters_.projectIterations);
    else
        stepNoReflect(dt, parameters_.buoyancy, parameters_.sourceStrength, parameters_.projectIterations);
}

FluidDiagnostics MacGridFluidSolver::diagnostics() const {
    return {kineticEnergy(), maxVelocity(), computeDivergenceNorm()};
}

FluidSolverControls MacGridFluidSolver::controls() const noexcept {
    return {
        parameters_.sourceStrength,
        parameters_.emitterRadius,
        parameters_.buoyancy,
        parameters_.smokeDecay,
        parameters_.projectIterations,
        parameters_.sorOmega,
        parameters_.stirStrength,
        parameters_.cflMc
    };
}

void MacGridFluidSolver::setControls(const FluidSolverControls& controls) {
    parameters_.sourceStrength = controls.sourceStrength;
    parameters_.emitterRadius = controls.emitterRadius;
    parameters_.buoyancy = controls.buoyancy;
    parameters_.smokeDecay = controls.smokeDecay;
    parameters_.projectIterations = controls.pressureIterations;
    parameters_.sorOmega = controls.pressureRelaxation;
    parameters_.stirStrength = controls.stirStrength;
    parameters_.cflMc = controls.advectionCfl;
}

PressureSolveDiagnostics MacGridFluidSolver::pressureDiagnostics() const noexcept {
    return {
        lastProjection_.residualAvailable,
        lastProjection_.linearSolve.converged,
        lastProjection_.linearSolve.iterations,
        lastProjection_.linearSolve.finalResidual
    };
}
