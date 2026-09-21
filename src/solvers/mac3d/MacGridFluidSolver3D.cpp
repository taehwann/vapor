#include "solvers/mac3d/MacGridFluidSolver3D.hpp"
#include "solvers/mac3d/MacGridOperators3D.hpp"
#include "common/Vec3.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <utility>
#include <stdexcept>

MacGridParameters3D::MacGridParameters3D(MacAlgorithm algorithm) {
    switch (algorithm) {
    case MacAlgorithm::SimpleSemiLagrangian: reflection = false; macCormackSmoke = false; break;
    case MacAlgorithm::SimpleMacCormack: reflection = false; macCormackSmoke = true; break;
    case MacAlgorithm::ReflectionSemiLagrangian: reflection = true; macCormackSmoke = false; break;
    case MacAlgorithm::ReflectionMacCormack: reflection = true; macCormackSmoke = true; break;
    default: throw std::invalid_argument("Unknown MAC algorithm");
    }
    macCormackVel = macCormackSmoke;
    emitterRadius = .08f;
    sourceStrength = 4.f;
    buoyancy = 3.f;
    smokeDecay = .08f;
    projectIterations = 600;
    sorOmega = 1.9f;
}

MacGridFluidSolver3D::MacGridFluidSolver3D(MacAlgorithm algorithm)
    : MacGridFluidSolver3D(MacGridState3D(32, 4.5f)) {
    parameters_ = MacGridParameters3D(algorithm);
    resetState();
}
float MacGridFluidSolver3D::randomUnit() {
    randomState_ = randomState_ * 1664525u + 1013904223u;
    return (randomState_ >> 8) * (1.0f / 16777216.0f);
}

MacGridFluidSolver3D::MacGridFluidSolver3D(int res)
    : MacGridFluidSolver3D(MacGridState3D(res)) {}

MacGridFluidSolver3D::MacGridFluidSolver3D(MacGridState3D state)
    : state_(std::move(state)) {
    parameters_.projectIterations = state_.resolution() * 2;
    allocateBuffers();
}

void MacGridFluidSolver3D::allocateBuffers() {
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

int MacGridFluidSolver3D::idC(int x, int y, int z) const { return static_cast<int>(state_.idC(x, y, z)); }

int MacGridFluidSolver3D::idX(int x, int y, int z) const { return static_cast<int>(state_.idX(x, y, z)); }

int MacGridFluidSolver3D::idY(int x, int y, int z) const { return static_cast<int>(state_.idY(x, y, z)); }

int MacGridFluidSolver3D::idZ(int x, int y, int z) const { return static_cast<int>(state_.idZ(x, y, z)); }

float MacGridFluidSolver3D::h() const { return state_.cellSize(); }

void MacGridFluidSolver3D::advectScalar(std::span<float> field, float dt) {
    if (advection_) {
        advection_->advectScalar({state_, velocityView(state_), dt, parameters_.cflMc}, field, advectionWorkspace_,
        parameters_.macCormackSmoke ? AdvectionMethod::MacCormack : AdvectionMethod::SemiLagrangian);
    } else {
        cpuAdvection_.advectScalar({state_, velocityView(state_), dt, parameters_.cflMc}, field, advectionWorkspace_,
        parameters_.macCormackSmoke ? AdvectionMethod::MacCormack : AdvectionMethod::SemiLagrangian);
    }
}

void MacGridFluidSolver3D::advectVx(float dt, std::span<const float> vxSrc,
    std::span<const float> vxVel, std::span<const float> vyVel, std::span<const float> vzVel,
    std::span<float> out) {
    if (advection_) {
        advection_->advectFace({state_, {vxVel, vyVel, vzVel}, dt, parameters_.cflMc}, FaceAxis::X,
        vxSrc, out, advectionWorkspace_, parameters_.macCormackVel ? AdvectionMethod::MacCormack : AdvectionMethod::SemiLagrangian);
    } else {
        cpuAdvection_.advectFace({state_, {vxVel, vyVel, vzVel}, dt, parameters_.cflMc}, FaceAxis::X,
        vxSrc, out, advectionWorkspace_, parameters_.macCormackVel ? AdvectionMethod::MacCormack : AdvectionMethod::SemiLagrangian);
    }
}

void MacGridFluidSolver3D::advectVy(float dt, std::span<const float> vySrc,
    std::span<const float> vxVel, std::span<const float> vyVel, std::span<const float> vzVel,
    std::span<float> out) {
    if (advection_) {
        advection_->advectFace({state_, {vxVel, vyVel, vzVel}, dt, parameters_.cflMc}, FaceAxis::Y,
        vySrc, out, advectionWorkspace_, parameters_.macCormackVel ? AdvectionMethod::MacCormack : AdvectionMethod::SemiLagrangian);
    } else {
        cpuAdvection_.advectFace({state_, {vxVel, vyVel, vzVel}, dt, parameters_.cflMc}, FaceAxis::Y,
        vySrc, out, advectionWorkspace_, parameters_.macCormackVel ? AdvectionMethod::MacCormack : AdvectionMethod::SemiLagrangian);
    }
}

void MacGridFluidSolver3D::advectVz(float dt, std::span<const float> vzSrc,
    std::span<const float> vxVel, std::span<const float> vyVel, std::span<const float> vzVel,
    std::span<float> out) {
    if (advection_) {
        advection_->advectFace({state_, {vxVel, vyVel, vzVel}, dt, parameters_.cflMc}, FaceAxis::Z,
        vzSrc, out, advectionWorkspace_, parameters_.macCormackVel ? AdvectionMethod::MacCormack : AdvectionMethod::SemiLagrangian);
    } else {
        cpuAdvection_.advectFace({state_, {vxVel, vyVel, vzVel}, dt, parameters_.cflMc}, FaceAxis::Z,
        vzSrc, out, advectionWorkspace_, parameters_.macCormackVel ? AdvectionMethod::MacCormack : AdvectionMethod::SemiLagrangian);
    }
}

void MacGridFluidSolver3D::enforceVelocityBoundaries() {
    MacGridOperators3D::enforceClosedBoundaries(state_);
}

void MacGridFluidSolver3D::project(float dt, int iterations) {
    ProjectionOptions options;
    options.dt = dt;
    options.relaxation = parameters_.sorOmega;
    options.linearSolve.maxIterations = iterations;
    options.linearSolve.fixedIterations = true;
    lastProjection_ = pressureSolver_ ? pressureSolver_->project(state_, options) : cpuPressure_.project(state_, options);
}

void MacGridFluidSolver3D::emit(float strength, float dt) {
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
                const float kick = parameters_.emitterKickMultiplier * (parameters_.emitterKickBase + parameters_.emitterKickScale * strength);
                if (parameters_.emitterKickMultiplier > 0.f) {
                state_.velocityZ()[idZ(x, y, z)] = std::max(state_.velocityZ()[idZ(x, y, z)], kick);
                state_.velocityZ()[idZ(x, y, z + 1)] = std::max(state_.velocityZ()[idZ(x, y, z + 1)], kick);
                }
                const float j = (randomUnit() - 0.5f) * parameters_.stirStrength * dtScale;
                state_.velocityX()[idX(x, y, z)] += j; state_.velocityX()[idX(x + 1, y, z)] += j;
                state_.velocityY()[idY(x, y, z)] += j; state_.velocityY()[idY(x, y + 1, z)] += j;
            }
}

void MacGridFluidSolver3D::applyBuoyancy(float dtStep, float buoyancyVal) {
#pragma omp parallel for
    for (int z = 1; z < state_.resolution(); ++z) {
        for (int y = 0; y < state_.resolution(); ++y) {
            for (int x = 0; x < state_.resolution(); ++x) {
                // Each face has one writer, including when OpenMP is enabled.
                // Gather both adjacent cells instead of racing cell scatters.
                const float lower = buoyancyVal * state_.density()[idC(x, y, z - 1)] * dtStep;
                const float upper = buoyancyVal * state_.density()[idC(x, y, z)] * dtStep;
                state_.velocityZ()[idZ(x, y, z)] += parameters_.buoyancySplit * lower;
                state_.velocityZ()[idZ(x, y, z)] += parameters_.buoyancySplit * upper;
            }
        }
    }
}

void MacGridFluidSolver3D::applyDissipation(float dt) {
    const float sFactor = std::clamp(1.0f - parameters_.smokeDecay * dt, 0.0f, 1.0f);

#pragma omp parallel for
    for (int i = 0; i < (int)state_.density().size(); ++i) {
        state_.density()[i] *= sFactor;
    }
}

void MacGridFluidSolver3D::step(float dt, float buoyancy, float sourceStrength, int projectIterations) {
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

void MacGridFluidSolver3D::stepNoReflect(float dt, float buoyancy, float sourceStrength, int projectIterations) {
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

float MacGridFluidSolver3D::kineticEnergy() const {
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

float MacGridFluidSolver3D::maxVelocity() const {
    float m = 0.f;
    for (float f : state_.velocityX()) m = std::max(m, std::abs(f));
    for (float f : state_.velocityY()) m = std::max(m, std::abs(f));
    for (float f : state_.velocityZ()) m = std::max(m, std::abs(f));
    return m;
}

float MacGridFluidSolver3D::computeDivergenceNorm() const {
    return MacGridOperators3D::maxDivergence(state_);
}

void MacGridFluidSolver3D::resetState() {
    randomState_ = parameters_.randomSeed;
    state_.reset();
    cpuPressure_.reset();
    if (pressureSolver_) pressureSolver_->reset();
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

int MacGridFluidSolver3D::resolution() const { return state_.resolution(); }

float MacGridFluidSolver3D::boxSize() const noexcept { return state_.boxSize(); }

const ProjectionResult& MacGridFluidSolver3D::lastProjection() const { return lastProjection_; }

void MacGridFluidSolver3D::setBoxSize(float size) {
    if (size != state_.boxSize()) {
        state_.setBoxSize(size);
        resetState();
    }
}

void MacGridFluidSolver3D::setPressureSolver(std::unique_ptr<Projection3D> solver) {
    if (!solver) throw std::invalid_argument("Pressure backend cannot be null");
    pressureSolver_ = std::move(solver);
    lastProjection_ = {};
}
void MacGridFluidSolver3D::setAdvection(std::unique_ptr<Advection3D> method) {
    if (!method) throw std::invalid_argument("Advection backend cannot be null");
    advection_ = std::move(method);
}
void MacGridFluidSolver3D::useDefaultPressureSolver() {
    pressureSolver_.reset();
    cpuPressure_.reset();
    lastProjection_ = {};
}
void MacGridFluidSolver3D::useDefaultAdvection() {
    advection_.reset();
}

void MacGridFluidSolver3D::advance(float dt) {
    if (!std::isfinite(parameters_.emitterKickMultiplier) || parameters_.emitterKickMultiplier < 0.f)
        throw std::invalid_argument("Invalid emitter kick multiplier");
    if (parameters_.reflection)
        step(dt, parameters_.buoyancy, parameters_.sourceStrength, parameters_.projectIterations);
    else
        stepNoReflect(dt, parameters_.buoyancy, parameters_.sourceStrength, parameters_.projectIterations);
}

FluidDiagnostics MacGridFluidSolver3D::diagnostics() const {
    return {kineticEnergy(), maxVelocity(), computeDivergenceNorm()};
}

