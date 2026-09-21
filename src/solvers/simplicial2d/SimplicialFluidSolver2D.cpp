#include "solvers/simplicial2d/SimplicialFluidSolver2D.hpp"

#include "solvers/simplicial2d/SimplicialOperators2D.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

SimplicialFluidSolver2D::SimplicialFluidSolver2D(int resolution, float boxSize)
    : SimplicialFluidSolver2D(SimplicialFluidState2D(resolution, boxSize)) {}

SimplicialFluidSolver2D::SimplicialFluidSolver2D(SimplicialFluidState2D state)
    : state_(std::move(state)), triangleVelocity_(state_.domain().triangleCount()),
      backtracked_(state_.domain().triangleCount()),
      backtrackedVelocity_(state_.domain().triangleCount()),
      forceFlux_(state_.domain().edgeCount()), forceVorticity_(state_.domain().vertexCount()),
      densitySource_(state_.domain().vertexCount()), densityNext_(state_.domain().vertexCount()) {
    randomState_ = parameters_.randomSeed;
    if (state_.domain().isTeapot()) {
        parameters_.emitterCenterS=.49f; parameters_.emitterCenterT=.28f;
        parameters_.emitterRadius=.045f;
    }
}

float SimplicialFluidSolver2D::randomUnit() noexcept {
    randomState_ = randomState_ * 1664525u + 1013904223u;
    return (randomState_ >> 8) * (1.0f / 16777216.0f);
}

void SimplicialFluidSolver2D::advance(float dt) {
    const auto& p = parameters_;
    if (!std::isfinite(dt) || dt <= 0.f || !std::isfinite(p.sourceStrength) || p.sourceStrength < 0.f ||
        !std::isfinite(p.emitterRadius) || p.emitterRadius <= 0.f || p.emitterRadius > 1.f ||
        !std::isfinite(p.emitterCenterS) || p.emitterCenterS < 0.f || p.emitterCenterS > 1.f ||
        !std::isfinite(p.emitterCenterT) || p.emitterCenterT < 0.f || p.emitterCenterT > 1.f ||
        !std::isfinite(p.buoyancy) || !std::isfinite(p.smokeDecay) || p.smokeDecay < 0.f ||
        p.cgIterations < 0 || !std::isfinite(p.cgAbsoluteTolerance) || p.cgAbsoluteTolerance < 0.0 ||
        !std::isfinite(p.cgRelativeTolerance) || p.cgRelativeTolerance < 0.0 ||
        !std::isfinite(p.advectionCfl) || p.advectionCfl <= 0.f || p.advectionCfl > 1.f)
        throw std::invalid_argument("Invalid 2D simplicial simulation time step or parameters");

    double remaining = dt;
    lastSubstepCount_ = 0;
    while (remaining > double(dt) * 1e-7) {
        const auto velocity = SimplicialOperators2D::reconstructTriangleVelocities(
            state_.domain(), state_.edgeFlux());
        double maximumVelocity = 0.0;
        for (const auto sample : velocity)
            maximumVelocity = std::max(maximumVelocity, std::sqrt(dot(sample, sample)));
        if (!std::isfinite(maximumVelocity))
            throw std::runtime_error("Simplicial velocity became non-finite before advection");
        // Include this substep's forcing: an initially resting but strongly
        // heated source must not bypass CFL subdivision.
        const double densityBound = std::max(double(*std::max_element(state_.density().begin(), state_.density().end())),
                                             double(std::min(1.f, 3.f * p.sourceStrength)));
        const double accelerationBound = std::abs(p.buoyancy) * densityBound;
        const double travelLimit = p.advectionCfl * state_.domain().edgeScale();
        const double denominator = maximumVelocity +
            std::sqrt(maximumVelocity * maximumVelocity + 4.0 * accelerationBound * travelLimit);
        const double cflStep = denominator > 0.0 ? 2.0 * travelLimit / denominator : remaining;
        const float substep = static_cast<float>(std::min(remaining, cflStep));
        if (!(substep > 0.f) || !std::isfinite(substep) || ++lastSubstepCount_ > 256)
            throw std::runtime_error("Simplicial CFL subdivision could not make stable progress");
        advanceSubstep(substep);
        remaining -= substep;
    }
}

void SimplicialFluidSolver2D::advanceSubstep(float dt) {
    emitDensity();
    advectVorticity(dt);
    addBuoyancy(dt);
    stirEmitter(dt);
    recoverFluxFromVorticity();
    advectDensity(dt);
    const float decay = std::exp(-dt * parameters_.smokeDecay);
    for (std::size_t i = 0; i < state_.density().size(); ++i)
        state_.density()[i] = decay * densityNext_[i];
}

void SimplicialFluidSolver2D::resetState() {
    state_.reset();
    std::fill(triangleVelocity_.begin(), triangleVelocity_.end(), SimplicialPoint2D{});
    std::fill(backtracked_.begin(), backtracked_.end(), SimplicialPoint2D{});
    std::fill(backtrackedVelocity_.begin(), backtrackedVelocity_.end(), SimplicialPoint2D{});
    std::fill(forceFlux_.begin(), forceFlux_.end(), 0.0);
    std::fill(forceVorticity_.begin(), forceVorticity_.end(), 0.0);
    std::fill(densitySource_.begin(), densitySource_.end(), 0.f);
    std::fill(densityNext_.begin(), densityNext_.end(), 0.f);
    lastRecovery_ = {};
    randomState_ = parameters_.randomSeed;
    lastSubstepCount_ = 0;
}

void SimplicialFluidSolver2D::setBoxSize(float size) {
    const float previous = state_.boxSize();
    state_.setBoxSize(size);
    if (previous != size) resetState();
}

FluidDiagnostics SimplicialFluidSolver2D::diagnostics() const {
    FluidDiagnostics result;
    double energy = 0.0;
    for (std::size_t i = 0; i < state_.domain().edgeCount(); ++i) {
        const double flux = state_.edgeFlux()[i];
        energy += state_.domain().edges()[i].hodge * flux * flux;
    }
    result.kineticEnergy = static_cast<float>(0.5 * energy);
    const auto velocity = SimplicialOperators2D::reconstructTriangleVelocities(
        state_.domain(), state_.edgeFlux());
    for (const auto sample : velocity)
        result.maxVelocity = std::max(result.maxVelocity,
            static_cast<float>(std::sqrt(dot(sample, sample))));
    result.maxDivergence = static_cast<float>(
        SimplicialOperators2D::maxDivergence(state_.domain(), state_.edgeFlux()));
    return result;
}

void SimplicialFluidSolver2D::synchronizeFromStreamFunction() {
    SimplicialOperators2D::fluxFromStreamFunction(
        state_.domain(), state_.streamFunction(), state_.edgeFlux());
    SimplicialOperators2D::vorticityFromFlux(
        state_.domain(), state_.edgeFlux(), state_.vorticity());
}

SimplicialPoint2D SimplicialFluidSolver2D::traceBack(
    SimplicialPoint2D point, float dt,
    std::span<const SimplicialPoint2D> triangleVelocity) const noexcept {
    const auto firstVelocity = SimplicialOperators2D::sampleVelocity(
        state_.domain(), triangleVelocity, point);
    const auto midpoint = state_.domain().clipSegment(point, point - (0.5 * dt) * firstVelocity);
    const auto midpointVelocity = SimplicialOperators2D::sampleVelocity(
        state_.domain(), triangleVelocity, midpoint);
    return state_.domain().clipSegment(point, point - double(dt) * midpointVelocity);
}

void SimplicialFluidSolver2D::advectVorticity(float dt) {
    if (!std::isfinite(dt) || dt < 0.f) throw std::invalid_argument("Invalid simplicial advection step");
    triangleVelocity_ = SimplicialOperators2D::reconstructTriangleVelocities(
        state_.domain(), state_.edgeFlux());
    for (std::size_t i = 0; i < state_.domain().triangleCount(); ++i) {
        const auto center = state_.domain().triangles()[i].circumcenter;
        backtracked_[i] = traceBack(center, dt, triangleVelocity_);
        backtrackedVelocity_[i] = SimplicialOperators2D::sampleVelocity(
            state_.domain(), triangleVelocity_, backtracked_[i]);
    }

    auto updated = state_.vorticity();
    for (std::size_t vertex = 0; vertex < state_.domain().vertexCount(); ++vertex) {
        const auto& meshVertex = state_.domain().vertices()[vertex];
        // Boundary circulation is determined by the recovered free-slip flow,
        // not an independently transported interior dual loop.
        if (meshVertex.boundary) continue;
        double circulation = 0.0;
        const auto& ring = meshVertex.incidentTriangles;
        for (std::size_t i = 0; i < ring.size(); ++i) {
            const int first = ring[i], second = ring[(i + 1) % ring.size()];
            const auto averageVelocity = 0.5 *
                (backtrackedVelocity_[first] + backtrackedVelocity_[second]);
            // d0^T star uses the clockwise orientation of each dual-cell loop.
            circulation += dot(averageVelocity, backtracked_[first] - backtracked_[second]);
        }
        updated[vertex] = circulation;
    }
}

void SimplicialFluidSolver2D::recoverFluxFromVorticity() {
    LinearSolveOptions options;
    options.maxIterations = parameters_.cgIterations;
    options.absoluteTolerance = parameters_.cgAbsoluteTolerance;
    options.relativeTolerance = parameters_.cgRelativeTolerance;
    lastRecovery_ = SimplicialOperators2D::recoverFlux(
        state_.domain(), state_.vorticity(), state_.streamFunction(), state_.edgeFlux(), options);
    SimplicialOperators2D::vorticityFromFlux(state_.domain(), state_.edgeFlux(), forceVorticity_);
    for (std::size_t i = 0; i < state_.domain().vertexCount(); ++i)
        if (state_.domain().vertices()[i].boundary) state_.vorticity()[i] = forceVorticity_[i];
}

void SimplicialFluidSolver2D::emitDensity() {
    if (parameters_.sourceStrength <= 0.f) return;
    const auto center = state_.domain().latticeToWorld(
        parameters_.emitterCenterS, parameters_.emitterCenterT);
    const double radius = parameters_.emitterRadius * state_.boxSize();
    for (std::size_t i = 0; i < state_.domain().vertexCount(); ++i) {
        const auto offset = state_.domain().vertices()[i].position - center;
        const double normalizedSquared = dot(offset, offset) / (radius * radius);
        if (normalizedSquared > 1.0) continue;
        // Same normalized source profile and density range as both 2D MAC solvers.
        const float weight = static_cast<float>(std::exp(-3.5 * normalizedSquared));
        state_.density()[i] = std::max(state_.density()[i],
            std::min(1.f, 3.f * parameters_.sourceStrength * weight));
    }
}

void SimplicialFluidSolver2D::stirEmitter(float dt) {
    if (parameters_.sourceStrength <= 0.f || parameters_.stirStrength <= 0.f) return;
    const auto center = state_.domain().latticeToWorld(
        parameters_.emitterCenterS, parameters_.emitterCenterT);
    const double radius = parameters_.emitterRadius * state_.boxSize();
    // Seed the shear-layer instability with a small random circulation at the
    // source. The emitter, mesh, and forcing are otherwise mirror-symmetric,
    // which keeps a rising plume laminar forever. Applied after vorticity
    // transport so it is not overwritten, and before flux recovery.
    const double dtScale = double(dt) * 60.0; // Per-frame normalization, matching 3D MAC stirring.
    const double scale = double(parameters_.stirStrength) * dtScale * state_.domain().edgeScale();
    for (std::size_t i = 0; i < state_.domain().vertexCount(); ++i) {
        if (state_.domain().vertices()[i].boundary) continue;
        const auto offset = state_.domain().vertices()[i].position - center;
        if (dot(offset, offset) > radius * radius) continue;
        state_.vorticity()[i] += (randomUnit() - 0.5f) * scale;
    }
}

void SimplicialFluidSolver2D::addBuoyancy(float dt) {
    for (std::size_t i = 0; i < state_.domain().edgeCount(); ++i) {
        const auto& edge = state_.domain().edges()[i];
        const auto tangent = state_.domain().vertices()[edge.second].position -
                             state_.domain().vertices()[edge.first].position;
        const double smoke = 0.5 * (state_.density()[edge.first] + state_.density()[edge.second]);
        // Integral of the vertical force through the edge's oriented left normal.
        forceFlux_[i] = parameters_.buoyancy * smoke * tangent.x;
    }
    SimplicialOperators2D::vorticityFromFlux(
        state_.domain(), forceFlux_, forceVorticity_);
    for (std::size_t i = 0; i < state_.domain().vertexCount(); ++i)
        if (!state_.domain().vertices()[i].boundary)
            state_.vorticity()[i] += dt * forceVorticity_[i];
}

void SimplicialFluidSolver2D::advectDensity(float dt) {
    densitySource_.assign(state_.density().begin(), state_.density().end());
    triangleVelocity_ = SimplicialOperators2D::reconstructTriangleVelocities(
        state_.domain(), state_.edgeFlux());
    for (std::size_t i = 0; i < state_.domain().vertexCount(); ++i) {
        const auto origin = traceBack(state_.domain().vertices()[i].position, dt, triangleVelocity_);
        densityNext_[i] = SimplicialOperators2D::sampleVertexScalar(
            state_.domain(), densitySource_, origin);
    }
}
