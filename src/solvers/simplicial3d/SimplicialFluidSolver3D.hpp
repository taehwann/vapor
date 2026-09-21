#pragma once
#include "SimplicialBoundary3D.hpp"
#include "SimplicialMesh3D.hpp"
#include "SimplicialCompute3D.hpp"
#include "common/LinearSolveOptions.hpp"
#include <functional>

struct SimplicialFluidParameters3D {
    bool emitterEnabled = true;
    float emitterX = .5f, emitterY = .5f, emitterZ = .2f;
    float advectionCfl = .4f;
    // Local rotating body force; strength is an acceleration bound (world units/s^2).
    float stirStrength = 0, stirRadius = .18f, stirFrequency = 1;
    float sourceStrength = 4, emitterRadius = .12f, buoyancy = 3, smokeDecay = .08f;
    int cgIterations = 600;
    double cgAbsoluteTolerance = 1e-10, cgRelativeTolerance = 1e-9;
};
struct SimplicialDiagnostics3D {
    double maxDivergence = 0, maxBoundaryFlux = 0, kineticEnergy = 0;
    double maxSpeed = 0, recoveryCurlError = 0, smokeMass = 0;
};
struct SimplicialTimings3D {
    double emittedMass = 0;
    int sourceVertices = 0;
    double emissionMs=0, advectionMs=0, forcesMs=0, recoveryMs=0, smokeMs=0;
    int substeps=0;
};

class SimplicialFluidSolver3D {
  public:
    void setComputeBackend(SimplicialCompute3D* backend) { compute_=backend; }
    // The path constructor retains the box fixture; pass Mesh::loadDomain for
    // a general well-centered tetrahedral domain such as bunny-fluid.tet.
    explicit SimplicialFluidSolver3D(const std::filesystem::path& meshPath);
    explicit SimplicialFluidSolver3D(simplicial3d::Mesh mesh);
    const simplicial3d::Mesh &domain() const { return mesh_; }
    const simplicial3d::DualMesh &dual() const { return dual_; }
    const std::vector<simplicial3d::Point> &dualVelocity() const { return dualVelocity_; }
    SimplicialFluidParameters3D &parameters() { return parameters_; }
    const SimplicialFluidParameters3D &parameters() const { return parameters_; }
    const std::vector<double> &flux() const { return flux_; }
    const std::vector<double> &vorticity() const { return omega_; }
    const std::vector<double> &potential() const { return phi_; }
    const std::vector<double> &boundaryCirculation() const { return boundaryCirculation_; }
    const std::vector<float> &density() const { return density_; }
    const std::vector<float> &vertexDensity() const { return vertexDensity_; }
    const LinearSolveResult &lastRecovery() const { return recovery_; }
    const LinearSolveResult &lastBoundaryReconstruction() const { return boundarySolve_; }
    float boxSize() const { return float(mesh_.size); }
    void setBoxSize(float size);
    void resetState();
    void advance(float dt);
    SimplicialDiagnostics3D diagnostics() const;
    const SimplicialTimings3D& timings() const { return timings_; }
    // Manufactured-field and embedding API, with explicit cochain dimensions.
    void setPotential(std::span<const double> potential);
    void setDensity(const std::function<double(simplicial3d::Point)> &field);
    std::vector<double> curl(std::span<const double> flux) const;
    std::vector<double> applyLaplacian(std::span<const double> potential) const;
    void recoverFluxFromVorticity();
    void advectVorticity(double dt);
    void addForce(const std::function<simplicial3d::Point(simplicial3d::Point)> &force, double dt);
    // Section 4.5: solve once for the harmonic extension of balanced prescribed
    // boundary fluxes; this field is then included in all velocity sampling.
    void setBoundaryFlux(std::span<const double> boundaryFlux);
    const std::vector<double> &harmonicFlux() const { return harmonic_; }
    simplicial3d::Point sampleVelocity(simplicial3d::Point point) const;
    double sampleDensity(simplicial3d::Point point) const;

  private:
    friend class GpuSimplicial3D;
    SimplicialCompute3D* compute_=nullptr;
    void buildPreconditioner();
    LinearSolveResult solve(std::span<const double> rhs, std::vector<double> &x) const;
    std::vector<double> laplacian(std::span<const double> potential) const;
    void reconstructVelocity(bool initializeBoundary = false);
    void reconstructDensity();
    simplicial3d::Point trace(simplicial3d::Point p, double dt) const;
    void substep(double dt);
    void emitSmoke();
    void applyBuoyancy(double dt);
    void stirEmitter(double dt);
    void advectAndDecaySmoke(double dt);
    simplicial3d::Mesh mesh_;
    simplicial3d::DualMesh dual_;
    simplicial3d::BoundaryCirculation boundary_;
    SimplicialFluidParameters3D parameters_;
    std::vector<double> flux_, omega_, phi_, boundaryCirculation_, diagonal_, harmonic_;
    std::vector<float> density_;
    std::vector<simplicial3d::Point> dualVelocity_;
    std::vector<int> boundarySeeds_;
    struct DensitySample { std::array<int,4> vertices; std::array<double,4> weights; };
    std::vector<DensitySample> faceSamples_,dualSamples_;
    std::vector<int> dualTets_;
    std::vector<int> primalTets_;
    std::vector<float> vertexDensity_; // Passive scalar, advected on primal vertices.
    LinearSolveResult recovery_, boundarySolve_;
    SimplicialTimings3D timings_;
    double stirringTime_ = 0;
    std::vector<std::vector<std::pair<int,double>>> incompleteLower_;
    std::vector<double> incompleteDiagonal_;
};
