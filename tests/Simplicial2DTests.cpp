#include "solvers/simplicial2d/SimplicialFluidSolver2D.hpp"
#include "solvers/simplicial2d/SimplicialOperators2D.hpp"
#include "solvers/simplicial2d/SimplicialRasterizer2D.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <numbers>
#include <numeric>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace {
void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}
void near(double actual, double expected, double tolerance, const char* message) {
    if (!std::isfinite(actual) || std::abs(actual - expected) > tolerance)
        throw std::runtime_error(std::string(message) + ": " + std::to_string(actual) +
                                 " vs " + std::to_string(expected));
}
template<class Exception = std::invalid_argument, class Function>
void rejects(Function&& function) {
    try { function(); } catch (const Exception&) { return; }
    throw std::runtime_error("Expected invalid input to be rejected");
}

void meshTopology() {
    for (int n : {2, 5, 12}) {
        SimplicialMesh2D mesh(n, 3.f);
        require(mesh.dimension() == Dimension::D2, "Simplicial dimension tag");
        require(mesh.vertexCount() == std::size_t((n + 1) * (n + 1)), "Vertex count");
        require(mesh.edgeCount() == std::size_t(3 * n * n + 2 * n), "Edge count");
        require(mesh.triangleCount() == std::size_t(2 * n * n), "Triangle count");
        int boundaryEdges = 0;
        for (const auto& edge : mesh.edges()) {
            require(edge.first < edge.second, "Canonical edge orientation");
            require(edge.incidentTriangles[0] >= 0, "Edge without incident triangle");
            require(edge.hodge > 0.0 && std::isfinite(edge.hodge), "Positive Hodge star");
            boundaryEdges += edge.boundary();
        }
        require(boundaryEdges == 4 * n, "Boundary edge count");
        double area = 0.0;
        near(mesh.bounds().max.x, 3.0, 1e-12, "Square width");
        near(mesh.bounds().max.y, 3.0, 1e-12, "Square height");
        for (const auto& triangle : mesh.triangles()) {
            require(triangle.area > 0, "Positive triangle area");
            area += triangle.area;
            const auto bary = mesh.barycentric(
                static_cast<int>(&triangle - mesh.triangles().data()), triangle.circumcenter);
            for (double coordinate : bary) require(coordinate >= -1e-12, "Non-obtuse triangle");
        }
        near(area, 9.0, 1e-12, "Mesh covers square");
        for (int j=0; j<=40; ++j) for (int i=0; i<=40; ++i) {
            const auto p=mesh.latticeToWorld(i/40.0,j/40.0);
            const auto bary=mesh.barycentric(mesh.locateTriangle(p),p);
            for (double w:bary) require(w>=-1e-12 && w<=1+1e-12,"Point location outside triangle");
        }

        std::vector<double> potential(mesh.vertexCount()), flux(mesh.edgeCount());
        for (std::size_t i = 0; i < potential.size(); ++i)
            potential[i] = std::sin(.37 * i) + .01 * i;
        SimplicialOperators2D::fluxFromStreamFunction(mesh, potential, flux);
        near(SimplicialOperators2D::maxDivergence(mesh, flux), 0.0, 2e-12,
             "d1*d0 must vanish");
    }
    rejects([] { SimplicialMesh2D mesh(1); });
    rejects([] { SimplicialMesh2D mesh(4, 0.f); });
    rejects([] { SimplicialMesh2D mesh(4, std::numeric_limits<float>::infinity()); });
}

void velocityReconstruction() {
    SimplicialMesh2D mesh(10, 2.f);
    const SimplicialPoint2D expected{.7, -.35};
    std::vector<double> flux(mesh.edgeCount());
    for (std::size_t i = 0; i < mesh.edgeCount(); ++i) {
        const auto& edge = mesh.edges()[i];
        const auto tangent = mesh.vertices()[edge.second].position -
                             mesh.vertices()[edge.first].position;
        flux[i] = dot(expected, {-tangent.y, tangent.x});
    }
    const auto velocity = SimplicialOperators2D::reconstructTriangleVelocities(mesh, flux);
    for (const auto value : velocity) {
        near(value.x, expected.x, 2e-14, "Reconstructed velocity x");
        near(value.y, expected.y, 2e-14, "Reconstructed velocity y");
    }
}

void recovery() {
    SimplicialFluidSolver2D solver(20, 2.f);
    const auto& mesh = solver.state().domain();
    for (int j = 0; j <= mesh.resolution(); ++j) for (int i = 0; i <= mesh.resolution(); ++i) {
        const auto position = mesh.vertices()[mesh.vertexIndex(i,j)].position;
        const double s = position.x / mesh.boxSize();
        const double t = position.y / mesh.boxSize();
        solver.state().streamFunction()[mesh.vertexIndex(i, j)] =
            std::sin(std::numbers::pi_v<double> * s) *
            std::sin(2.0 * std::numbers::pi_v<double> * t);
    }
    solver.synchronizeFromStreamFunction();
    const std::vector<double> expectedFlux(
        solver.state().edgeFlux().begin(), solver.state().edgeFlux().end());
    const std::vector<double> expectedVorticity(
        solver.state().vorticity().begin(), solver.state().vorticity().end());
    std::fill(solver.state().streamFunction().begin(), solver.state().streamFunction().end(), 0.0);
    std::fill(solver.state().edgeFlux().begin(), solver.state().edgeFlux().end(), 0.0);
    solver.recoverFluxFromVorticity();
    require(solver.lastRecovery().converged, "CG flux recovery did not converge");
    require(solver.lastRecovery().finalResidual < 1e-8, "CG residual too large");
    for (std::size_t i = 0; i < expectedFlux.size(); ++i)
        near(solver.state().edgeFlux()[i], expectedFlux[i], 2e-8, "Recovered flux");
    std::vector<double> recoveredVorticity(mesh.vertexCount());
    SimplicialOperators2D::vorticityFromFlux(
        mesh, solver.state().edgeFlux(), recoveredVorticity);
    for (std::size_t i = 0; i < mesh.vertexCount(); ++i)
        if (!mesh.vertices()[i].boundary)
            near(recoveredVorticity[i], expectedVorticity[i], 2e-8, "Curl of recovered flux");
    near(SimplicialOperators2D::maxDivergence(mesh, solver.state().edgeFlux()),
         0.0, 2e-12, "Recovered flux divergence");
    for (const auto& edge : mesh.edges()) if (edge.boundary())
        near(solver.state().edgeFlux()[&edge - mesh.edges().data()], 0.0, 0.0,
             "No-transfer boundary flux");
}

void circulationIdentity() {
    SimplicialFluidSolver2D solver(18, 2.f);
    const auto& mesh = solver.state().domain();
    for (int j = 0; j <= mesh.resolution(); ++j) for (int i = 0; i <= mesh.resolution(); ++i) {
        const auto position = mesh.vertices()[mesh.vertexIndex(i,j)].position;
        const double s = position.x / mesh.boxSize();
        const double t = position.y / mesh.boxSize();
        solver.state().streamFunction()[mesh.vertexIndex(i, j)] =
            std::sin(std::numbers::pi_v<double> * s) *
            std::sin(std::numbers::pi_v<double> * t) *
            (1.0 + .2 * std::sin(2.0 * std::numbers::pi_v<double> * s));
    }
    solver.synchronizeFromStreamFunction();
    const std::vector<double> before(
        solver.state().vorticity().begin(), solver.state().vorticity().end());
    solver.advectVorticity(0.f);
    for (std::size_t i = 0; i < mesh.vertexCount(); ++i) {
        if (mesh.vertices()[i].boundary) near(solver.state().vorticity()[i], before[i], 0.0, "Boundary circulation remains consistent");
            else near(solver.state().vorticity()[i], before[i], 2e-13,
                  "Backtracked dual-loop circulation at dt=0");
    }

    const auto integralMagnitude = [&] {
        double sum = 0.0;
        for (std::size_t i = 0; i < mesh.vertexCount(); ++i)
            if (!mesh.vertices()[i].boundary) sum += std::abs(solver.state().vorticity()[i]);
        return sum;
    };
    // Small CFL must still transport a non-steady vortex. Endpoint quadrature
    // with piecewise-constant samples can telescope to exactly zero change.
    solver.advectVorticity(1e-4f);
    double change = 0.0;
    for (std::size_t i = 0; i < mesh.vertexCount(); ++i)
        if (!mesh.vertices()[i].boundary)
            change = std::max(change, std::abs(solver.state().vorticity()[i] - before[i]));
    std::cout << "small-step circulation change=" << change << '\n';
    require(change > 1e-8, "Small-CFL circulation transport is frozen");
    // Use a one-signed interior vortex for the retention benchmark. The sine
    // stream function above has mixed-sign curl extending to truncated wall
    // cells; its sum(abs(Omega)) is not a discrete Kelvin invariant.
    solver.resetState();
    const auto center = mesh.latticeToWorld(.5, .5);
    for (std::size_t i = 0; i < mesh.vertexCount(); ++i) {
        const auto delta = mesh.vertices()[i].position - center;
        if (!mesh.vertices()[i].boundary)
            solver.state().vorticity()[i] = mesh.edgeScale() * mesh.edgeScale() *
                std::exp(-dot(delta, delta) / .04);
    }
    // Linear precision of the continuous interpolant inside complete dual cells.
    std::vector<SimplicialPoint2D> affine(mesh.triangleCount());
    for (std::size_t i = 0; i < affine.size(); ++i) {
        const auto p = mesh.triangles()[i].circumcenter;
        affine[i] = {1.0 + .3 * p.x - .2 * p.y, -.4 + .1 * p.x + .5 * p.y};
    }
    for (int j = 2; j < 8; ++j) for (int i = 2; i < 8; ++i) {
        const auto p = mesh.latticeToWorld((i + .37) / 10., (j + .23) / 10.);
        const auto sampled = SimplicialOperators2D::sampleVelocity(mesh, affine, p);
        near(sampled.x, 1.0 + .3 * p.x - .2 * p.y, 1e-12, "Dual interpolation affine x");
        near(sampled.y, -.4 + .1 * p.x + .5 * p.y, 1e-12, "Dual interpolation affine y");
    }
    solver.recoverFluxFromVorticity();
    const double initialMagnitude = integralMagnitude();
    for (int step = 0; step < 40; ++step) {
        solver.advectVorticity(.01f);
        solver.recoverFluxFromVorticity();
        require(solver.lastRecovery().converged, "Circulation recovery did not converge");
    }
    const double retained = integralMagnitude() / initialMagnitude;
    std::cout << "circulation magnitude retained=" << retained << '\n';
    require(retained > .98 && retained < 1.02, "Interior circulation drift exceeds two percent");
}

void simulation() {
    SimplicialFluidSolver2D solver(24, 4.f);
    solver.parameters().sourceStrength = 0.f;
    for (int i = 0; i < 4; ++i) solver.advance(1.f / 60.f);
    near(solver.diagnostics().kineticEnergy, 0.0, 0.0, "Rest energy");
    solver.parameters().sourceStrength = 4.f;
    for (int i = 0; i < 60; ++i) {
        solver.advance(1.f / 60.f);
        require(solver.lastRecovery().converged, "Simulation recovery did not converge");
        require(solver.diagnostics().maxDivergence < 1e-10f, "Simulation lost incompressibility");
        for (float density : solver.state().density())
            require(std::isfinite(density) && density >= 0.f, "Invalid advected dye");
    }
    require(solver.diagnostics().kineticEnergy > 1e-5f, "Buoyancy did not produce simplicial flow");
    require(solver.diagnostics().maxVelocity > 1e-3f, "Simplicial flow has no velocity");

    SimplicialRasterizer2D presentation(solver.state());
    const auto image = presentation.renderData();
    require(image.density.data() != solver.state().density().data(), "Presentation must use a physical raster");
    require(image.width == 25 && image.height == 25, "Physical raster image extent");
    solver.resetState();
    for (float density : solver.state().density()) near(density, 0.0, 0.0, "Reset dye");
    rejects([&] { solver.advance(0.f); });
}

void matchedDensitySource() {
    SimplicialFluidSolver2D solver(32,1.f);
    auto& p=solver.parameters();
    p.buoyancy=0; p.smokeDecay=0; p.sourceStrength=.2f; p.emitterRadius=.15f; p.stirStrength=0.f;
    const auto& mesh=solver.domain();
    for(int step=0;step<12;++step) solver.advance(1.f/60);
    for(size_t i=0;i<mesh.vertexCount();++i) {
        const auto position=mesh.vertices()[i].position;
        const double x=position.x-.5,y=position.y-.2;
        const double q=(x*x+y*y)/(double(p.emitterRadius)*p.emitterRadius);
        const double expected=q>1?0:std::min(1.0,3.0*p.sourceStrength*std::exp(-3.5*q));
        near(solver.state().density()[i],expected,2e-6,"Matched Gaussian source must not accumulate");
    }
    p.sourceStrength=4;solver.advance(1.f/60);
    for(float density:solver.state().density()) require(density>=0 && density<=1,"Source density outside MAC range");
    const auto before=std::vector<float>(solver.state().density().begin(),solver.state().density().end());
    p.sourceStrength=0;p.smokeDecay=.3f;solver.advance(1.f/60);
    for(size_t i=0;i<before.size();++i)
        near(solver.state().density()[i],before[i]*std::exp(-.3/60),2e-6,"Source-off density decay differs from MAC");
    SimplicialRasterizer2D raster(solver.state());
    require(raster.renderData().historicalPalette,"Simplicial density must use the 2D MAC palette");
}

void presentationCoordinates() {
    SimplicialFluidState2D state(24, 4.f);
    const auto& mesh = state.domain();
    for (std::size_t i = 0; i < mesh.vertexCount(); ++i) {
        const double normalizedX = mesh.vertices()[i].position.x / state.boxSize();
        state.density()[i] = static_cast<float>(
            std::exp(-200.0 * (normalizedX - .75) * (normalizedX - .75)));
    }
    SimplicialRasterizer2D presentation(state);
    const auto image = presentation.renderData();
    const auto maximumColumn = [&](int row) {
        const auto first = image.density.begin() + std::size_t(row) * image.width;
        return static_cast<int>(std::max_element(first, first + image.width) - first);
    };
    const int lower = maximumColumn(image.height / 4);
    const int upper = maximumColumn(3 * image.height / 4);
    require(std::abs(lower - upper) <= 1,
            "A physical vertical feature was sheared diagonally by presentation");
    near(double(image.worldWidth) / image.worldHeight, 1.0,
         2e-7, "Physical raster aspect ratio");
}

void cflSubdivision() {
    SimplicialFluidSolver2D solver(64, 4.f);
    solver.parameters().sourceStrength = 0.f;
    solver.parameters().buoyancy = 0.f;
    const auto& mesh = solver.state().domain();
    for (int j = 0; j <= mesh.resolution(); ++j) for (int i = 0; i <= mesh.resolution(); ++i) {
        const auto position = mesh.vertices()[mesh.vertexIndex(i,j)].position;
        const double s = position.x / mesh.boxSize();
        const double t = position.y / mesh.boxSize();
        solver.state().streamFunction()[mesh.vertexIndex(i, j)] =
            4.0 * std::sin(std::numbers::pi_v<double> * s) *
            std::sin(std::numbers::pi_v<double> * t);
    }
    solver.synchronizeFromStreamFunction();
    solver.advance(1.f / 30.f);
    require(solver.lastSubstepCount() > 1, "High-CFL frame was not subdivided");
    require(std::isfinite(solver.diagnostics().kineticEnergy), "CFL-subdivided energy became non-finite");
    require(solver.lastRecovery().converged, "CFL-subdivided recovery failed");
    solver.resetState();
    solver.parameters().sourceStrength=10;solver.parameters().buoyancy=100;
    solver.advance(1.f/30);
    require(solver.lastSubstepCount()>1,"Buoyancy from rest bypasses CFL subdivision");
}

void wallAndRefinement() {
    double previousError = 0;
    for (int n : {16,32,64}) {
        SimplicialFluidSolver2D solver(n,1.f);
        const auto& mesh=solver.domain();
        constexpr double pi=std::numbers::pi_v<double>;
        for (size_t i=0;i<mesh.vertexCount();++i) {
            const auto p=mesh.vertices()[i].position;
            solver.state().streamFunction()[i]=mesh.vertices()[i].boundary ? 0.0 : std::sin(pi*p.x)*std::sin(pi*p.y);
        }
        solver.synchronizeFromStreamFunction();
        const auto velocity=SimplicialOperators2D::reconstructTriangleVelocities(mesh,solver.state().edgeFlux());
        double error=0;
        for (int j=0;j<101;++j) {
            const double t=(j+.37)/102;
            for (const auto p : {SimplicialPoint2D{0,t}, {1,t}, {t,0}, {t,1}}) {
                const auto v=SimplicialOperators2D::sampleVelocity(mesh,velocity,p);
                if(p.x==0 || p.x==1) near(v.x,0,1e-11,"Vertical wall leakage");
                else near(v.y,0,1e-11,"Horizontal wall leakage");
                const SimplicialPoint2D expected{-pi*std::sin(pi*p.x)*std::cos(pi*p.y),pi*std::cos(pi*p.x)*std::sin(pi*p.y)};
                error+=dot(v-expected,v-expected);
            }
        }
        error=std::sqrt(error/404);
        std::cout<<"wall velocity n="<<n<<" rms_error="<<error<<'\n';
        if(previousError>0) require(error<previousError*.8,"Wall velocity does not improve under refinement");
        previousError=error;
        const auto slip=SimplicialOperators2D::sampleVelocity(mesh,velocity,{0,.5});
        require(std::abs(slip.y)>2.5,"Wall incorrectly clamps tangential velocity");
        // Probe almost exactly on every interior dual edge (previous 0/0 case).
        for(const auto& e:mesh.edges()) if(!e.boundary()) {
            auto p=.5*(mesh.triangles()[e.incidentTriangles[0]].circumcenter+mesh.triangles()[e.incidentTriangles[1]].circumcenter);
            p.x+=1e-13;
            const auto v=SimplicialOperators2D::sampleVelocity(mesh,velocity,p);
            require(std::isfinite(v.x)&&std::isfinite(v.y),"Nonfinite velocity beside a dual edge");
        }
        solver.recoverFluxFromVorticity();
        std::vector<double> curl(mesh.vertexCount());
        SimplicialOperators2D::vorticityFromFlux(mesh,solver.state().edgeFlux(),curl);
        for(size_t i=0;i<curl.size();++i) if(mesh.vertices()[i].boundary)
            near(solver.state().vorticity()[i],curl[i],1e-12,"Boundary circulation inconsistent with flux");
    }
    // Fine meshes/small physical boxes must not fail an absolute determinant threshold.
    SimplicialMesh2D fine(64,.01f);
    std::vector<double> zero(fine.edgeCount());
    const auto velocities=SimplicialOperators2D::reconstructTriangleVelocities(fine,zero);
    for(auto v:velocities) require(v.x==0 && v.y==0,"Fine-grid rest velocity");
}

void temporalRefinement() {
    const auto run=[](double dt, double centerX) {
        SimplicialFluidSolver2D solver(32,1.f);
        auto& p=solver.parameters();p.sourceStrength=0;p.buoyancy=0;p.smokeDecay=0;
        const auto& mesh=solver.domain();
        for(size_t i=0;i<mesh.vertexCount();++i) if(!mesh.vertices()[i].boundary) {
            const auto a=mesh.vertices()[i].position-SimplicialPoint2D{centerX,.45};
            const auto b=mesh.vertices()[i].position-SimplicialPoint2D{.65,.55};
            solver.state().vorticity()[i]=mesh.edgeScale()*mesh.edgeScale()*
                (10*std::exp(-dot(a,a)/.008)-8*std::exp(-dot(b,b)/.014));
        }
        solver.recoverFluxFromVorticity();
        for(int i=0;i<int(std::lround(.2/dt));++i) {
            solver.advance(float(dt));
            require(solver.lastRecovery().converged,"Vortex CG did not converge");
        }
        return std::vector<double>(solver.state().edgeFlux().begin(),solver.state().edgeFlux().end());
    };
    for(double centerX : {.3,.1}) {
    const auto coarse=run(.02,centerX), medium=run(.01,centerX), fine=run(.005,centerX), reference=run(.0025,centerX);
    const auto error=[&](const auto& field) {double sum=0;for(size_t i=0;i<field.size();++i)sum+=(field[i]-reference[i])*(field[i]-reference[i]);return std::sqrt(sum/field.size());};
    const double a=error(coarse),b=error(medium),c=error(fine);
    std::cout<<"vortex center x="<<centerX<<" temporal errors="<<a<<", "<<b<<", "<<c<<'\n';
    require(b<a && c<b,"Vortex transport fails timestep refinement");
    }
}

void stability(int frames, const char* snapshotPrefix) {
    SimplicialFluidSolver2D solver(64, 4.f);
    constexpr float dt = 1.f / 30.f;
    for (int frame = 0; frame < frames; ++frame) {
        solver.advance(dt);
        const auto diagnostics = solver.diagnostics();
        const auto density = solver.state().density();
        const auto maximumDensity = *std::max_element(density.begin(), density.end());
        const auto maximumVorticity = *std::max_element(
            solver.state().vorticity().begin(), solver.state().vorticity().end(),
            [](double a, double b) { return std::abs(a) < std::abs(b); });
        if ((frame + 1) % 30 == 0) {
            std::cout << "t=" << (frame + 1) * dt
                      << " energy=" << diagnostics.kineticEnergy
                      << " velocity=" << diagnostics.maxVelocity
                      << " density=" << maximumDensity
                      << " vorticity=" << std::abs(maximumVorticity)
                      << " residual=" << solver.lastRecovery().finalResidual << '\n';
        }
        require(std::isfinite(diagnostics.kineticEnergy), "Long-run energy became non-finite");
        require(std::isfinite(diagnostics.maxVelocity), "Long-run velocity became non-finite");
        require(std::isfinite(maximumDensity), "Long-run dye became non-finite");
        require(std::isfinite(maximumVorticity), "Long-run vorticity became non-finite");
        require(solver.lastRecovery().converged, "Long-run CG recovery failed");
        require(diagnostics.maxDivergence < 1e-10f, "Long-run flux lost incompressibility");
        if (snapshotPrefix && ((frame + 1) == 30 || (frame + 1) == 90 || (frame + 1) == frames)) {
            SimplicialRasterizer2D raster(solver.state());
            const auto image = raster.renderData();
            std::ofstream file(std::string(snapshotPrefix) + "-" + std::to_string(frame + 1) + ".ppm", std::ios::binary);
            file << "P6\n" << image.width << ' ' << image.height << "\n255\n";
            for (int y = image.height - 1; y >= 0; --y) for (int x = 0; x < image.width; ++x) {
                const unsigned char value = static_cast<unsigned char>(255.f * (1.f - std::exp(-2.f * image.density[x + image.width * y])));
                for (int c = 0; c < 3; ++c) file.put(static_cast<char>(value));
            }
            require(bool(file), "Could not write plume snapshot");
        }
    }
}
}

int main(int argc, char** argv) {
    try {
        const std::string_view name = argc > 1 ? argv[1] : "";
        if (name == "simplicial2d_mesh") meshTopology();
        else if (name == "simplicial2d_velocity") velocityReconstruction();
        else if (name == "simplicial2d_recovery") recovery();
        else if (name == "simplicial2d_circulation") circulationIdentity();
        else if (name == "simplicial2d_simulation") simulation();
        else if (name == "simplicial2d_presentation") { presentationCoordinates(); matchedDensitySource(); }
        else if (name == "simplicial2d_cfl") cflSubdivision();
        else if (name == "simplicial2d_walls") wallAndRefinement();
        else if (name == "simplicial2d_refinement") temporalRefinement();
        else if (name == "simplicial2d_stability")
            stability(argc > 2 ? std::stoi(argv[2]) : 30, argc > 3 ? argv[3] : nullptr);
        else throw std::invalid_argument("Unknown 2D simplicial test case");
        std::cout << "PASS " << name << '\n';
    } catch (const std::exception& error) {
        std::cerr << "FAIL 2D simplicial: " << error.what() << '\n';
        return 1;
    }
}
