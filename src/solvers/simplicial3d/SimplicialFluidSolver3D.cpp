#include "SimplicialFluidSolver3D.hpp"
#include <algorithm>
#include <cmath>
#include <chrono>
#include <map>
#include "common/ParallelSamples.hpp"
#include <numeric>
#include <stdexcept>

using namespace simplicial3d;
namespace {
double inner(std::span<const double> a, std::span<const double> b) {
    return std::inner_product(a.begin(), a.end(), b.begin(), 0.);
}
double norm(std::span<const double> a) {
    double s = 0;
    for (double x : a)
        s = std::max(s, std::abs(x));
    return s;
}
void validDt(double dt) {
    if (!std::isfinite(dt) || dt < 0)
        throw std::invalid_argument("Timestep must be finite and nonnegative");
}
} // namespace
SimplicialFluidSolver3D::SimplicialFluidSolver3D(const std::filesystem::path& meshPath)
    : SimplicialFluidSolver3D(Mesh::loadBox(meshPath, 2.f)) {
    parameters_.emitterRadius = .15f;
    parameters_.cgIterations = 1000;
}
SimplicialFluidSolver3D::SimplicialFluidSolver3D(Mesh mesh)
    : mesh_(std::move(mesh)), dual_(mesh_), boundary_(mesh_, dual_) {
    if(!mesh_.boxDomain) {
        parameters_.sourceStrength=1;
        parameters_.cgIterations=1500;
    }
    diagonal_.assign(mesh_.edges.size(), 0);
    for (size_t f = 0; f < mesh_.faces.size(); ++f)
        for (int e : mesh_.faces[f].edges)
            diagonal_[e] += mesh_.star2[f];
    for (size_t e = 0; e < mesh_.edges.size(); ++e) {
        for (int v : mesh_.edges[e].vertices)
            if (!mesh_.boundaryVertices[v])
                diagonal_[e] += mesh_.star1[e] * mesh_.star1[e] / mesh_.star0[v];
    }
    if(!mesh_.boxDomain)buildPreconditioner();
    boundarySeeds_.resize(dual_.vertices.size());
    for(size_t i=0;i<dual_.vertices.size();++i) {
        if(dual_.vertices[i].tet>=0){boundarySeeds_[i]=dual_.vertices[i].tet;continue;}
        double best=1e300;int nearest=0;
        for(int t=0;t<int(mesh_.tets.size());++t) {
            auto d=dual_.vertices[i].position-mesh_.tets[t].center;
            if(dot(d,d)<best){best=dot(d,d);nearest=t;}
        }
        boundarySeeds_[i]=nearest;
    }
    if(!mesh_.boxDomain) {
        auto sampler=[&](Point point,int hint=-1) {
            int ti=mesh_.locate(point,hint);
            if(ti<0)throw std::invalid_argument("Fixed density sample outside mesh");
            auto w=mesh_.barycentric(ti,point);double total=0;
            for(auto& x:w){x=std::max(0.,x);total+=x;}
            for(auto& x:w)x/=total;
            return DensitySample{mesh_.tets[ti].vertices,w};
        };
        for(const auto& f:mesh_.faces)faceSamples_.push_back(sampler(f.center,f.tets[0]>=0?f.tets[0]:f.tets[1]));
        for(const auto& v:dual_.vertices) {
            int ti=mesh_.locate(v.position,v.tet);
            dualTets_.push_back(ti);dualSamples_.push_back(sampler(v.position,ti));
        }
        primalTets_.resize(mesh_.vertices.size());
        for(int t=0;t<int(mesh_.tets.size());++t)for(int v:mesh_.tets[t].vertices)primalTets_[v]=t;
    }
    resetState();
}
void SimplicialFluidSolver3D::resetState() {
    timings_ = {};
    stirringTime_ = 0;
    flux_.assign(mesh_.faces.size(), 0);
    phi_.assign(mesh_.edges.size(), 0);
    omega_ = phi_;
    harmonic_ = flux_;
    vertexDensity_.assign(mesh_.vertices.size(), 0);
    boundaryCirculation_ = phi_;
    density_.assign(mesh_.tets.size(), 0);
    dualVelocity_.assign(dual_.vertices.size(), {});
    recovery_ = {};
    boundarySolve_ = {};
    dual_.prepareInterpolation(dualVelocity_);
}
void SimplicialFluidSolver3D::setBoxSize(float s) {
    auto p = parameters_;
    auto* backend=compute_;
    *this = SimplicialFluidSolver3D(mesh_.scaledBox(s));
    parameters_ = p;
    compute_=backend;
}
std::vector<double> SimplicialFluidSolver3D::curl(std::span<const double> u) const {
    if (u.size() != mesh_.faces.size())
        throw std::invalid_argument("Flux size mismatch");
    std::vector<double> c(u.begin(), u.end());
    for (size_t f = 0; f < c.size(); ++f)
        c[f] *= mesh_.star2[f];
    return mesh_.d1Transpose(c);
}
std::vector<double> SimplicialFluidSolver3D::applyLaplacian(std::span<const double> input) const {
    if (input.size() != mesh_.edges.size())
        throw std::invalid_argument("Potential size mismatch");
    std::vector<double> x(input.begin(), input.end());
    for (size_t e = 0; e < x.size(); ++e)
        if (mesh_.edges[e].boundary)
            x[e] = 0;
    auto y = laplacian(x);
    for (size_t e = 0; e < x.size(); ++e)
        if (mesh_.edges[e].boundary)
            y[e] = input[e];
    return y;
}
std::vector<double> SimplicialFluidSolver3D::laplacian(std::span<const double> x) const {
    auto y = curl(mesh_.d1(x));
    std::vector<double> divergence(mesh_.vertices.size(), 0);
    for (size_t e = 0; e < x.size(); ++e) {
        auto v = mesh_.edges[e].vertices;
        const double value = mesh_.star1[e] * x[e];
        divergence[v[0]] -= value;
        divergence[v[1]] += value;
    }
    for (size_t v = 0; v < divergence.size(); ++v)
        divergence[v] = mesh_.boundaryVertices[v] ? 0 : divergence[v] / mesh_.star0[v];
    auto gradient = mesh_.d0(divergence);
    for (size_t e = 0; e < x.size(); ++e)
        y[e] += mesh_.star1[e] * gradient[e];
    return y;
}
LinearSolveResult SimplicialFluidSolver3D::solve(std::span<const double> b, std::vector<double> &x) const {
    if (parameters_.cgIterations < 0 || !std::isfinite(parameters_.cgAbsoluteTolerance) ||
        parameters_.cgAbsoluteTolerance < 0 || !std::isfinite(parameters_.cgRelativeTolerance) ||
        parameters_.cgRelativeTolerance < 0)
        throw std::invalid_argument("Invalid CG parameters");
    auto apply = [&](std::span<const double> p) {
        return applyLaplacian(p);
    };
    std::vector<double> rhs(b.begin(), b.end());
    x.assign(b.size(), 0);
    std::vector<double> r(rhs.begin(), rhs.end()), z(b.size()), p(b.size());
    if(!mesh_.boxDomain) {
        x=phi_;
        auto ax=apply(x);
        for(size_t i=0;i<r.size();++i)r[i]-=ax[i];
    }
    for (size_t e = 0; e < r.size(); ++e) {
        if (!std::isfinite(r[e]))
            throw std::runtime_error("Non-finite CG RHS");
        if (mesh_.edges[e].boundary) {
            r[e] = 0;
            rhs[e] = 0;
        }
    }
    LinearSolveResult result;
    result.initialResidual = result.finalResidual = norm(r);
    const double target =
        std::max(parameters_.cgAbsoluteTolerance, parameters_.cgRelativeTolerance * (mesh_.boxDomain ? result.initialResidual : norm(rhs)));
    auto precondition = [&] {
        if(!incompleteLower_.empty()) {
            z=r;
            for(size_t i=0;i<z.size();++i) {
                for(auto [j,value]:incompleteLower_[i])z[i]-=value*z[j];
                z[i]/=incompleteDiagonal_[i];
            }
            for(size_t i=z.size();i-- >0;) {
                z[i]/=incompleteDiagonal_[i];
                for(auto [j,value]:incompleteLower_[i])z[j]-=value*z[i];
            }
            return;
        }
        for (size_t e = 0; e < z.size(); ++e)
            z[e] = r[e] / diagonal_[e];
    };
    if (result.finalResidual <= target) {
        result.converged = true;
        return result;
    }
    precondition();
    p = z;
    double rz = inner(r, z);
    for (int i = 0; i < parameters_.cgIterations; ++i) {
        auto ap = apply(p);
        double pap = inner(p, ap);
        if (!std::isfinite(pap) || pap <= 0)
            throw std::runtime_error("CG lost positive definiteness");
        double alpha = rz / pap;
        for (size_t e = 0; e < r.size(); ++e) {
            x[e] += alpha * p[e];
            r[e] -= alpha * ap[e];
        }
        result.iterations = i + 1;
        result.finalResidual = norm(r);
        if (result.finalResidual <= target)
            break;
        precondition();
        double next = inner(r, z), beta = next / rz;
        for (size_t e = 0; e < p.size(); ++e)
            p[e] = z[e] + beta * p[e];
        rz = next;
    }
    // Report an actual residual, not only CG's recursively accumulated estimate.
    auto ax = apply(x);
    for (size_t e = 0; e < r.size(); ++e)
        r[e] = (mesh_.edges[e].boundary) ? 0 : rhs[e] - ax[e];
    result.finalResidual = norm(r);
    result.converged = result.finalResidual <= target;
    return result;
}
void SimplicialFluidSolver3D::setPotential(std::span<const double> p) {
    if (p.size() != phi_.size())
        throw std::invalid_argument("Potential size mismatch");
    for (double v : p)
        if (!std::isfinite(v))
            throw std::invalid_argument("Non-finite potential");
    phi_.assign(p.begin(), p.end());
    for (size_t e = 0; e < phi_.size(); ++e)
        if (mesh_.edges[e].boundary)
            phi_[e] = 0;
    flux_ = mesh_.d1(phi_);
    for (size_t f = 0; f < flux_.size(); ++f)
        flux_[f] += harmonic_[f];
    omega_ = curl(flux_);
    std::fill(boundaryCirculation_.begin(), boundaryCirculation_.end(), 0);
    reconstructVelocity(true);
}
void SimplicialFluidSolver3D::recoverFluxFromVorticity() {
    std::vector<double> next;
    recovery_ = compute_ ? compute_->recover(*this,omega_,next) : solve(omega_, next);
    if (!recovery_.converged)
        throw std::runtime_error("3D simplicial flux recovery did not converge; increase solver_iterations");
    phi_ = std::move(next);
    flux_ = mesh_.d1(phi_);
    for (size_t f = 0; f < flux_.size(); ++f)
        flux_[f] += harmonic_[f];
    const auto inside = curl(flux_);
    // Section 4.5: missing boundary circulation = total integral - interior sum.
    for (size_t e = 0; e < omega_.size(); ++e)
        boundaryCirculation_[e] = mesh_.edges[e].boundary ? omega_[e] - inside[e] : 0;
    reconstructVelocity();
}
void SimplicialFluidSolver3D::setBoundaryFlux(std::span<const double> boundary) {
    if (boundary.size() != mesh_.faces.size())
        throw std::invalid_argument("Boundary flux size mismatch");
    std::vector<double> h(boundary.begin(), boundary.end());
    for (size_t f = 0; f < h.size(); ++f) {
        if (!std::isfinite(h[f]) || (!mesh_.faces[f].boundary() && h[f] != 0))
            throw std::invalid_argument("Supply finite fluxes on boundary faces only");
    }
    auto b = mesh_.d2(h);
    double sum = std::accumulate(b.begin(), b.end(), 0.), scale = 0;
    for (double v : h)
        scale += std::abs(v);
    if (std::abs(sum) > 1e-10 * std::max(1., scale))
        throw std::invalid_argument("Boundary flux must have zero net outflow");
    // Cell Neumann Poisson: A=d2 star2^-1 d2^T on interior faces. Its
    // constant nullspace is removed by keeping every iterate zero-mean.
    auto meanZero = [](std::vector<double> &v) {
        double mean = std::accumulate(v.begin(), v.end(), 0.) / v.size();
        for (auto &x : v)
            x -= mean;
    };
    auto apply = [&](const std::vector<double> &p) {
        std::vector<double> y(p.size(), 0);
        for (size_t f = 0; f < h.size(); ++f)
            if (!mesh_.faces[f].boundary()) {
                auto t = mesh_.faces[f].tets;
                double flow = (p[t[0]] - p[t[1]]) / mesh_.star2[f];
                y[t[0]] += flow;
                y[t[1]] -= flow;
            }
        return y;
    };
    for (auto &v : b)
        v = -v;
    meanZero(b);
    std::vector<double> diagonal(b.size(), 0);
    for (size_t f = 0; f < h.size(); ++f)
        if (!mesh_.faces[f].boundary())
            for (int t : mesh_.faces[f].tets)
                diagonal[t] += 1 / mesh_.star2[f];
    std::vector<double> p(b.size(), 0), r = b, z(b.size()), direction;
    for (size_t t = 0; t < z.size(); ++t)
        z[t] = r[t] / diagonal[t];
    meanZero(z);
    direction = z;
    double rr = inner(r, z),
           target = std::max(parameters_.cgAbsoluteTolerance, parameters_.cgRelativeTolerance * norm(b));
    for (int i = 0; i < parameters_.cgIterations && norm(r) > target; ++i) {
        auto ad = apply(direction);
        double denom = inner(direction, ad);
        if (denom <= 0 || !std::isfinite(denom))
            throw std::runtime_error("Harmonic boundary solve breakdown");
        double alpha = rr / denom;
        for (size_t t = 0; t < p.size(); ++t) {
            p[t] += alpha * direction[t];
            r[t] -= alpha * ad[t];
        }
        meanZero(p);
        meanZero(r);
        for (size_t t = 0; t < z.size(); ++t)
            z[t] = r[t] / diagonal[t];
        meanZero(z);
        double next = inner(r, z);
        for (size_t t = 0; t < p.size(); ++t)
            direction[t] = z[t] + (next / rr) * direction[t];
        rr = next;
    }
    auto ap = apply(p);
    for (size_t t = 0; t < r.size(); ++t)
        r[t] = b[t] - ap[t];
    if (norm(r) > target)
        throw std::runtime_error("Harmonic boundary solve did not converge");
    for (size_t f = 0; f < h.size(); ++f)
        if (!mesh_.faces[f].boundary()) {
            auto t = mesh_.faces[f].tets;
            h[f] = (p[t[0]] - p[t[1]]) / mesh_.star2[f];
        }
    for (size_t f = 0; f < h.size(); ++f)
        flux_[f] += h[f] - harmonic_[f];
    harmonic_ = std::move(h);
    reconstructVelocity();
}
// Section 4.4: one flux-consistent vector per tet circumcenter.
void SimplicialFluidSolver3D::reconstructVelocity(bool initializeBoundary) {
    for (size_t i = 0; i < dual_.vertices.size(); ++i) {
        int ti = dual_.vertices[i].tet;
        if (ti < 0)
            continue;
        const auto &t = mesh_.tets[ti];
        Point velocity{};
        for (int k = 0; k < 4; ++k)
            velocity = velocity + (t.center - mesh_.vertices[t.vertices[k]]) *
                                      (t.signs[k] * flux_[t.faces[k]] / (3 * t.volume));
        dualVelocity_[i] = velocity;
    }
    // Initialize the wall trace from an interior vector. On the nonconvex
    // domain, retain the previous trace thereafter: repeatedly replacing its
    // unconstrained tangential component can inject spurious boundary energy.
    // The circulation solve determines the constrained tangential correction.
    for (size_t i = 0; i < dual_.vertices.size(); ++i) {
        const auto &v = dual_.vertices[i];
        if (!v.walls)
            continue;
        if(!mesh_.boxDomain && !initializeBoundary)continue;
        int nearest=boundarySeeds_[i];
        Point velocity = mesh_.boxDomain ? BoundaryCirculation::tangent(dualVelocity_[nearest], v.walls)
                                        : tangentToNormals(dualVelocity_[nearest],v.normals);
        // Nonzero normal boundary flow is supplied by the precomputed harmonic field.
        for (int axis = 0; mesh_.boxDomain && axis < 3; ++axis)
            if (v.walls & (3u << (2 * axis))) {
                double best = 1e300, normalVelocity = 0;
                for (size_t f = 0; f < mesh_.faces.size(); ++f) {
                    const auto &face = mesh_.faces[f];
                    if (!face.boundary())
                        continue;
                    Point n = face.areaNormal;
                    double component = axis == 0 ? n.x : axis == 1 ? n.y : n.z;
                    if (std::abs(component) < .9 * length(n))
                        continue;
                    Point delta = face.center - v.position;
                    double d = dot(delta, delta);
                    if (d < best) {
                        best = d;
                        normalVelocity = harmonic_[f] / component;
                    }
                }
                if (axis == 0)
                    velocity.x = normalVelocity;
                if (axis == 1)
                    velocity.y = normalVelocity;
                if (axis == 2)
                    velocity.z = normalVelocity;
            }
        dualVelocity_[i] = velocity;
    }
    auto inside = curl(flux_);
    if (initializeBoundary) {
        std::vector<double> wall(dual_.edges.size(), 0);
        for (size_t i = 0; i < wall.size(); ++i) {
            auto e = dual_.edges[i];
            if (e.wall)
                wall[i] = .5 * dot(dualVelocity_[e.a] + dualVelocity_[e.b],
                                   dual_.vertices[e.b].position - dual_.vertices[e.a].position);
        }
        auto extra = dual_.circulation(wall);
        for (size_t e = 0; e < omega_.size(); ++e)
            if (mesh_.edges[e].boundary)
                omega_[e] = inside[e] + extra[e];
    }
    // Section 4.5: complementary circulation closes every truncated dual face.
    for (size_t e = 0; e < omega_.size(); ++e)
        boundaryCirculation_[e] = mesh_.edges[e].boundary ? omega_[e] - inside[e] : 0;
    boundarySolve_ = boundary_.reconstruct(dual_, boundaryCirculation_, dualVelocity_, parameters_.cgIterations);
    dual_.prepareInterpolation(dualVelocity_);
}
Point SimplicialFluidSolver3D::sampleVelocity(Point p) const {
    return dual_.interpolate(mesh_, p, dualVelocity_);
}
Point SimplicialFluidSolver3D::trace(Point p, double dt) const {
    if(mesh_.boxDomain) {
        auto clip=[&](Point q){return Point{std::clamp(q.x,0.,mesh_.size),std::clamp(q.y,0.,mesh_.size),std::clamp(q.z,0.,mesh_.size)};};
        Point midpoint=clip(p-sampleVelocity(p)*(dt*.5));
        return clip(p-sampleVelocity(midpoint)*dt);
    }
    Point midpoint = mesh_.clipSegment(p,p - sampleVelocity(p) * (dt * .5));
    return mesh_.clipSegment(p,p - sampleVelocity(midpoint) * dt);
}
// Figure 6, first two loops: trace dual vertices, then integrate dual edges.
void SimplicialFluidSolver3D::advectVorticity(double dt) {
    validDt(dt);
    if (dt == 0)
        return;
    reconstructVelocity();
    if(compute_){compute_->advectCirculation(*this,dt,omega_);return;}
    std::vector<Point> traced(dual_.vertices.size()), velocity(traced.size());
    parallelSamples(int(traced.size()), [&](int i) {
        if(mesh_.boxDomain) {
            traced[i] = trace(dual_.vertices[i].position, dt);
            velocity[i] = sampleVelocity(traced[i]);
        } else {
            Point p=dual_.vertices[i].position;
            Point midpoint=mesh_.clipSegment(p,p-dualVelocity_[i]*(dt*.5),dualTets_[i]);
            Point middleVelocity=dual_.interpolate(mesh_,midpoint,dualVelocity_,dualTets_[i]);
            traced[i]=mesh_.clipSegment(p,p-middleVelocity*dt,dualTets_[i]);
            velocity[i]=dual_.interpolate(mesh_,traced[i],dualVelocity_,dualTets_[i]);
        }
    });
    std::vector<double> integral(dual_.edges.size());
    for (size_t i = 0; i < integral.size(); ++i) {
        const auto &e = dual_.edges[i];
        integral[i] = .5 * dot(velocity[e.a] + velocity[e.b], traced[e.b] - traced[e.a]);
    }
    omega_ = dual_.circulation(integral);
}
void SimplicialFluidSolver3D::addForce(const std::function<Point(Point)> &force, double dt) {
    validDt(dt);
    std::vector<double> f(mesh_.faces.size());
    for (size_t i = 0; i < f.size(); ++i) {
        f[i] = dot(force(mesh_.faces[i].center), mesh_.faces[i].areaNormal);
        if (!std::isfinite(f[i]))
            throw std::invalid_argument("Non-finite body force");
    }
    auto c = curl(f);
    // Complete the truncated dual-face circulation on Euler walls.
    {
        std::vector<double> wall(dual_.edges.size(), 0);
        for (size_t i = 0; i < wall.size(); ++i) {
            const auto e = dual_.edges[i];
            if (!e.wall)
                continue;
            Point a = dual_.vertices[e.a].position, b = dual_.vertices[e.b].position;
            wall[i] = .5 * dot(force(a) + force(b), b - a);
            if (!std::isfinite(wall[i]))
                throw std::invalid_argument("Non-finite wall force");
        }
        auto extra = dual_.circulation(wall);
        for (size_t e = 0; e < c.size(); ++e)
            c[e] += extra[e];
    }
    for (size_t e = 0; e < omega_.size(); ++e)
        omega_[e] += dt * c[e];
}
void SimplicialFluidSolver3D::setDensity(const std::function<double(Point)> &f) {
    auto next = vertexDensity_;
    for (size_t i = 0; i < next.size(); ++i) {
        double value = f(mesh_.vertices[i]);
        if (!std::isfinite(value) || value < 0 || value > 1e30)
            throw std::invalid_argument("Density must be finite and nonnegative");
        next[i] = float(value);
    }
    vertexDensity_ = std::move(next);
    reconstructDensity();
}
void SimplicialFluidSolver3D::reconstructDensity() {
    for (size_t i = 0; i < density_.size(); ++i) {
        const auto &t = mesh_.tets[i];
        double value = 0;
        for (int v : t.vertices)
            value += vertexDensity_[v] * .25;
        density_[i] = float(value); // Exact cell mean of a primal P1 scalar.
    }
}
double SimplicialFluidSolver3D::sampleDensity(Point p) const {
    int t = mesh_.locate(p);
    if (t < 0)
        return 0;
    auto w = mesh_.barycentric(t, p);
    double value = 0, total = 0;
    for (int k = 0; k < 4; ++k) {
        double weight = std::max(0., w[k]);
        value += weight * vertexDensity_[mesh_.tets[t].vertices[k]];
        total += weight;
    }
    return value / total;
}
void SimplicialFluidSolver3D::emitSmoke() {
    // Smoke is a passive visualization scalar, not the constant fluid mass density.
    if (!parameters_.emitterEnabled || parameters_.sourceStrength == 0) return;
    Point emitter{mesh_.size * parameters_.emitterX, mesh_.size * parameters_.emitterY,
                  mesh_.size * parameters_.emitterZ};
    double radius = parameters_.emitterRadius * mesh_.size;
    std::vector<double> added(vertexDensity_.size(), 0);
    timings_.sourceVertices = 0;
    for (size_t i = 0; i < vertexDensity_.size(); ++i)
        if (length(mesh_.vertices[i] - emitter) < radius) {
            ++timings_.sourceVertices;
            added[i] = std::max(0., double(parameters_.sourceStrength - vertexDensity_[i]));
            vertexDensity_[i] = std::max(vertexDensity_[i], parameters_.sourceStrength);
        }
    // Exact mass increment for the piecewise-linear concentration field.
    for (const auto& tet : mesh_.tets)
        for (int v : tet.vertices) timings_.emittedMass += .25 * tet.volume * added[v];
    reconstructDensity();
}
void SimplicialFluidSolver3D::applyBuoyancy(double dt) {
    if(mesh_.boxDomain) {
        addForce([&](Point p) { return Point{0, 0, parameters_.buoyancy * sampleDensity(p)}; }, dt);
        return;
    }
    auto density=[&](const DensitySample& sample) {
        double value=0;for(int k=0;k<4;++k)value+=sample.weights[k]*vertexDensity_[sample.vertices[k]];return value;
    };
    std::vector<double> force(mesh_.faces.size()),dualDensity(dual_.vertices.size());
    for(size_t f=0;f<force.size();++f)force[f]=parameters_.buoyancy*density(faceSamples_[f])*mesh_.faces[f].areaNormal.z;
    auto circulation=curl(force);
    for(size_t i=0;i<dualDensity.size();++i)dualDensity[i]=density(dualSamples_[i]);
    std::vector<double> wall(dual_.edges.size());
    for(size_t i=0;i<wall.size();++i) {
        auto e=dual_.edges[i];if(e.wall)wall[i]=.5*parameters_.buoyancy*(dualDensity[e.a]+dualDensity[e.b])*(dual_.vertices[e.b].position.z-dual_.vertices[e.a].position.z);
    }
    auto extra=dual_.circulation(wall);
    for(size_t e=0;e<omega_.size();++e)omega_[e]+=dt*(circulation[e]+extra[e]);
}
void SimplicialFluidSolver3D::stirEmitter(double dt) {
    if (!parameters_.emitterEnabled || parameters_.sourceStrength <= 0 || parameters_.stirStrength == 0)
        return;
    const Point center{mesh_.size*parameters_.emitterX, mesh_.size*parameters_.emitterY,
                       mesh_.size*parameters_.emitterZ};
    const double radius=mesh_.size*parameters_.stirRadius;
    const double phase=6.283185307179586*parameters_.stirFrequency*(stirringTime_+.5*dt);
    const Point axis=Point{std::cos(phase),std::sin(phase),.5}*(1/std::sqrt(1.25));
    // Smooth compact support. Apply through force curl (including wall terms),
    // then recover flux normally: stirring does not bypass incompressibility.
    addForce([&](Point p) {
        const Point q=(p-center)*(1/radius);
        const double r2=dot(q,q);
        if(r2>=1)return Point{};
        return cross(axis,q)*(parameters_.stirStrength*(1-r2)*(1-r2));
    },dt);
}
void SimplicialFluidSolver3D::substep(double dt) {
    auto previous=std::chrono::steady_clock::now();
    auto elapsed=[&](){auto now=std::chrono::steady_clock::now();double ms=std::chrono::duration<double,std::milli>(now-previous).count();previous=now;return ms;};
    ++timings_.substeps;
    emitSmoke();
    timings_.emissionMs+=elapsed();
    // Figure 6, in order. No energy scaling, circulation correction, or
    // velocity projection is inserted between these inviscid paper operations.
    advectVorticity(dt);
    timings_.advectionMs+=elapsed();
    applyBuoyancy(dt);
    stirEmitter(dt);
    timings_.forcesMs+=elapsed();
    recoverFluxFromVorticity();
    timings_.recoveryMs+=elapsed();
    advectAndDecaySmoke(dt);
    stirringTime_ += dt;
    timings_.smokeMs+=elapsed();
}
void SimplicialFluidSolver3D::advectAndDecaySmoke(double dt) {
    if(compute_){std::vector<float> next;compute_->advectDensity(*this,dt,next);vertexDensity_=std::move(next);reconstructDensity();return;}
    auto next = vertexDensity_;
    double decay=std::exp(-parameters_.smokeDecay*dt);
    parallelSamples(int(next.size()), [&](int i) {
        if(mesh_.boxDomain)next[i]=float(sampleDensity(trace(mesh_.vertices[i],dt))*decay);
        else {
            Point p=mesh_.vertices[i];int hint=primalTets_[i];
            Point middle=mesh_.clipSegment(p,p-dual_.primalVertexVelocity(int(i))*(dt*.5),hint);
            Point velocity=dual_.interpolate(mesh_,middle,dualVelocity_,hint);
            Point q=mesh_.clipSegment(p,p-velocity*dt,hint);
            int ti=mesh_.locate(q,hint);double value=0,total=0;
            if(ti>=0) {
                auto w=mesh_.barycentric(ti,q);
                for(int k=0;k<4;++k){double weight=std::max(0.,w[k]);value+=weight*vertexDensity_[mesh_.tets[ti].vertices[k]];total+=weight;}
            }
            next[i]=total>0?float(value/total*decay):0;
        }
    });
    vertexDensity_ = std::move(next);
    reconstructDensity();
}
void SimplicialFluidSolver3D::advance(float dt) {
    timings_={};
    validDt(dt);
    if (dt == 0)
        return;
    if (!std::isfinite(parameters_.sourceStrength) || parameters_.sourceStrength < 0 ||
        !std::isfinite(parameters_.emitterRadius) || parameters_.emitterRadius <= 0 ||
        !std::isfinite(parameters_.buoyancy) || !std::isfinite(parameters_.smokeDecay) ||
        parameters_.smokeDecay < 0 || !std::isfinite(parameters_.advectionCfl) ||
        parameters_.advectionCfl <= 0 || parameters_.advectionCfl > 1 ||
        !std::isfinite(parameters_.emitterX) || !std::isfinite(parameters_.emitterY) ||
        !std::isfinite(parameters_.emitterZ) ||
        !std::isfinite(parameters_.stirStrength) || parameters_.stirStrength < 0 ||
        !std::isfinite(parameters_.stirRadius) || parameters_.stirRadius <= 0 ||
        !std::isfinite(parameters_.stirFrequency) || parameters_.stirFrequency < 0)
        throw std::invalid_argument("Invalid simplicial fluid parameters");
    if(mesh_.boxDomain){substep(dt);return;}
    double remaining=dt;
    int count=0;
    while(remaining>1e-10) {
        double speed=0;for(auto v:dualVelocity_)speed=std::max(speed,length(v));
        double acceleration=std::abs(parameters_.buoyancy)*parameters_.sourceStrength;
        if(parameters_.emitterEnabled && parameters_.sourceStrength>0)
            acceleration += parameters_.stirStrength;
        double h=parameters_.advectionCfl*mesh_.minAltitude;
        double step=std::min(remaining,2*h/(speed+std::sqrt(speed*speed+4*acceleration*h)+1e-12));
        if(++count>2048)throw std::runtime_error("Bunny CFL subdivision exceeded 2048 steps");
        substep(step);remaining-=step;
    }
}
SimplicialDiagnostics3D SimplicialFluidSolver3D::diagnostics() const {
    SimplicialDiagnostics3D d;
    auto divergence = mesh_.d2(flux_), c = curl(flux_);
    for (size_t t = 0; t < divergence.size(); ++t) {
        d.maxDivergence = std::max(d.maxDivergence, std::abs(divergence[t]) / mesh_.tets[t].volume);
        d.smokeMass += density_[t] * mesh_.tets[t].volume;
    }
    for (size_t f = 0; f < flux_.size(); ++f) {
        d.kineticEnergy += .5 * mesh_.star2[f] * flux_[f] * flux_[f];
        if (mesh_.faces[f].boundary())
            d.maxBoundaryFlux = std::max(d.maxBoundaryFlux, std::abs(flux_[f]));
    }
    for (size_t e = 0; e < c.size(); ++e)
        if (!mesh_.edges[e].boundary)
            d.recoveryCurlError = std::max(d.recoveryCurlError, std::abs(c[e] - omega_[e]));
    for (auto v : dualVelocity_)
        d.maxSpeed = std::max(d.maxSpeed, length(v));
    return d;
}

void SimplicialFluidSolver3D::buildPreconditioner() {
    // IC(0) of the fixed potential operator. A diagonal shift, if needed for
    // factorization, changes only this preconditioner, never the physical solve.
    std::vector<std::map<int,double>> rows(mesh_.edges.size());
    auto add=[&](int i,int j,double value){if(i>j&&!mesh_.edges[i].boundary&&!mesh_.edges[j].boundary)rows[i][j]+=value;};
    for(size_t f=0;f<mesh_.faces.size();++f) {
        const auto& face=mesh_.faces[f];
        for(int a=0;a<3;++a)for(int b=0;b<3;++b)add(face.edges[a],face.edges[b],face.signs[a]*face.signs[b]*mesh_.star2[f]);
    }
    std::vector<std::vector<std::pair<int,double>>> incident(mesh_.vertices.size());
    for(int e=0;e<int(mesh_.edges.size());++e)if(!mesh_.edges[e].boundary) {
        auto v=mesh_.edges[e].vertices;
        incident[v[0]].push_back({e,-mesh_.star1[e]});incident[v[1]].push_back({e,mesh_.star1[e]});
    }
    for(size_t v=0;v<incident.size();++v)if(!mesh_.boundaryVertices[v])
        for(auto [i,a]:incident[v])for(auto [j,b]:incident[v])add(i,j,a*b/mesh_.star0[v]);
    const auto count=rows.size();
    for(double shift:{0.,.001,.01,.1,1.}) {
        incompleteLower_.assign(count,{});incompleteDiagonal_.resize(count);
        bool valid=true;
        for(size_t i=0;i<count&&valid;++i) {
            auto& row=incompleteLower_[i];row.assign(rows[i].begin(),rows[i].end());
            double diagonal=mesh_.edges[i].boundary?1:diagonal_[i]*(1+shift);
            for(size_t k=0;k<row.size();++k) {
                int j=row[k].first;double value=row[k].second;
                const auto& other=incompleteLower_[j];size_t a=0,b=0;
                while(a<k&&b<other.size()) {
                    if(row[a].first==other[b].first){value-=row[a].second*other[b].second;++a;++b;}
                    else if(row[a].first<other[b].first)++a;else ++b;
                }
                row[k].second=value/incompleteDiagonal_[j];diagonal-=row[k].second*row[k].second;
            }
            if(!std::isfinite(diagonal)||diagonal<=1e-14*diagonal_[i])valid=false;
            else incompleteDiagonal_[i]=std::sqrt(diagonal);
        }
        if(valid)return;
    }
    incompleteLower_.clear();incompleteDiagonal_.clear();
}
