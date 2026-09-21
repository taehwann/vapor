#include "solvers/simplicial3d/SimplicialFluidSolver3D.hpp"
#include "renderer/SimplicialRenderData3D.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
using namespace simplicial3d;
namespace {
Mesh box(double size=2) {return Mesh::loadBox("examples/simplicial3d/box.tet",size);}
void require(bool x,const char* message) {if(!x) throw std::runtime_error(message);}
double norm(std::span<const double> x) {double n=0;for(double a:x)n=std::max(n,std::abs(a));return n;}
double error(std::span<const double> a,std::span<const double> b) {double n=0;for(size_t i=0;i<a.size();++i)n=std::max(n,std::abs(a[i]-b[i]));return n;}
template<class F> void rejects(F f) {try{f();}catch(const std::invalid_argument&){return;}throw std::runtime_error("Invalid input accepted");}
void seed(SimplicialFluidSolver3D& s) {
    const auto& m=s.domain(); std::vector<double> p(m.edges.size());
    for(size_t e=0;e<p.size();++e) {
        auto v=m.edges[e].vertices; Point a=m.vertices[v[0]],b=m.vertices[v[1]],c=(a+b)*.5;
        p[e]=std::sin(3.141592653589793*c.x/m.size)*std::sin(3.141592653589793*c.y/m.size)*std::sin(3.141592653589793*c.z/m.size)*(b.z-a.z);
    }
    s.setPotential(p);
}
void topology() {
    Mesh m=box();
    std::vector<double> v(m.vertices.size()); for(size_t i=0;i<v.size();++i)v[i]=double(i*i%97);
    require(norm(m.d1(m.d0(v)))==0,"d1 d0 != 0");
    std::vector<double> e(m.edges.size());for(size_t i=0;i<e.size();++i)e[i]=double(i*17%101);
    require(norm(m.d2(m.d1(e)))==0,"d2 d1 != 0");
    double volume=0;for(auto t:m.tets)volume+=t.volume;
    require(std::abs(volume-8)<1e-12,"Mesh does not fill box");
    for(size_t i=0;i<m.tets.size();++i)require(m.locate(m.tets[i].center)>=0,"Tet point location");
    require(int(m.vertices.size())-int(m.edges.size())+int(m.faces.size())-int(m.tets.size())==1,"Euler characteristic");
}
void metric() {
    Mesh m=box();
    for(const auto* s:{&m.star0,&m.star1,&m.star2,&m.star3})for(double x:*s)require(x>0&&std::isfinite(x),"Invalid Hodge star");
    require(std::abs(std::accumulate(m.star0.begin(),m.star0.end(),0.)-8)<1e-12,"Dual volume partition");
    // Regular tetrahedron is well-centered: validate actual circumcentric mode.
    Mesh c({{1,1,1},{1,-1,-1},{-1,1,-1},{-1,-1,1}},{{0,1,2,3}});
    require(length(c.tets[0].center)<1e-14,"Tet circumcenter");
    for(double s:c.star2)require(std::abs(s-1./6)<1e-14,"Circumcentric star2");
    require(std::abs(std::accumulate(c.star0.begin(),c.star0.end(),0.)-c.tets[0].volume)<1e-13,"Circumcentric dual volume");
}
void recovery() {
    SimplicialFluidSolver3D s(box(2));seed(s);auto u=s.flux();
    s.recoverFluxFromVorticity();
    require(s.lastRecovery().converged,"CG failed");
    require(error(u,s.flux())<1e-8,"Manufactured flux recovery");
    auto d=s.diagnostics();require(d.maxDivergence<1e-11&&d.maxBoundaryFlux==0,"Recovered flux violates incompressibility/walls");
    require(d.recoveryCurlError<1e-8,"Recovery curl mismatch");
    std::vector<double> a(s.potential().size()),b(a.size());
    for(size_t i=0;i<a.size();++i){a[i]=std::sin(double(i));b[i]=std::cos(double(i));}
    auto la=s.applyLaplacian(a),lb=s.applyLaplacian(b);
    double ab=std::inner_product(a.begin(),a.end(),lb.begin(),0.),ba=std::inner_product(b.begin(),b.end(),la.begin(),0.);
    require(std::abs(ab-ba)<1e-9,"L is not symmetric");
    require(std::inner_product(a.begin(),a.end(),la.begin(),0.)>0,"L is not positive");
    std::cout<<"recovery flux error="<<error(u,s.flux())<<" residual="<<s.lastRecovery().finalResidual<<'\n';
}
void circulation() {
    SimplicialFluidSolver3D s(box(2));seed(s);auto old=s.vorticity();s.advectVorticity(0);
    require(old==s.vorticity(),"Zero-step changed circulation");
    const auto& dual=s.dual();
    auto clip=[](Point p){return Point{std::clamp(p.x,0.,2.),std::clamp(p.y,0.,2.),std::clamp(p.z,0.,2.)};};
    std::vector<Point> traced,velocity;
    for(auto v:dual.vertices) {
        Point midpoint=clip(v.position-s.sampleVelocity(v.position)*.005);
        Point p=clip(v.position-s.sampleVelocity(midpoint)*.01);
        traced.push_back(p);velocity.push_back(s.sampleVelocity(p));
    }
    std::vector<double> expected(s.vorticity().size(),0);
    for(size_t f=0;f<dual.loops.size();++f)for(auto side:dual.loops[f]) {
        const auto e=dual.edges[side.edge];
        expected[f]+=.5*side.sign*dot(velocity[e.a]+velocity[e.b],traced[e.b]-traced[e.a]);
    }
    s.advectVorticity(.01);
    require(error(expected,s.vorticity())<1e-10,"Figure 6 endpoint quadrature mismatch");
    const auto& m=s.domain();const auto& o=s.vorticity();
    // Boundary of a closed dual volume: d0^T Omega=0 at every interior vertex.
    std::vector<double> closed(m.vertices.size(),0);
    for(size_t e=0;e<o.size();++e){closed[m.edges[e].vertices[0]]-=o[e];closed[m.edges[e].vertices[1]]+=o[e];}
    for(size_t v=0;v<closed.size();++v)if(!m.boundaryVertices[v])require(std::abs(closed[v])<1e-12,"Advected dual loops do not cancel");

    s.recoverFluxFromVorticity(); require(s.diagnostics().recoveryCurlError<1e-8,"Advected vorticity not recoverable");
}
void forces() {
    SimplicialFluidSolver3D s(box(2));const auto& m=s.domain();
    auto force=[](Point p){return Point{-p.y,p.x,p.x*p.y};};
    std::vector<double> f(m.faces.size());for(size_t i=0;i<f.size();++i)f[i]=dot(force(m.faces[i].center),m.faces[i].areaNormal);
    auto expected=s.curl(f);s.addForce(force,.125);
    std::vector<double> wall(s.dual().edges.size(),0);
    for(size_t i=0;i<wall.size();++i) {
        auto e=s.dual().edges[i];if(!e.wall)continue;
        auto a=s.dual().vertices[e.a].position,b=s.dual().vertices[e.b].position;
        wall[i]=.5*dot(force(a)+force(b),b-a);
    }
    auto extra=s.dual().circulation(wall);
    for(size_t e=0;e<expected.size();++e)expected[e]=.125*(expected[e]+extra[e]);
    require(error(expected,s.vorticity())<1e-14,"Force is not dt C F");
    s.recoverFluxFromVorticity();require(s.diagnostics().maxDivergence<1e-10,"Force generated divergence");
}
void boundaries() {
    SimplicialFluidSolver3D s(box(2));const auto& m=s.domain();
    seed(s);s.advectVorticity(.001);s.recoverFluxFromVorticity();
    std::vector<double> integral(s.dual().edges.size(),0);
    for(size_t i=0;i<integral.size();++i) {
        auto edge=s.dual().edges[i];if(edge.wall)
            integral[i]=.5*dot(s.dualVelocity()[edge.a]+s.dualVelocity()[edge.b],s.dual().vertices[edge.b].position-s.dual().vertices[edge.a].position);
    }
    require(error(s.dual().circulation(integral),s.boundaryCirculation())<1e-9,"Missing tangential circulation not reconstructed");
    s.resetState();
    std::vector<double> b(m.faces.size());
    for(size_t f=0;f<b.size();++f) if(m.faces[f].boundary()) b[f]=dot(Point{.3,0,0},m.faces[f].areaNormal);
    s.setBoundaryFlux(b);auto h=s.harmonicFlux();
    require(norm(m.d2(h))<1e-9,"Harmonic extension has divergence");
    auto curl=s.curl(h);for(size_t e=0;e<curl.size();++e) if(!m.edges[e].boundary)require(std::abs(curl[e])<1e-10,"Boundary extension is not harmonic");
    for(size_t f=0;f<b.size();++f) if(m.faces[f].boundary())require(h[f]==b[f],"Boundary flux not respected");
    s.recoverFluxFromVorticity(); require(error(h,s.flux())<1e-8,"Recovery lost harmonic boundary field");
    for(size_t f=0;f<b.size();++f)if(m.faces[f].boundary()){b[f]+=1;break;}
    rejects([&]{s.setBoundaryFlux(b);});
}
void simulation() {
    SimplicialFluidSolver3D s(box(2));s.parameters().emitterRadius=.2f;
    for(int i=0;i<60;++i)s.advance(1.f/60);
    auto d=s.diagnostics();
    require(d.kineticEnergy>0&&std::isfinite(d.kineticEnergy)&&d.smokeMass>0,"Plume failed");
    require(d.maxDivergence<1e-9&&d.maxBoundaryFlux==0,"Plume divergence/walls");
    require(d.recoveryCurlError<1e-7,"Plume recovery mismatch");
    SimplicialRenderData3D raster(16);raster.update(s);auto data=raster.renderData();
    require(data.depth==16&&data.density.size()==4096,"Volume shape");
    require(*std::max_element(data.density.begin(),data.density.end())>0,"Empty volume");
    auto flux=s.flux();auto rho=s.density();s.resetState();for(int i=0;i<60;++i)s.advance(1.f/60);
    require(flux==s.flux()&&rho==s.density(),"Reset replay mismatch");
    std::cout<<"plume energy="<<d.kineticEnergy<<" divergence="<<d.maxDivergence<<" smoke mass="<<d.smokeMass<<'\n';
}
void contracts() {
    rejects([]{auto m=box(-1);});
    rejects([]{Mesh m({{0,0,0},{1,0,0},{0,1,0},{1,1,0}},{{0,1,2,3}});});
    SimplicialFluidSolver3D s(box(1));
    rejects([&]{s.advance(-1);});rejects([&]{s.setPotential(std::vector<double>(1));});
    seed(s);s.parameters().cgIterations=0;
    bool failed=false;try{s.recoverFluxFromVorticity();}catch(const std::runtime_error&){failed=true;}
    require(failed,"Failed recovery silently accepted");
    s.setBoxSize(2);require(s.boxSize()==2&&s.diagnostics().kineticEnergy==0,"Resize did not reset");
}
void stability() {
    SimplicialFluidSolver3D s(box(2));s.parameters().emitterRadius=.15f;
    for(int i=0;i<600;++i) {
        s.advance(1.f/60);
        if(i%120==119)std::cout<<"time="<<(i+1)/60.<<" energy="<<s.diagnostics().kineticEnergy<<" speed="<<s.diagnostics().maxSpeed<<" boundaryResidual="<<s.lastBoundaryReconstruction().finalResidual<<std::endl;
    }
    require(s.diagnostics().maxDivergence<1e-8,"Sustained plume divergence");
}
void transport() {
    SimplicialFluidSolver3D s(box(2));s.parameters().sourceStrength=0;s.parameters().buoyancy=0;s.parameters().smokeDecay=0;
    s.setDensity([](Point p){return std::exp(-dot(p-Point{1,1,1},p-Point{1,1,1}));});
    auto rho=s.density();s.advance(.1f);require(rho==s.density(),"Stationary smoke changes without decay");
    seed(s);double initial=s.diagnostics().kineticEnergy;
    for(int i=0;i<30;++i)s.advance(.02f);
    std::cout<<"unforced energy "<<initial<<" -> "<<s.diagnostics().kineticEnergy<<std::endl;
    require(std::isfinite(s.diagnostics().kineticEnergy)&&s.diagnostics().kineticEnergy<initial*2,"Unforced energy diverged");
    require(s.diagnostics().maxDivergence<1e-10,"Unforced flow broke incompressibility");
}
}
int main(int argc,char** argv) {
    try {
        std::string c=argc>1?argv[1]:"";
        if(c=="topology")topology();else if(c=="metric")metric();else if(c=="recovery")recovery();
        else if(c=="circulation")circulation();else if(c=="forces")forces();
        else if(c=="simulation")simulation();else if(c=="contracts")contracts();else if(c=="boundaries")boundaries();else if(c=="stability")stability();else if(c=="transport")transport();else throw std::runtime_error("Unknown case");
        std::cout<<"PASS "<<c<<'\n';return 0;
    }catch(const std::exception& e){std::cerr<<"FAIL "<<e.what()<<'\n';return 1;}
}
