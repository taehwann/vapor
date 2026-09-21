#include "solvers/simplicial2d/SimplicialFluidSolver2D.hpp"
#include "solvers/simplicial2d/SimplicialOperators2D.hpp"
#include "solvers/simplicial2d/SimplicialRasterizer2D.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <queue>
#include <stdexcept>
#include <string>

void require(bool condition,const char* message) {if(!condition)throw std::runtime_error(message);}
void geometry() {
    for(int n:{32,48,80,128}) {
        auto mesh=SimplicialMesh2D::teapot(n,4);
        require(mesh.vertexCount()>100,"Empty teapot");
        require(mesh.vertexCount()+mesh.triangleCount()==mesh.edgeCount(),"Expected one connected teapot with one handle hole");
        std::vector<bool> seen(mesh.triangleCount(),false);std::queue<int> queue;
        queue.push(0);seen[0]=true;
        while(!queue.empty()) {int t=queue.front();queue.pop();for(int next:mesh.triangles()[t].neighbors)
            if(next>=0&&!seen[next]){seen[next]=true;queue.push(next);}}
        require(std::all_of(seen.begin(),seen.end(),[](bool v){return v;}),"Disconnected teapot region");
        for(size_t t=0;t<mesh.triangleCount();++t) {
            auto weights=mesh.barycentric(int(t),mesh.triangles()[t].circumcenter);
            require(*std::min_element(weights.begin(),weights.end())>.32,"Teapot triangles must be acute/equilateral");
            require(mesh.locateTriangle(mesh.triangles()[t].circumcenter)==int(t),"Teapot point location");
        }
        for(auto e:mesh.edges())require(e.hodge>0&&std::isfinite(e.hodge),"Invalid circumcentric metric");
        require(mesh.locateTriangle({.8,1.88})<0,"Handle hole filled with fluid");
        require(mesh.locateTriangle({.1,.1})<0,"Bounding square treated as fluid");
        require(mesh.locateTriangle({1.96,1.5})>=0,"Missing teapot body");
        std::cout<<"mesh "<<n<<": "<<mesh.vertexCount()<<" vertices, "<<mesh.triangleCount()<<" triangles\n";
        auto first=mesh.vertices()[0].position;mesh.setBoxSize(2);
        require(std::abs(mesh.vertices()[0].position.x-.5*first.x)<1e-12,"Resize distorted mesh");
        require(mesh.locateTriangle({.4,.94})<0,"Resize filled handle hole");
    }
}
void walls() {
    auto mesh=SimplicialMesh2D::teapot(80,4);
    // Both endpoints lie inside fluid, with the handle hole between them.
    SimplicialPoint2D a{.32,1.88},b{1.40,1.88};
    require(mesh.locateTriangle(a)>=0&&mesh.locateTriangle(b)>=0,"Invalid crossing fixture");
    auto hit=mesh.clipSegment(a,b);
    require(hit.x<b.x-.5&&mesh.locateTriangle(hit)>=0,"Trace crossed handle hole and re-entered");
    auto outward=mesh.clipSegment({1.96,1.5},{1.96,4});
    require(mesh.locateTriangle(outward)>=0&&outward.y<3.1,"Trace escaped lid");
    std::vector<SimplicialPoint2D> velocity(mesh.triangleCount(),{.7,-.3});
    double slip=0;
    for(auto e:mesh.edges())if(e.boundary()) {
        auto a=mesh.vertices()[e.first].position,b=mesh.vertices()[e.second].position;
        auto v=SimplicialOperators2D::sampleVelocity(mesh,velocity,.5*(a+b));
        require(std::isfinite(v.x)&&std::isfinite(v.y)&&std::abs(dot(v,e.outwardNormal))<1e-9,"Wall normal velocity");
        slip=std::max(slip,std::sqrt(dot(v,v)));
        auto q=mesh.clipSegment(.5*(a+b),.5*(a+b)+.1*e.outwardNormal);
        require(dot(q-.5*(a+b),q-.5*(a+b))<1e-15,"Wall-starting trace escaped");
    }
    require(slip>.1,"Free-slip walls incorrectly force zero tangential velocity");
}
void recovery() {
    SimplicialFluidSolver2D solver{SimplicialFluidState2D(SimplicialMesh2D::teapot())};
    const auto& mesh=solver.domain();
    for(size_t i=0;i<mesh.vertexCount();++i) {
        auto v=mesh.vertices()[i];solver.state().streamFunction()[i]=v.boundary?0:.01*std::sin(3*v.position.x)*std::sin(4*v.position.y);
    }
    solver.synchronizeFromStreamFunction();
    std::vector<double> expected(solver.state().edgeFlux().begin(),solver.state().edgeFlux().end());
    auto old=std::vector<double>(solver.state().vorticity().begin(),solver.state().vorticity().end());
    solver.advectVorticity(0);
    for(size_t i=0;i<old.size();++i)if(!mesh.vertices()[i].boundary)
        require(std::abs(old[i]-solver.state().vorticity()[i])<1e-10,"Zero-step circulation changed");
    std::fill(solver.state().streamFunction().begin(),solver.state().streamFunction().end(),0);
    solver.recoverFluxFromVorticity();require(solver.lastRecovery().converged,"Teapot recovery failed");
    for(size_t i=0;i<expected.size();++i)require(std::abs(expected[i]-solver.state().edgeFlux()[i])<1e-8,"Recovered flux mismatch");
    require(solver.diagnostics().maxDivergence<1e-10,"Teapot divergence");
    std::fill(solver.state().density().begin(),solver.state().density().end(),1.f);
    SimplicialRasterizer2D raster(solver.state());auto data=raster.renderData();
    for(int y=0;y<data.height;++y)for(int x=0;x<data.width;++x) {
        auto p=mesh.latticeToWorld((x+.5)/data.width,(y+.5)/data.height);
        require(data.density[x+data.width*y]==(mesh.locateTriangle(p)>=0?1.f:0.f),"Raster leaked beyond mesh");
    }
}
void stability() {
    SimplicialFluidSolver2D solver{SimplicialFluidState2D(SimplicialMesh2D::teapot())};
    for(int frame=0;frame<600;++frame) {
        solver.advance(1.f/30);
        auto d=solver.diagnostics();
        require(solver.lastRecovery().converged,"Teapot CG failed");
        require(std::isfinite(d.kineticEnergy)&&std::isfinite(d.maxVelocity)&&d.maxVelocity<100,"Teapot unstable");
        require(d.maxDivergence<1e-9,"Teapot lost incompressibility");
        for(size_t i=0;i<solver.domain().edgeCount();++i)if(solver.domain().edges()[i].boundary())
            require(solver.state().edgeFlux()[i]==0,"Teapot wall flux");
        for(float rho:solver.state().density())require(std::isfinite(rho)&&rho>=0&&rho<=1.000001f,"Invalid smoke");
        if((frame+1)%120==0)std::cout<<"t="<<(frame+1)/30.<<" energy="<<d.kineticEnergy<<" speed="<<d.maxVelocity<<" divergence="<<d.maxDivergence<<std::endl;
    }
    require(solver.diagnostics().kineticEnergy>0,"Teapot flow did not start");
}
int main(int argc,char** argv) {
    try {
        std::string mode=argc>1?argv[1]:"geometry";
        if(mode=="geometry")geometry();else if(mode=="walls")walls();else if(mode=="recovery")recovery();else if(mode=="stability")stability();else throw std::invalid_argument("Unknown test");
        std::cout<<"PASS teapot "<<mode<<'\n';return 0;
    }catch(const std::exception& e){std::cerr<<"FAIL "<<e.what()<<'\n';return 1;}
}

