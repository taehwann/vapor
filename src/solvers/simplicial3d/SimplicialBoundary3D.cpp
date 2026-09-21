#include "SimplicialBoundary3D.hpp"
#include <algorithm>
#include <cmath>
#include <map>
#include <numeric>
#include <stdexcept>

namespace simplicial3d {
Point BoundaryCirculation::tangent(Point p,unsigned walls) {
    if(walls&3)p.x=0;if(walls&12)p.y=0;if(walls&48)p.z=0;return p;
}
BoundaryCirculation::BoundaryCirculation(const Mesh& mesh,const DualMesh& dual) {
    for(int e=0;e<int(mesh.edges.size());++e)if(mesh.edges[e].boundary) {
        std::map<int,Point> coefficients;
        for(auto side:dual.loops[e]) {
            const auto& edge=dual.edges[side.edge];if(!edge.wall)continue;
            Point delta=(dual.vertices[edge.b].position-dual.vertices[edge.a].position)*(.5*side.sign);
            for(int v:{edge.a,edge.b})coefficients[v]=coefficients[v]+(mesh.boxDomain ? tangent(delta,dual.vertices[v].walls) : tangentToNormals(delta,dual.vertices[v].normals));
        }
        Row row;row.primalEdge=e;
        for(auto [v,c]:coefficients) {row.coefficients.push_back({v,c});row.diagonal+=dot(c,c);}
        if(row.diagonal<=0)throw std::invalid_argument("Boundary dual face has no tangential circulation degrees of freedom");
        rows_.push_back(std::move(row));
    }
}
LinearSolveResult BoundaryCirculation::reconstruct(const DualMesh& dual,std::span<const double> missing,
                                                 std::vector<Point>& velocity,int maxIterations) const {
    std::vector<double> wallIntegrals(dual.edges.size(),0);
    for(size_t e=0;e<wallIntegrals.size();++e) {
        auto edge=dual.edges[e];if(edge.wall)
            wallIntegrals[e]=.5*dot(velocity[edge.a]+velocity[edge.b],dual.vertices[edge.b].position-dual.vertices[edge.a].position);
    }
    auto current=dual.circulation(wallIntegrals);
    auto transpose=[&](std::span<const double> x) {
        std::vector<Point> y(dual.vertices.size());
        for(size_t r=0;r<rows_.size();++r)for(auto c:rows_[r].coefficients)y[c.vertex]=y[c.vertex]+c.value*x[r];
        return y;
    };
    auto multiply=[&](std::span<const double> x) {
        auto y=transpose(x);std::vector<double> result(rows_.size(),0);
        for(size_t r=0;r<rows_.size();++r)for(auto c:rows_[r].coefficients)result[r]+=dot(c.value,y[c.vertex]);
        return result;
    };
    auto inner=[](const auto& a,const auto& b){return std::inner_product(a.begin(),a.end(),b.begin(),0.);};
    auto norm=[](const auto& a){double value=0;for(double x:a)value=std::max(value,std::abs(x));return value;};
    std::vector<double> rhs(rows_.size()),r,z(rhs.size()),p(rhs.size()),lambda(rhs.size(),0);
    for(size_t i=0;i<rhs.size();++i)rhs[i]=missing[rows_[i].primalEdge]-current[rows_[i].primalEdge];
    r=rhs;LinearSolveResult result;result.initialResidual=result.finalResidual=norm(rhs);
    const double target=std::max(1e-11,result.initialResidual*1e-9);
    for(size_t i=0;i<r.size();++i)p[i]=z[i]=r[i]/rows_[i].diagonal;
    double rz=inner(r,z);
    for(int i=0;i<maxIterations&&result.finalResidual>target;++i) {
        auto ap=multiply(p);double pap=inner(p,ap);
        if(pap<=0||!std::isfinite(pap))throw std::runtime_error("Boundary circulation constraints are inconsistent");
        double alpha=rz/pap;
        for(size_t j=0;j<r.size();++j){lambda[j]+=alpha*p[j];r[j]-=alpha*ap[j];z[j]=r[j]/rows_[j].diagonal;}
        double next=inner(r,z);for(size_t j=0;j<p.size();++j)p[j]=z[j]+next/rz*p[j];rz=next;
        result.iterations=i+1;result.finalResidual=norm(r);
    }
    auto actual=multiply(lambda);for(size_t i=0;i<r.size();++i)r[i]=rhs[i]-actual[i];
    result.finalResidual=norm(r);result.converged=result.finalResidual<=target;
    if(!result.converged)throw std::runtime_error("Boundary circulation reconstruction did not converge");
    auto change=transpose(lambda);for(size_t v=0;v<velocity.size();++v)velocity[v]=velocity[v]+change[v];
    return result;
}
}
