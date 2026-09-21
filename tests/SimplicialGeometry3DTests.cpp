#include "solvers/simplicial3d/SimplicialDual3D.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <stdexcept>
using namespace simplicial3d;
void require(bool value,const char* why){if(!value)throw std::runtime_error(why);}
int main() {
    try {
        auto mesh=Mesh::loadBox("examples/simplicial3d/box.tet",2);
        DualMesh dual(mesh);
        for(const auto& t:mesh.tets) {
            Point d=mesh.vertices[t.vertices[0]]-t.center;double radius2=dot(d,d);
            for(auto p:mesh.vertices){Point q=p-t.center;require(dot(q,q)>=radius2-1e-10,"Non-Delaunay tetrahedron");}
        }
        double worst=0,minWeight=1;
        for(int i=0;i<int(dual.cells.size());++i) {
            Point center{};for(auto& c:dual.cells[i].corners)center=center+dual.vertices[c.vertex].position;
            center=center*(1./dual.cells[i].corners.size());
            for(auto& corner:dual.cells[i].corners) {
                Point query=center*.37+dual.vertices[corner.vertex].position*.63;
                auto w=dual.coordinates(i,query);Point reproduced{};
                for(size_t j=0;j<w.size();++j){reproduced=reproduced+dual.vertices[dual.cells[i].corners[j].vertex].position*w[j];minWeight=std::min(minWeight,w[j]);}
                worst=std::max(worst,length(query-reproduced));
                require(std::abs(std::accumulate(w.begin(),w.end(),0.)-1)<1e-12,"Partition of unity");
            }
        }
        require(worst<1e-8,"Generalized barycentric coordinates lack linear precision");
        require(minWeight>=0,"Negative generalized coordinate");
        std::vector<Point> arbitrary;for(size_t i=0;i<dual.vertices.size();++i)arbitrary.push_back({std::sin(double(i)),std::cos(double(i)),std::sin(.17*i)});
        for(size_t e=0;e<mesh.edges.size();++e) {
            Point area{};
            for(auto side:dual.loops[e]) {
                auto segment=dual.edges[side.edge];
                area=area+cross(dual.vertices[segment.a].position,dual.vertices[segment.b].position)*(.5*side.sign);
            }
            auto sites=mesh.edges[e].vertices;Point direction=mesh.vertices[sites[1]]-mesh.vertices[sites[0]];
            require(std::abs(dot(area,direction)/dot(direction,direction)-mesh.star1[e])<1e-10,"Hodge star1 differs from clipped Voronoi face measure");
            Point p{};for(auto side:dual.loops[e])p=p+dual.vertices[dual.edges[side.edge].a].position+dual.vertices[dual.edges[side.edge].b].position;
            p=p*(.5/dual.loops[e].size());
            Point values[2];for(int k=0;k<2;++k) {
                int site=mesh.edges[e].vertices[k];auto w=dual.coordinates(site,p);
                for(size_t j=0;j<w.size();++j)values[k]=values[k]+arbitrary[dual.cells[site].corners[j].vertex]*w[j];
            }
            require(length(values[0]-values[1])<1e-7,"Discontinuous interpolation across Voronoi face");
        }
        std::vector<Point> affine;for(auto v:dual.vertices)affine.push_back({1+2*v.position.x-v.position.y,3*v.position.y-v.position.z,1-v.position.x});
        for(auto p:mesh.vertices) {
            Point exact{1+2*p.x-p.y,3*p.y-p.z,1-p.x};
            require(length(dual.interpolate(mesh,p,affine)-exact)<1e-8,"Affine velocity interpolation, including walls");
        }
        std::vector<double> potential(mesh.vertices.size());for(size_t i=0;i<potential.size();++i)potential[i]=double(i%97);
        auto exact=mesh.d1(mesh.d0(potential));for(double x:exact)require(x==0,"d1 d0");
        std::cout<<"PASS mesh/dual: "<<mesh.vertices.size()<<" vertices, "<<mesh.tets.size()<<" tets, "<<dual.vertices.size()<<" dual vertices, linear error="<<worst<<'\n';
    }catch(const std::exception& e){std::cerr<<"FAIL "<<e.what()<<'\n';return 1;}
}
