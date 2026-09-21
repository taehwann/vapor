#include "renderer/TetrahedralWireData.hpp"

#include <iostream>
#include <sstream>

void require(bool ok,const char* message){if(!ok)throw std::runtime_error(message);}
int main() {
    try {
        const auto bunny=TetrahedralWireData::load("examples/simplicial3d-bunny/bunny.tet",2);
        require(bunny.vertexCount==443&&bunny.tetCount==1996,"Bunny fixture size");
        require(bunny.circumcentric.size()==4*bunny.tetCount,"Missing dual half edges");
        require(bunny.barycentric.size()==bunny.circumcentric.size(),"Dual modes disagree");
        // Euler characteristic of the solid bunny (a topological ball).
        const auto faces=(4*bunny.tetCount+bunny.boundaryFaces)/2;
        require(bunny.vertexCount+faces==bunny.primal.size()+bunny.tetCount+1,"Bunny topology");
        using P=TetrahedralWireData::Point;
        const std::vector<P> vertices{{0,0,0},{1,0,0},{0,1,0},{0,0,1}};
        auto one=TetrahedralWireData::build(vertices,{{0,1,2,3}});
        require(one.primal.size()==6&&one.boundaryFaces==4,"Single tetrahedron topology");
        for(auto s:one.circumcentric)
            require(std::abs(s.a.x-.5)<1e-12&&std::abs(s.a.y-.5)<1e-12&&std::abs(s.a.z-.5)<1e-12,"Circumcenter coordinates");
        for(auto s:one.barycentric)
            require(std::abs(s.a.x-.25)<1e-12&&std::abs(s.a.y-.25)<1e-12&&std::abs(s.a.z-.25)<1e-12,"Barycenter coordinates");
        auto reject=[&](std::vector<std::array<int,4>> cells) {
            bool rejected=false;try{TetrahedralWireData::build(vertices,cells);}catch(const std::invalid_argument&){rejected=true;}
            require(rejected,"Invalid connectivity accepted");
        };
        reject({{0,1,2,4}});reject({{0,1,1,3}});reject({{0,1,2,3},{3,2,1,0}});

        std::cout<<"PASS bunny geometry: "<<bunny.primal.size()<<" edges, "<<bunny.boundaryFaces<<" boundary faces\n";
    }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
