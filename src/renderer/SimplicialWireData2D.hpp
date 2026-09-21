#pragma once
#include "renderer/TetrahedralWireData.hpp"
#include "solvers/simplicial2d/SimplicialMesh2D.hpp"

inline TetrahedralWireData simplicialWires2D(const SimplicialMesh2D& mesh, bool outlineOnly = false) {
    TetrahedralWireData data;
    auto point=[](SimplicialPoint2D p){return TetrahedralWireData::Point{p.x,p.y,0};};
    for(const auto& edge:mesh.edges()) {
        auto a=mesh.vertices()[edge.first].position,b=mesh.vertices()[edge.second].position;
        if(!outlineOnly||edge.boundary())data.primal.push_back({point(a),point(b)});
        if(outlineOnly)continue;
        auto first=mesh.triangles()[edge.incidentTriangles[0]].circumcenter;
        auto second=edge.boundary()?.5*(a+b):mesh.triangles()[edge.incidentTriangles[1]].circumcenter;
        data.circumcentric.push_back({point(first),point(second)});
        if(edge.boundary()) {
            data.circumcentric.push_back({point(a),point(second)});
            data.circumcentric.push_back({point(second),point(b)});
        }
    }
    return data;
}

inline void simplicialMvp2D(float* matrix, int width, int height, float size) {
    std::fill(matrix,matrix+16,0.f);
    float fit=float(std::min(width,height));
    float x=fit/std::max(1,width),y=fit/std::max(1,height);
    matrix[0]=2*x/size;matrix[5]=2*y/size;matrix[10]=1;matrix[15]=1;
    matrix[12]=-x;matrix[13]=-y;
}
