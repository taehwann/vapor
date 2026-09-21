#pragma once
#include "TetrahedralWireData.hpp"
#include "solvers/simplicial3d/SimplicialDual3D.hpp"

inline TetrahedralWireData simplicialWires3D(const simplicial3d::Mesh& mesh,
                                            const simplicial3d::DualMesh& dual,
                                            bool boundaryOnly = false) {
    TetrahedralWireData data;
    if (boundaryOnly)
        for (const auto& face : mesh.faces)
            if (face.boundary())
                for (int v : face.vertices) data.surfaceTriangles.push_back(mesh.vertices[v]);
    for (const auto& e : mesh.edges)
        if (!boundaryOnly || e.boundary)
            data.primal.push_back({mesh.vertices[e.vertices[0]], mesh.vertices[e.vertices[1]]});
    for (const auto& e : dual.edges)
        if (!boundaryOnly || e.wall)
            data.circumcentric.push_back({dual.vertices[e.a].position, dual.vertices[e.b].position});
    return data;
}
