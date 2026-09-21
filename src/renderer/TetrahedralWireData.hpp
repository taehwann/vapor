#pragma once
#include "solvers/simplicial3d/SimplicialMesh3D.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <map>
#include <set>
#include <stdexcept>

// Geometry-only importer: deliberately independent of positive DEC metrics.
// Arbitrary non-Delaunay meshes are useful for inspection but not for simulation.
struct TetrahedralWireData {
    using Point = simplicial3d::Point;
    struct Segment { Point a, b; };
    std::vector<Segment> primal, circumcentric, barycentric;
    std::vector<Point> surfaceTriangles; // Optional depth-only surface, three vertices per triangle.
    size_t vertexCount = 0, tetCount = 0, boundaryFaces = 0;

    static TetrahedralWireData load(const std::filesystem::path& path, double size) {
        std::ifstream in(path);
        std::string magic; int version, nv, nt;
        if (!(in >> magic >> version >> nv >> nt) || magic != "VAPOR_TET" || version != 1 ||
            nv < 4 || nv > 1000000 || nt < 1 || nt > 6000000 || !std::isfinite(size) || size <= 0)
            throw std::invalid_argument("Invalid tetrahedral visualization mesh: " + path.string());
        std::vector<Point> points(nv);
        std::vector<std::array<int,4>> cells(nt);
        for (auto& p : points) {
            if (!(in >> p.x >> p.y >> p.z) || !std::isfinite(p.x) || !std::isfinite(p.y) || !std::isfinite(p.z) ||
                p.x < 0 || p.x > 1 || p.y < 0 || p.y > 1 || p.z < 0 || p.z > 1)
                throw std::invalid_argument("Expected finite unit-box mesh coordinates");
            p = p * size;
        }
        for (auto& c : cells) for (auto& i : c)
            if (!(in >> i)) throw std::invalid_argument("Truncated tetrahedral connectivity");
        if (in >> magic) throw std::invalid_argument("Trailing tetrahedral mesh data");
        return build(points, cells);
    }

    static TetrahedralWireData build(const std::vector<Point>& p, const std::vector<std::array<int,4>>& cells) {
        if (p.empty() || cells.empty()) throw std::invalid_argument("Empty tetrahedral visualization mesh");
        for (auto v : p) if (!std::isfinite(v.x) || !std::isfinite(v.y) || !std::isfinite(v.z))
            throw std::invalid_argument("Non-finite tetrahedral vertex");
        TetrahedralWireData result;
        result.vertexCount = p.size(); result.tetCount = cells.size();
        struct Face { Point circum, bary; int count = 0; double side = 0; };
        std::map<std::array<int,3>, Face> faces;
        std::set<std::array<int,2>> edges;
        std::set<std::array<int,4>> unique;
        auto solve = [](Point a, Point b, Point c, Point rhs) {
            return (cross(b,c)*rhs.x + cross(c,a)*rhs.y + cross(a,b)*rhs.z) * (1 / dot(a,cross(b,c)));
        };
        for (auto ids : cells) {
            std::sort(ids.begin(), ids.end());
            if (ids[0] < 0 || ids[3] >= int(p.size()) || std::adjacent_find(ids.begin(),ids.end()) != ids.end() ||
                !unique.insert(ids).second) throw std::invalid_argument("Invalid or duplicate tetrahedron");
            auto a=p[ids[0]], b=p[ids[1]]-a, c=p[ids[2]]-a, d=p[ids[3]]-a;
            double scale2=std::max({dot(b,b),dot(c,c),dot(d,d)});
            if (!std::isfinite(scale2) || std::abs(dot(b,cross(c,d))) <= 1e-12*std::pow(scale2,1.5))
                throw std::invalid_argument("Degenerate tetrahedron");
            auto circum=a+solve(b,c,d,{dot(b,b)/2,dot(c,c)/2,dot(d,d)/2});
            auto bary=(p[ids[0]]+p[ids[1]]+p[ids[2]]+p[ids[3]])*.25;
            for (int i=0;i<4;++i) {
                for (int j=i+1;j<4;++j) edges.insert({ids[i],ids[j]});
                std::array<int,3> key; int k=0;
                for (int j=0;j<4;++j) if (i!=j) key[k++]=ids[j];
                auto x=p[key[0]], u=p[key[1]]-x, v=p[key[2]]-x, n=cross(u,v);
                auto& face=faces[key];
                double side=dot(n,p[ids[i]]-x);
                if (face.count>=2 || (face.count==1 && side*face.side>=0))
                    throw std::invalid_argument("Non-manifold tetrahedral face");
                if (face.count++==0) {
                    face.circum=x+solve(u,v,n,{dot(u,u)/2,dot(v,v)/2,0});
                    face.bary=(p[key[0]]+p[key[1]]+p[key[2]])*(1./3);
                    face.side=side;
                }
                // Dual half-edges meet at each shared face; boundary halves end on the face.
                result.circumcentric.push_back({circum,face.circum});
                result.barycentric.push_back({bary,face.bary});
            }
        }
        for (auto e:edges) result.primal.push_back({p[e[0]],p[e[1]]});
        for (auto& [key,face]:faces) if (face.count==1) ++result.boundaryFaces;
        return result;
    }
};
