#pragma once
#include "SimplicialMesh3D.hpp"

class GpuSimplicial3D;
namespace simplicial3d {
struct DualVertex { Point position; int tet=-1; unsigned walls=0; std::vector<Point> normals; };
Point tangentToNormals(Point p, std::span<const Point> normals);
struct DualEdge { int a,b; bool wall=false; };
struct OrientedDualEdge { int edge,sign; };
struct DualPlane { Point normal; double offset; int primalEdge=-1; };
struct DualCorner {
    int vertex;
    // Triangulation of the normal cone; one triple for a simple vertex.
    std::vector<std::array<int,3>> cones;
};
struct DualCell { std::vector<DualPlane> planes; std::vector<DualCorner> corners; };

// Circumcentric dual. Box interpolation uses Warren et al. coordinates;
// nonconvex domains use a continuous piecewise-linear circumcentric subdivision.
// Geometry and cached interpolation samples; no time integration or stabilization.
class DualMesh {
public:
    explicit DualMesh(const Mesh& mesh);
    std::vector<DualVertex> vertices;
    std::vector<DualEdge> edges;
    std::vector<DualCell> cells; // Box-only convex Voronoi cells.
    std::vector<std::vector<OrientedDualEdge>> loops; // One dual face per primal edge.
    std::vector<int> faceEdges; // Primal face -> normal dual edge.
    std::vector<double> coordinates(int cell,Point p) const;
    Point interpolate(const Mesh& mesh,Point p,std::span<const Point> values,int tetHint=-1) const;
    std::vector<double> circulation(std::span<const double> edgeIntegrals) const;
    // Refresh shared subdivision samples whenever domain dual velocities change.
    // Box interpolation evaluates values directly and needs no preparation.
    void prepareInterpolation(std::span<const Point> values);
    Point primalVertexVelocity(int vertex) const { return samples_.at(vertex); }
private:
    friend class ::GpuSimplicial3D;
    double size_;
    void buildDomain(const Mesh& mesh);
    struct SampleNode { Point position; std::vector<std::pair<int,double>> weights; std::vector<Point> normals; };
    struct Subtet { std::array<int,4> nodes; std::array<Point,3> gradients; };
    std::vector<SampleNode> nodes_;
    std::vector<std::vector<Subtet>> subdivisions_;
    std::vector<Point> samples_;
};
}
