#pragma once
#include <array>
#include <span>
#include <vector>
#include <filesystem>

class GpuSimplicial3D;
namespace simplicial3d {
struct Point {
    double x=0, y=0, z=0;
    Point operator+(Point b) const { return {x+b.x,y+b.y,z+b.z}; }
    Point operator-(Point b) const { return {x-b.x,y-b.y,z-b.z}; }
    Point operator*(double s) const { return {x*s,y*s,z*s}; }
};
inline double dot(Point a, Point b) { return a.x*b.x+a.y*b.y+a.z*b.z; }
inline Point cross(Point a, Point b) { return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x}; }
double length(Point a);
struct Edge { std::array<int,2> vertices; bool boundary=false; };
struct Face {
    std::array<int,3> vertices, edges, signs{1,-1,1};
    std::array<int,2> tets{-1,-1}; // negative/positive side of oriented normal
    Point center, areaNormal;
    bool boundary() const { return tets[0]<0 || tets[1]<0; }
};
struct Tet {
    std::array<int,4> vertices, faces, signs;
    std::array<Point,4> gradients;
    Point center;
    double volume=0;
};
// Paper Sections 2.2 and 3.2: oriented Delaunay tetrahedra and signed
// circumcentric Hodge stars. Invalid metrics are rejected, never regularized.
class Mesh {
public:
    Mesh(std::vector<Point> vertices, std::vector<std::array<int,4>> tets);
    static Mesh loadBox(const std::filesystem::path& file, double size=1);
    static Mesh loadDomain(const std::filesystem::path& file, double size=1);
    Mesh scaledBox(double size) const;
    Point clipSegment(Point start, Point end, int tetHint=-1) const;
    std::vector<Point> vertices;
    std::vector<Edge> edges;
    std::vector<Face> faces;
    std::vector<Tet> tets;
    std::vector<double> star0, star1, star2, star3;
    std::vector<bool> boundaryVertices;
    bool boxDomain=false;
    double size=1, minAltitude=1;
    std::vector<double> d0(std::span<const double>) const;
    std::vector<double> d1(std::span<const double>) const;
    std::vector<double> d2(std::span<const double>) const;
    std::vector<double> d1Transpose(std::span<const double>) const;
    std::array<double,4> barycentric(int tet, Point p) const;
    int locate(Point p, int hint=-1) const;
private:
    friend class ::GpuSimplicial3D;
    void build(std::vector<std::array<int,4>> cells);
    void buildSpatialBins();
    int binsPerAxis_=0;
    std::vector<std::vector<int>> bins_;
    std::vector<std::vector<int>> boundaryBins_;
};
}
