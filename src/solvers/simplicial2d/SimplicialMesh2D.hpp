#pragma once

#include "common/DomainInfo.hpp"

#include <array>
#include <cstddef>
#include <span>
#include <vector>

struct SimplicialPoint2D {
    double x = 0.0;
    double y = 0.0;
};

[[nodiscard]] constexpr SimplicialPoint2D operator+(SimplicialPoint2D a, SimplicialPoint2D b) noexcept {
    return {a.x + b.x, a.y + b.y};
}
[[nodiscard]] constexpr SimplicialPoint2D operator-(SimplicialPoint2D a, SimplicialPoint2D b) noexcept {
    return {a.x - b.x, a.y - b.y};
}
[[nodiscard]] constexpr SimplicialPoint2D operator*(double scale, SimplicialPoint2D value) noexcept {
    return {scale * value.x, scale * value.y};
}
[[nodiscard]] constexpr double dot(SimplicialPoint2D a, SimplicialPoint2D b) noexcept {
    return a.x * b.x + a.y * b.y;
}

// Non-obtuse square or acute teapot triangulation with positive dual weights. Edges have a
// canonical low-to-high vertex orientation and triangles are counter-clockwise.
class SimplicialMesh2D final {
public:
    struct Vertex {
        SimplicialPoint2D position;
        bool boundary = false;
        std::vector<int> incidentTriangles; // Counter-clockwise for interior vertices.
        std::vector<int> boundaryEdges;
    };
    struct Edge {
        int first = -1;
        int second = -1;
        std::array<int, 2> incidentTriangles{-1, -1};
        double primalLength = 0.0;
        double dualLength = 0.0;
        double hodge = 0.0;
        SimplicialPoint2D outwardNormal{}; // Boundary edges only.
        [[nodiscard]] bool boundary() const noexcept { return incidentTriangles[1] < 0; }
    };
    struct Triangle {
        std::array<int, 3> vertices{};
        // Local edges are (v0,v1), (v1,v2), (v2,v0).
        std::array<int, 3> edges{};
        std::array<int, 3> edgeSigns{};
        std::array<int, 3> neighbors{-1, -1, -1};
        SimplicialPoint2D circumcenter;
        double area = 0.0;
    };

    explicit SimplicialMesh2D(int resolution = 32, float boxSize = 4.f);
    [[nodiscard]] static SimplicialMesh2D teapot(int resolution = 80, float size = 4.f);
    [[nodiscard]] bool isTeapot() const noexcept { return teapot_; }

    [[nodiscard]] Dimension dimension() const noexcept { return Dimension::D2; }
    [[nodiscard]] DomainBounds bounds() const noexcept;
    [[nodiscard]] int resolution() const noexcept { return resolution_; }
    [[nodiscard]] float boxSize() const noexcept { return boxSize_; }
    [[nodiscard]] double edgeScale() const noexcept { return double(boxSize_) / resolution_; }
    [[nodiscard]] std::span<const Vertex> vertices() const noexcept { return vertices_; }
    [[nodiscard]] std::span<const Edge> edges() const noexcept { return edges_; }
    [[nodiscard]] std::span<const Triangle> triangles() const noexcept { return triangles_; }
    [[nodiscard]] std::size_t vertexCount() const noexcept { return vertices_.size(); }
    [[nodiscard]] std::size_t edgeCount() const noexcept { return edges_.size(); }
    [[nodiscard]] std::size_t triangleCount() const noexcept { return triangles_.size(); }
    [[nodiscard]] std::size_t vertexIndex(int i, int j) const noexcept;
    [[nodiscard]] SimplicialPoint2D latticeToWorld(double s, double t) const noexcept;
    [[nodiscard]] SimplicialPoint2D clampToDomain(SimplicialPoint2D point) const noexcept;
    // Stops at the first wall, even if the endpoint re-enters across a hole or concavity.
    [[nodiscard]] SimplicialPoint2D clipSegment(SimplicialPoint2D start, SimplicialPoint2D end) const noexcept;
    [[nodiscard]] SimplicialPoint2D wallTangent(int vertex, SimplicialPoint2D velocity) const noexcept;
    [[nodiscard]] int locateTriangle(SimplicialPoint2D point) const noexcept;
    [[nodiscard]] std::array<double, 3> barycentric(int triangle, SimplicialPoint2D point) const noexcept;
    void setBoxSize(float size);

private:
    void buildTopology();
    void buildIncidence();
    void buildBins();
    void updateGeometry();

    int resolution_;
    float boxSize_;
    bool teapot_ = false;
    std::vector<SimplicialPoint2D> unitPositions_;
    std::vector<int> boundaryEdges_;
    std::vector<std::vector<int>> bins_;
    std::vector<std::vector<int>> boundaryBins_;
    std::vector<Vertex> vertices_;
    std::vector<Edge> edges_;
    std::vector<Triangle> triangles_;
};
