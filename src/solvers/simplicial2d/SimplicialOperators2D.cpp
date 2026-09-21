#include "solvers/simplicial2d/SimplicialOperators2D.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

namespace {
void requireSize(std::size_t actual, std::size_t expected, const char* field) {
    if (actual != expected) throw std::invalid_argument(std::string("Invalid simplicial ") + field + " extent");
}

void applyLaplacian(const SimplicialMesh2D& mesh, std::span<const double> input,
                    std::span<double> output) {
    std::fill(output.begin(), output.end(), 0.0);
    const auto vertices = mesh.vertices();
    for (const auto& edge : mesh.edges()) {
        const bool firstInterior = !vertices[edge.first].boundary;
        const bool secondInterior = !vertices[edge.second].boundary;
        if (firstInterior && secondInterior) {
            const double difference = input[edge.first] - input[edge.second];
            output[edge.first] += edge.hodge * difference;
            output[edge.second] -= edge.hodge * difference;
        } else if (firstInterior) {
            output[edge.first] += edge.hodge * input[edge.first];
        } else if (secondInterior) {
            output[edge.second] += edge.hodge * input[edge.second];
        }
    }
    for (std::size_t i = 0; i < vertices.size(); ++i)
        if (vertices[i].boundary) output[i] = input[i];
}
}

void SimplicialOperators2D::fluxFromStreamFunction(
    const SimplicialMesh2D& mesh, std::span<const double> streamFunction,
    std::span<double> flux) {
    requireSize(streamFunction.size(), mesh.vertexCount(), "stream function");
    requireSize(flux.size(), mesh.edgeCount(), "flux");
    for (std::size_t i = 0; i < mesh.edgeCount(); ++i) {
        const auto& edge = mesh.edges()[i];
        flux[i] = streamFunction[edge.second] - streamFunction[edge.first];
    }
}

void SimplicialOperators2D::vorticityFromFlux(
    const SimplicialMesh2D& mesh, std::span<const double> flux,
    std::span<double> vorticity) {
    requireSize(flux.size(), mesh.edgeCount(), "flux");
    requireSize(vorticity.size(), mesh.vertexCount(), "vorticity");
    std::fill(vorticity.begin(), vorticity.end(), 0.0);
    for (std::size_t i = 0; i < mesh.edgeCount(); ++i) {
        const auto& edge = mesh.edges()[i];
        const double circulation = edge.hodge * flux[i];
        vorticity[edge.first] -= circulation;
        vorticity[edge.second] += circulation;
    }
}

std::vector<SimplicialPoint2D> SimplicialOperators2D::reconstructTriangleVelocities(
    const SimplicialMesh2D& mesh, std::span<const double> flux) {
    requireSize(flux.size(), mesh.edgeCount(), "flux");
    std::vector<SimplicialPoint2D> result(mesh.triangleCount());
    for (std::size_t ti = 0; ti < mesh.triangleCount(); ++ti) {
        const auto& triangle = mesh.triangles()[ti];
        double a00 = 0.0, a01 = 0.0, a11 = 0.0, b0 = 0.0, b1 = 0.0;
        for (const int edgeIndex : triangle.edges) {
            const auto& edge = mesh.edges()[edgeIndex];
            const auto tangent = mesh.vertices()[edge.second].position - mesh.vertices()[edge.first].position;
            const double nx = -tangent.y, ny = tangent.x;
            a00 += nx * nx;
            a01 += nx * ny;
            a11 += ny * ny;
            b0 += nx * flux[edgeIndex];
            b1 += ny * flux[edgeIndex];
        }
        const double determinant = a00 * a11 - a01 * a01;
        if (!(determinant > std::numeric_limits<double>::epsilon() * (a00 + a11) * (a00 + a11)))
            throw std::runtime_error("Could not reconstruct velocity on a degenerate triangle");
        result[ti] = {(a11 * b0 - a01 * b1) / determinant,
                      (a00 * b1 - a01 * b0) / determinant};
    }
    return result;
}

SimplicialPoint2D SimplicialOperators2D::sampleVelocity(
    const SimplicialMesh2D& mesh, std::span<const SimplicialPoint2D> triangleVelocity,
    SimplicialPoint2D point) noexcept {
    point = mesh.clampToDomain(point);
    const int index = mesh.locateTriangle(point);
    if (index < 0) return {};
    const auto& triangle = mesh.triangles()[index];
    // The nearest primal vertex owns the containing Voronoi cell.
    int nearest = triangle.vertices[0];
    for (int candidate : triangle.vertices) {
        const auto a = mesh.vertices()[candidate].position - point;
        const auto b = mesh.vertices()[nearest].position - point;
        if (dot(a, a) < dot(b, b)) nearest = candidate;
    }
    // Interpolate on the actual (possibly wall-clipped) circumcentric dual.
    // The old boundary fan assumed circumcenters were triangle centroids,
    // which is only true for equilateral triangles.
    struct Node { SimplicialPoint2D position, velocity; };
    std::array<Node, 16> nodes{};
    std::size_t count = 0;
    const double epsilon = 1e-12 * mesh.edgeScale();
    const auto wallVelocity = [&](SimplicialPoint2D p, SimplicialPoint2D v) {
        bool constrained=false;
        SimplicialPoint2D normal{};
        for(int ei:mesh.vertices()[nearest].boundaryEdges) {
            const auto& e=mesh.edges()[ei];
            auto a=mesh.vertices()[e.first].position,d=mesh.vertices()[e.second].position-a;
            double t=dot(p-a,d)/dot(d,d);
            if(t < -1e-10 || t > 1+1e-10 || std::abs(dot(p-a,e.outwardNormal))>epsilon)continue;
            if(constrained && std::abs(normal.x*e.outwardNormal.y-normal.y*e.outwardNormal.x)>1e-8)return SimplicialPoint2D{};
            normal=e.outwardNormal;constrained=true;
        }
        return constrained?v-dot(v,normal)*normal:v;
    };
    const auto append = [&](SimplicialPoint2D p, SimplicialPoint2D v) {
        for (std::size_t i = 0; i < count; ++i) {
            const auto delta = nodes[i].position - p;
            if (dot(delta, delta) < epsilon * epsilon) return;
        }
        nodes[count++] = {p, wallVelocity(p, v)};
    };
    const auto& vertex = mesh.vertices()[nearest];
    if (!vertex.boundary) {
        for (int ti:vertex.incidentTriangles) append(mesh.triangles()[ti].circumcenter,triangleVelocity[ti]);
    } else {
        // Walk the actual triangle fan from one wall edge to the other. Sorting
        // around an average point is incorrect at a nonconvex boundary corner.
        int incoming=vertex.boundaryEdges.front();
        int ti=mesh.edges()[incoming].incidentTriangles[0];
        auto midpoint=[&](int ei) {const auto& e=mesh.edges()[ei];return .5*(mesh.vertices()[e.first].position+mesh.vertices()[e.second].position);};
        append(midpoint(incoming),triangleVelocity[ti]);
        SimplicialPoint2D average{};
        for(size_t k=0;k<vertex.incidentTriangles.size();++k) {
            append(mesh.triangles()[ti].circumcenter,triangleVelocity[ti]);
            average=average+triangleVelocity[ti];
            int outgoing=-1;
            for(int ei:mesh.triangles()[ti].edges) {
                const auto& e=mesh.edges()[ei];
                if(ei!=incoming&&(e.first==nearest||e.second==nearest)){outgoing=ei;break;}
            }
            const auto& edge=mesh.edges()[outgoing];
            if(edge.boundary()){append(midpoint(outgoing),triangleVelocity[ti]);break;}
            ti=edge.incidentTriangles[0]==ti?edge.incidentTriangles[1]:edge.incidentTriangles[0];
            incoming=outgoing;
        }
        append(vertex.position,(1./vertex.incidentTriangles.size())*average);
        double area=0;
        for(size_t i=0;i<count;++i){auto a=nodes[i].position,b=nodes[(i+1)%count].position;area+=a.x*b.y-a.y*b.x;}
        if(area<0)std::reverse(nodes.begin(),nodes.begin()+count);
    }
    std::array<SimplicialPoint2D, 16> offset;
    std::array<double, 16> distance, halfAngle;
    for (std::size_t i = 0; i < count; ++i) {
        offset[i] = nodes[i].position - point;
        distance[i] = std::sqrt(dot(offset[i], offset[i]));
        if (distance[i] < epsilon) return nodes[i].velocity;
    }
    for (std::size_t i = 0; i < count; ++i) {
        const auto j = (i + 1) % count;
        const double cross = offset[i].x * offset[j].y - offset[i].y * offset[j].x;
        const double product = dot(offset[i], offset[j]);
        if (std::abs(cross) < epsilon * mesh.edgeScale() && product < 0.0)
            return (1.0 / (distance[i] + distance[j])) *
                (distance[j] * nodes[i].velocity + distance[i] * nodes[j].velocity);
        // Rationalize near a dual edge: adding nearly opposite vectors can
        // round the usual denominator to zero even for a finite query point.
        halfAngle[i] = product < 0.0
            ? (distance[i] * distance[j] - product) / cross
            : cross / (distance[i] * distance[j] + product);
    }
    double sum = 0.0;
    SimplicialPoint2D velocity{};
    for (std::size_t i = 0; i < count; ++i) {
        const double weight = (halfAngle[(i + count - 1) % count] + halfAngle[i]) / distance[i];
        velocity = velocity + weight * nodes[i].velocity;
        sum += weight;
    }
    return (1.0 / sum) * velocity;
}

float SimplicialOperators2D::sampleVertexScalar(
    const SimplicialMesh2D& mesh, std::span<const float> values,
    SimplicialPoint2D point) noexcept {
    const int index = mesh.locateTriangle(point);
    if (index < 0) return 0.f;
    auto coordinates = mesh.barycentric(index, point);
    for (double& coordinate : coordinates) coordinate = std::clamp(coordinate, 0.0, 1.0);
    const double sum = coordinates[0] + coordinates[1] + coordinates[2];
    const auto& triangle = mesh.triangles()[index];
    return static_cast<float>((coordinates[0] * values[triangle.vertices[0]] +
                               coordinates[1] * values[triangle.vertices[1]] +
                               coordinates[2] * values[triangle.vertices[2]]) / sum);
}

double SimplicialOperators2D::maxDivergence(
    const SimplicialMesh2D& mesh, std::span<const double> flux) {
    requireSize(flux.size(), mesh.edgeCount(), "flux");
    double result = 0.0;
    for (const auto& triangle : mesh.triangles()) {
        double netFlux = 0.0;
        for (int local = 0; local < 3; ++local)
            netFlux += triangle.edgeSigns[local] * flux[triangle.edges[local]];
        result = std::max(result, std::abs(netFlux / triangle.area));
    }
    return result;
}

LinearSolveResult SimplicialOperators2D::recoverFlux(
    const SimplicialMesh2D& mesh, std::span<const double> vorticity,
    std::span<double> streamFunction, std::span<double> flux,
    const LinearSolveOptions& options) {
    requireSize(vorticity.size(), mesh.vertexCount(), "vorticity");
    requireSize(streamFunction.size(), mesh.vertexCount(), "stream function");
    requireSize(flux.size(), mesh.edgeCount(), "flux");
    if (options.maxIterations < 0 || !std::isfinite(options.absoluteTolerance) ||
        options.absoluteTolerance < 0.0 || !std::isfinite(options.relativeTolerance) ||
        options.relativeTolerance < 0.0)
        throw std::invalid_argument("Invalid simplicial CG options");

    const std::size_t count = mesh.vertexCount();
    std::vector<double> residual(count), direction(count), applied(count), inverseDiagonal(count), z(count);
    for (const auto& edge : mesh.edges()) {
        if (!mesh.vertices()[edge.first].boundary) inverseDiagonal[edge.first] += edge.hodge;
        if (!mesh.vertices()[edge.second].boundary) inverseDiagonal[edge.second] += edge.hodge;
    }
    for (std::size_t i = 0; i < count; ++i) {
        if (mesh.vertices()[i].boundary) {
            streamFunction[i] = 0.0;
            inverseDiagonal[i] = 1.0;
        } else {
            inverseDiagonal[i] = 1.0 / inverseDiagonal[i];
        }
    }

    applyLaplacian(mesh, streamFunction, applied);
    double residualSquared = 0.0;
    for (std::size_t i = 0; i < count; ++i) {
        residual[i] = mesh.vertices()[i].boundary ? 0.0 : vorticity[i] - applied[i];
        z[i] = inverseDiagonal[i] * residual[i];
        direction[i] = z[i];
        residualSquared += residual[i] * residual[i];
    }
    LinearSolveResult result;
    result.initialResidual = std::sqrt(residualSquared);
    result.finalResidual = result.initialResidual;
    const double tolerance = std::max(options.absoluteTolerance,
                                      options.relativeTolerance * result.initialResidual);
    if (result.initialResidual <= tolerance) {
        result.converged = true;
        fluxFromStreamFunction(mesh, streamFunction, flux);
        return result;
    }

    double rz = 0.0;
    for (std::size_t i = 0; i < count; ++i) rz += residual[i] * z[i];
    for (int iteration = 0; iteration < options.maxIterations; ++iteration) {
        applyLaplacian(mesh, direction, applied);
        double denominator = 0.0;
        for (std::size_t i = 0; i < count; ++i) denominator += direction[i] * applied[i];
        if (!(denominator > 0.0) || !std::isfinite(denominator)) break;
        const double alpha = rz / denominator;
        residualSquared = 0.0;
        for (std::size_t i = 0; i < count; ++i) {
            streamFunction[i] += alpha * direction[i];
            residual[i] -= alpha * applied[i];
            residualSquared += residual[i] * residual[i];
        }
        result.iterations = iteration + 1;
        result.finalResidual = std::sqrt(residualSquared);
        if (!options.fixedIterations && result.finalResidual <= tolerance) {
            result.converged = true;
            break;
        }
        double nextRz = 0.0;
        for (std::size_t i = 0; i < count; ++i) {
            z[i] = inverseDiagonal[i] * residual[i];
            nextRz += residual[i] * z[i];
        }
        if (!(rz > 0.0) || !std::isfinite(nextRz)) break;
        const double beta = nextRz / rz;
        for (std::size_t i = 0; i < count; ++i) direction[i] = z[i] + beta * direction[i];
        rz = nextRz;
    }
    if (options.fixedIterations) result.converged = result.finalResidual <= tolerance;
    fluxFromStreamFunction(mesh, streamFunction, flux);
    return result;
}
