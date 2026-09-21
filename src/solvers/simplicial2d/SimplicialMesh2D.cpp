#include "solvers/simplicial2d/SimplicialMesh2D.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <unordered_map>

namespace {

double cross(SimplicialPoint2D a, SimplicialPoint2D b) noexcept {
    return a.x * b.y - a.y * b.x;
}

} // namespace

SimplicialMesh2D SimplicialMesh2D::teapot(int resolution, float size) {
    if (resolution < 32 || resolution > 256)
        throw std::invalid_argument("Teapot resolution must be 32..256");
    SimplicialMesh2D mesh(2,size);
    mesh.resolution_=resolution; mesh.teapot_=true;
    mesh.vertices_.clear(); mesh.triangles_.clear(); mesh.unitPositions_.clear();
    // A stylized silhouette, not a projection of the Utah teapot asset. Keep whole
    // equilateral elements: all circumcenters remain inside their triangles and
    // all circumcentric weights remain strictly positive, including at the walls.
    auto inside=[](SimplicialPoint2D p) {
        auto ellipse=[&](double x,double y,double rx,double ry) {
            double a=(p.x-x)/rx,b=(p.y-y)/ry;return a*a+b*b<=1;
        };
        const std::array<SimplicialPoint2D,8> spout{{{.67,.31},{.80,.38},{.92,.60},{.97,.68},
                                                   {.89,.68},{.82,.56},{.74,.52},{.67,.49}}};
        bool inSpout=false;
        for(size_t i=0,j=spout.size()-1;i<spout.size();j=i++) {
            auto a=spout[i],b=spout[j];
            if((a.y>p.y)!=(b.y>p.y) && p.x<(b.x-a.x)*(p.y-a.y)/(b.y-a.y)+a.x) inSpout=!inSpout;
        }
        bool body=ellipse(.49,.43,.255,.23)&&p.y<=.64;
        bool foot=p.x>=.31&&p.x<=.67&&p.y>=.195&&p.y<=.235;
        bool lid=ellipse(.49,.65,.20,.052)||ellipse(.49,.715,.045,.043);
        bool handle=ellipse(.22,.47,.17,.18);
        bool hole=ellipse(.20,.47,.09,.095);
        return (body||foot||lid||handle||inSpout)&&!hole;
    };
    const int n=resolution,rows=int(n*2/std::sqrt(3.));
    std::vector<SimplicialPoint2D> points;
    for(int j=0;j<=rows;++j)for(int i=0;i<=n;++i)
        points.push_back({(i+.5*(j%2))/n,j*std::sqrt(3.)/(2*n)});
    std::vector<std::array<int,3>> cells;
    auto add=[&](std::array<int,3> c) {
        bool keep=true;
        for(int i=0;i<3;++i) {
            keep=keep&&inside(points[c[i]])&&inside(.5*(points[c[i]]+points[c[(i+1)%3]]));
        }
        if(keep&&inside((1./3)*(points[c[0]]+points[c[1]]+points[c[2]])))cells.push_back(c);
    };
    for(int j=0;j<rows;++j)for(int i=0;i<n;++i) {
        int a=i+(n+1)*j,b=a+1,c=a+n+1,d=c+1;
        if(j%2==0){add({a,b,c});add({b,d,c});}
        else{add({a,b,d});add({a,d,c});}
    }
    // Remove point-only contacts introduced by sampling the outline. A wall
    // vertex must have one connected fan, not two fluid regions touching at a point.
    for(int pass=0;pass<16;++pass) {
        std::vector<std::vector<int>> stars(points.size());
        for(int t=0;t<int(cells.size());++t)for(int v:cells[t])stars[v].push_back(t);
        std::vector<bool> remove(cells.size(),false);
        bool changed=false;
        for(const auto& star:stars) {
            std::vector<std::vector<int>> groups;
            std::vector<int> remaining=star;
            while(!remaining.empty()) {
                std::vector<int> group{remaining.back()};remaining.pop_back();
                for(size_t k=0;k<group.size();++k)for(size_t j=0;j<remaining.size();) {
                    int shared=0;for(int a:cells[group[k]])for(int b:cells[remaining[j]])shared+=a==b;
                    if(shared==2){group.push_back(remaining[j]);remaining.erase(remaining.begin()+j);}else ++j;
                }
                groups.push_back(std::move(group));
            }
            if(groups.size()>1) {
                const auto best=std::max_element(groups.begin(),groups.end(),[](auto& a,auto& b){return a.size()<b.size();});
                for(auto it=groups.begin();it!=groups.end();++it)if(it!=best)for(int t:*it)remove[t]=true;
                changed=true;
            }
        }
        if(!changed)break;
        std::vector<std::array<int,3>> kept;
        for(size_t t=0;t<cells.size();++t)if(!remove[t])kept.push_back(cells[t]);
        cells=std::move(kept);
    }
    std::vector<int> remap(points.size(),-1);
    for(auto c:cells) {
        for(auto& v:c) {
            if(remap[v]<0){remap[v]=int(mesh.unitPositions_.size());mesh.unitPositions_.push_back(points[v]);}
            v=remap[v];
        }
        mesh.triangles_.push_back({c});
    }
    mesh.vertices_.resize(mesh.unitPositions_.size());
    mesh.buildIncidence();
    for(const auto& v:mesh.vertices_)if(v.boundary&&v.boundaryEdges.size()!=2)
        throw std::runtime_error("Teapot wall is not a closed manifold curve");
    mesh.updateGeometry();
    return mesh;
}

void SimplicialMesh2D::buildBins() {
    const int n=resolution_;
    bins_.assign(n*n,{});
    boundaryBins_.assign(n*n,{});
    auto bin=[&](double x){return std::clamp(int(x/boxSize_*n),0,n-1);};
    for(int ti=0;ti<int(triangles_.size());++ti) {
        SimplicialPoint2D lo{boxSize_,boxSize_},hi{};
        for(int v:triangles_[ti].vertices) {
            auto p=vertices_[v].position;lo={std::min(lo.x,p.x),std::min(lo.y,p.y)};
            hi={std::max(hi.x,p.x),std::max(hi.y,p.y)};
        }
        for(int y=bin(lo.y);y<=bin(hi.y);++y)for(int x=bin(lo.x);x<=bin(hi.x);++x)bins_[x+n*y].push_back(ti);
    }
    for(int ei:boundaryEdges_) {
        auto a=vertices_[edges_[ei].first].position,b=vertices_[edges_[ei].second].position;
        for(int y=bin(std::min(a.y,b.y));y<=bin(std::max(a.y,b.y));++y)
            for(int x=bin(std::min(a.x,b.x));x<=bin(std::max(a.x,b.x));++x)boundaryBins_[x+n*y].push_back(ei);
    }
}

SimplicialPoint2D SimplicialMesh2D::clipSegment(SimplicialPoint2D start,SimplicialPoint2D end) const noexcept {
    if(!teapot_)return clampToDomain(end);
    const auto r=end-start;
    double first=1;
    auto intersect=[&](int ei) {
        const auto& e=edges_[ei];
        if(dot(r,e.outwardNormal)<=1e-14*boxSize_)return;
        const auto a=vertices_[e.first].position,s=vertices_[e.second].position-a;
        const double denominator=cross(r,s);
        if(std::abs(denominator)<1e-16*boxSize_*boxSize_)return;
        const double t=cross(a-start,s)/denominator,u=cross(a-start,r)/denominator;
        if(t>=-1e-10&&t<=first&&u>=-1e-10&&u<=1+1e-10)first=std::max(0.,t);
    };
    auto bin=[&](double x){return std::clamp(int(x/boxSize_*resolution_),0,resolution_-1);};
    // CFL-sized paths usually cross one or two bins, avoiding a scan of the entire outline.
    const double eps=1e-10*boxSize_;
    int x0=bin(std::min(start.x,end.x)-eps),x1=bin(std::max(start.x,end.x)+eps);
    int y0=bin(std::min(start.y,end.y)-eps),y1=bin(std::max(start.y,end.y)+eps);
    if((x1-x0+1)*(y1-y0+1)>64) {for(int ei:boundaryEdges_)intersect(ei);}
    else for(int y=y0;y<=y1;++y)for(int x=x0;x<=x1;++x)
        for(int ei:boundaryBins_[x+resolution_*y])intersect(ei);
    auto p=start+first*r;
    // Pull toward the originating fluid point to avoid leaving the domain by roundoff.
    if(first<1)p=start+(first*(1-1e-9))*r;
    return locateTriangle(p)>=0?p:clampToDomain(p);
}

SimplicialPoint2D SimplicialMesh2D::wallTangent(int vertex,SimplicialPoint2D velocity) const noexcept {
    const auto& boundary=vertices_[vertex].boundaryEdges;
    if(boundary.empty())return velocity;
    const auto n=edges_[boundary[0]].outwardNormal;
    for(int ei:boundary)if(std::abs(cross(n,edges_[ei].outwardNormal))>1e-8)return {};
    return velocity-dot(velocity,n)*n;
}

namespace {
std::uint64_t edgeKey(int a, int b) noexcept {
    const auto low = static_cast<std::uint32_t>(std::min(a, b));
    const auto high = static_cast<std::uint32_t>(std::max(a, b));
    return (std::uint64_t(low) << 32) | high;
}

SimplicialPoint2D circumcenter(
    SimplicialPoint2D a, SimplicialPoint2D b, SimplicialPoint2D c) {
    const double denominator = 2.0 * cross(b - a, c - a);
    if (!(denominator > 0.0)) throw std::invalid_argument("Simplicial triangle is inverted or degenerate");
    const double aa = dot(a, a), bb = dot(b, b), cc = dot(c, c);
    return {
        (aa * (b.y - c.y) + bb * (c.y - a.y) + cc * (a.y - b.y)) / denominator,
        (aa * (c.x - b.x) + bb * (a.x - c.x) + cc * (b.x - a.x)) / denominator
    };
}
}

SimplicialMesh2D::SimplicialMesh2D(int resolution, float boxSize)
    : resolution_(resolution), boxSize_(boxSize) {
    if (resolution < 2) throw std::invalid_argument("2D simplicial resolution must be at least two");
    const auto side = static_cast<std::size_t>(resolution) + 1;
    if (side > static_cast<std::size_t>(std::numeric_limits<int>::max()) / side)
        throw std::length_error("2D simplicial mesh exceeds the supported index range");
    setBoxSize(boxSize);
    buildTopology();
    updateGeometry();
}

DomainBounds SimplicialMesh2D::bounds() const noexcept {
    return {{0.f, 0.f, 0.f}, {boxSize_, boxSize_, 0.f}};
}

std::size_t SimplicialMesh2D::vertexIndex(int i, int j) const noexcept {
    return static_cast<std::size_t>(i) + (static_cast<std::size_t>(resolution_) + 1) * j;
}

SimplicialPoint2D SimplicialMesh2D::latticeToWorld(double s, double t) const noexcept {
    return {double(boxSize_) * s, double(boxSize_) * t};
}

SimplicialPoint2D SimplicialMesh2D::clampToDomain(SimplicialPoint2D point) const noexcept {
    if (!teapot_) return {std::clamp(point.x, 0.0, double(boxSize_)),
                           std::clamp(point.y, 0.0, double(boxSize_))};
    if (locateTriangle(point) >= 0) return point;
    double best = std::numeric_limits<double>::max();
    SimplicialPoint2D closest{};
    for (int ei : boundaryEdges_) {
        const auto& e = edges_[ei];
        const auto a = vertices_[e.first].position, d = vertices_[e.second].position - a;
        const auto q = a + std::clamp(dot(point-a,d)/dot(d,d),0.,1.) * d;
        const double distance = dot(point-q,point-q);
        if (distance < best) { best = distance; closest = q; }
    }
    return closest;
}

int SimplicialMesh2D::locateTriangle(SimplicialPoint2D point) const noexcept {
    if (!std::isfinite(point.x) || !std::isfinite(point.y) || point.x < 0 || point.y < 0 ||
        point.x > boxSize_ || point.y > boxSize_) return -1;
    if (teapot_) {
        const int x = std::min(resolution_-1, int(point.x / boxSize_ * resolution_));
        const int y = std::min(resolution_-1, int(point.y / boxSize_ * resolution_));
        for (int ti : bins_[x + resolution_ * y]) {
            const auto w = barycentric(ti, point);
            if (*std::min_element(w.begin(), w.end()) >= -1e-9) return ti;
        }
        return -1;
    }
    const double h = edgeScale();
    const int row = std::min(int(point.y / h), resolution_ - 1);
    const int column = std::min(int(point.x / h), resolution_ - 1);
    int best = 2 * (column + resolution_ * row);
    double bestMinimum = -std::numeric_limits<double>::infinity();
    // Odd rows are offset by half a cell; the containing quad is at most
    // one column away. Test actual triangle geometry, including wall strips.
    for (int i = std::max(0, column - 1); i <= std::min(resolution_ - 1, column + 1); ++i)
        for (int half = 0; half < 2; ++half) {
            const int triangle = 2 * (i + resolution_ * row) + half;
            const auto weights = barycentric(triangle, point);
            const double minimum = *std::min_element(weights.begin(), weights.end());
            if (minimum >= -1e-12) return triangle;
            if (minimum > bestMinimum) { bestMinimum = minimum; best = triangle; }
        }
    return best;
}

std::array<double, 3> SimplicialMesh2D::barycentric(
    int triangle, SimplicialPoint2D point) const noexcept {
    const auto& cell = triangles_[static_cast<std::size_t>(triangle)];
    const auto a = vertices_[cell.vertices[0]].position;
    const auto b = vertices_[cell.vertices[1]].position;
    const auto c = vertices_[cell.vertices[2]].position;
    const double denominator = cross(b - a, c - a);
    const double second = cross(point - a, c - a) / denominator;
    const double third = cross(b - a, point - a) / denominator;
    return {1.0 - second - third, second, third};
}

void SimplicialMesh2D::setBoxSize(float size) {
    if (!std::isfinite(size) || size <= 0.f || size / resolution_ <= 0.f)
        throw std::invalid_argument("2D simplicial mesh requires finite positive size and spacing");
    boxSize_ = size;
    if (!vertices_.empty()) updateGeometry();
}

void SimplicialMesh2D::buildTopology() {
    const int n = resolution_;
    vertices_.resize(std::size_t(n + 1) * (n + 1));
    triangles_.reserve(std::size_t(2) * n * n);
    for (int j = 0; j < n; ++j) for (int i = 0; i < n; ++i) {
        const int v00 = static_cast<int>(vertexIndex(i, j));
        const int v10 = static_cast<int>(vertexIndex(i + 1, j));
        const int v01 = static_cast<int>(vertexIndex(i, j + 1));
        const int v11 = static_cast<int>(vertexIndex(i + 1, j + 1));
        if (j % 2 == 0) {
            triangles_.push_back({{v00, v10, v01}});
            triangles_.push_back({{v10, v11, v01}});
        } else {
            triangles_.push_back({{v00, v10, v11}});
            triangles_.push_back({{v00, v11, v01}});
        }
    }

    unitPositions_.resize(vertices_.size());
    for (int j=0;j<=n;++j) for (int i=0;i<=n;++i)
        unitPositions_[vertexIndex(i,j)] = {(i+((j%2 && i>0 && i<n)?0.5:0.0))/n, double(j)/n};
    buildIncidence();
}

void SimplicialMesh2D::buildIncidence() {
    edges_.clear(); boundaryEdges_.clear();
    for (auto& vertex : vertices_) { vertex.boundary=false; vertex.incidentTriangles.clear(); vertex.boundaryEdges.clear(); }
    std::unordered_map<std::uint64_t, int> edgeLookup;
    for (int ti = 0; ti < static_cast<int>(triangles_.size()); ++ti) {
        auto& triangle = triangles_[ti];
        for (int local = 0; local < 3; ++local) {
            const int from = triangle.vertices[local];
            const int to = triangle.vertices[(local + 1) % 3];
            const auto key = edgeKey(from, to);
            auto [found, inserted] = edgeLookup.emplace(key, static_cast<int>(edges_.size()));
            if (inserted) edges_.push_back({std::min(from, to), std::max(from, to)});
            const int edge = found->second;
            triangle.edges[local] = edge;
            triangle.edgeSigns[local] = from < to ? 1 : -1;
            auto& incidence = edges_[edge].incidentTriangles;
            if (incidence[0] < 0) incidence[0] = ti;
            else if (incidence[1] < 0) incidence[1] = ti;
            else throw std::logic_error("Non-manifold edge in generated simplicial mesh");
        }
        for (const int vertex : triangle.vertices) vertices_[vertex].incidentTriangles.push_back(ti);
    }

    for (std::size_t edgeIndex = 0; edgeIndex < edges_.size(); ++edgeIndex) {
        auto& edge = edges_[edgeIndex];
        if (edge.boundary()) {
            vertices_[edge.first].boundary = true;
            vertices_[edge.second].boundary = true;
            vertices_[edge.first].boundaryEdges.push_back(int(edgeIndex));
            vertices_[edge.second].boundaryEdges.push_back(int(edgeIndex));
            boundaryEdges_.push_back(int(edgeIndex));
        } else {
            for (int side = 0; side < 2; ++side) {
                auto& triangle = triangles_[edge.incidentTriangles[side]];
                const auto local = std::find(triangle.edges.begin(), triangle.edges.end(),
                                             static_cast<int>(edgeIndex));
                if (local == triangle.edges.end())
                    throw std::logic_error("Broken triangle-edge incidence");
                triangle.neighbors[std::distance(triangle.edges.begin(), local)] =
                    edge.incidentTriangles[1 - side];
            }
        }
    }
}

void SimplicialMesh2D::updateGeometry() {
    for (size_t i=0;i<vertices_.size();++i) vertices_[i].position = double(boxSize_) * unitPositions_[i];

    for (auto& triangle : triangles_) {
        const auto a = vertices_[triangle.vertices[0]].position;
        const auto b = vertices_[triangle.vertices[1]].position;
        const auto c = vertices_[triangle.vertices[2]].position;
        triangle.area = 0.5 * cross(b - a, c - a);
        if (!(triangle.area > std::numeric_limits<double>::epsilon() * boxSize_ * boxSize_))
            throw std::invalid_argument("Simplicial triangle is near-degenerate");
        triangle.circumcenter = circumcenter(a, b, c);
    }

    for (auto& edge : edges_) {
        const auto a = vertices_[edge.first].position;
        const auto b = vertices_[edge.second].position;
        edge.primalLength = std::sqrt(dot(b - a, b - a));
        if (edge.boundary()) {
            const auto midpoint = 0.5 * (a + b);
            const auto delta = triangles_[edge.incidentTriangles[0]].circumcenter - midpoint;
            edge.dualLength = std::sqrt(dot(delta, delta));
        } else {
            const auto delta = triangles_[edge.incidentTriangles[0]].circumcenter -
                               triangles_[edge.incidentTriangles[1]].circumcenter;
            edge.dualLength = std::sqrt(dot(delta, delta));
        }
        if (edge.boundary()) {
            const auto& t = triangles_[edge.incidentTriangles[0]];
            auto center = (1./3) * (vertices_[t.vertices[0]].position + vertices_[t.vertices[1]].position + vertices_[t.vertices[2]].position);
            auto n = (1./edge.primalLength) * SimplicialPoint2D{b.y-a.y,a.x-b.x};
            edge.outwardNormal = dot(n,center-a)>0 ? -1.*n : n;
        }
        edge.hodge = edge.dualLength / edge.primalLength;
        if (!(edge.hodge > 1e-12) || !std::isfinite(edge.hodge))
            throw std::invalid_argument("Simplicial mesh has a non-positive or tiny dual edge");
    }

    for (auto& vertex : vertices_) if (!vertex.boundary) {
        std::sort(vertex.incidentTriangles.begin(), vertex.incidentTriangles.end(),
            [&](int first, int second) {
                const auto a = triangles_[first].circumcenter - vertex.position;
                const auto b = triangles_[second].circumcenter - vertex.position;
                return std::atan2(a.y, a.x) < std::atan2(b.y, b.x);
            });
    }
    if (teapot_) buildBins();
}
