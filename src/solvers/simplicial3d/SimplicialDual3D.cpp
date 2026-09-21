#include "SimplicialDual3D.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>

namespace simplicial3d {
namespace {
Point unit(Point p) {return p*(1/length(p));}
std::pair<Point,Point> tangentBasis(Point n) {
    Point u=unit(cross(n,std::abs(n.x)<.8?Point{1,0,0}:Point{0,1,0}));
    return {u,cross(n,u)};
}
}
DualMesh::DualMesh(const Mesh& mesh):size_(mesh.size) {
    if(!mesh.boxDomain){buildDomain(mesh);return;}
    const double tolerance=1e-8*size_;
    auto vertexId=[&](Point p,int tet=-1) {
        for(int i=0;i<int(vertices.size());++i) if(length(vertices[i].position-p)<tolerance) {
            if(tet>=0 && vertices[i].tet>=0 && vertices[i].tet!=tet)
                throw std::invalid_argument("Coincident circumcenters: regenerate a nondegenerate Delaunay mesh");
            if(tet>=0)vertices[i].tet=tet;return i;
        }
        unsigned walls=0;
        if(std::abs(p.x)<tolerance)walls|=1;if(std::abs(p.x-size_)<tolerance)walls|=2;
        if(std::abs(p.y)<tolerance)walls|=4;if(std::abs(p.y-size_)<tolerance)walls|=8;
        if(std::abs(p.z)<tolerance)walls|=16;if(std::abs(p.z-size_)<tolerance)walls|=32;
        vertices.push_back({p,tet,walls});return int(vertices.size())-1;
    };
    for(int t=0;t<int(mesh.tets.size());++t) vertexId(mesh.tets[t].center,t);
    cells.resize(mesh.vertices.size());loops.resize(mesh.edges.size());
    for(int e=0;e<int(mesh.edges.size());++e) {
        auto v=mesh.edges[e].vertices;Point a=mesh.vertices[v[0]],b=mesh.vertices[v[1]],n=unit(b-a);
        double offset=dot(n,(a+b)*.5);
        cells[v[0]].planes.push_back({n,offset,e});cells[v[1]].planes.push_back({n*(-1),-offset,e});
    }
    const std::array<DualPlane,6> walls{{{{-1,0,0},0},{{1,0,0},size_},{{0,-1,0},0},{{0,1,0},size_},{{0,0,-1},0},{{0,0,1},size_}}};
    std::map<std::array<int,2>,int> edgeIds;
    auto addEdge=[&](int a,int b) {
        std::array<int,2> key{std::min(a,b),std::max(a,b)};
        auto [it,inserted]=edgeIds.emplace(key,int(edges.size()));
        if(inserted) edges.push_back({key[0],key[1],(vertices[a].walls&vertices[b].walls)!=0});
        return OrientedDualEdge{it->second,a<b?1:-1};
    };
    for(int site=0;site<int(cells.size());++site) {
        auto& cell=cells[site];cell.planes.insert(cell.planes.end(),walls.begin(),walls.end());
        int count=int(cell.planes.size());
        for(int i=0;i<count;++i)for(int j=i+1;j<count;++j)for(int k=j+1;k<count;++k) {
            const auto& a=cell.planes[i];const auto& b=cell.planes[j];const auto& c=cell.planes[k];
            double det=dot(a.normal,cross(b.normal,c.normal));if(std::abs(det)<1e-10)continue;
            Point p=(cross(b.normal,c.normal)*a.offset+cross(c.normal,a.normal)*b.offset+cross(a.normal,b.normal)*c.offset)*(1/det);
            bool inside=true;for(const auto& plane:cell.planes)if(dot(plane.normal,p)>plane.offset+tolerance){inside=false;break;}
            if(!inside)continue;
            int id=vertexId(p);
            if(std::any_of(cell.corners.begin(),cell.corners.end(),[&](auto& corner){return corner.vertex==id;}))continue;
            std::vector<int> active;Point axis{};
            for(int f=0;f<count;++f)if(std::abs(dot(cell.planes[f].normal,p)-cell.planes[f].offset)<tolerance) {
                active.push_back(f);axis=axis+cell.planes[f].normal;
            }
            if(active.size()<3)throw std::runtime_error("Incomplete dual vertex");
            axis=unit(axis);auto [u,v]=tangentBasis(axis);
            std::sort(active.begin(),active.end(),[&](int a,int b){
                auto na=cell.planes[a].normal,nb=cell.planes[b].normal;
                return std::atan2(dot(na,v),dot(na,u))<std::atan2(dot(nb,v),dot(nb,u));
            });
            DualCorner corner{id,{}};
            for(size_t k=1;k+1<active.size();++k)corner.cones.push_back({active[0],active[k],active[k+1]});
            cell.corners.push_back(std::move(corner));
        }
        if(cell.corners.size()<4)throw std::invalid_argument("Empty clipped Voronoi cell");
        for(const auto& plane:cell.planes) {
            int e=plane.primalEdge;if(e<0||mesh.edges[e].vertices[0]!=site)continue;
            std::vector<int> polygon;Point center{};
            for(const auto& c:cell.corners)if(std::abs(dot(plane.normal,vertices[c.vertex].position)-plane.offset)<tolerance) {
                polygon.push_back(c.vertex);center=center+vertices[c.vertex].position;
            }
            if(polygon.size()<3)throw std::invalid_argument("Degenerate Voronoi face");
            center=center*(1./polygon.size());auto [u,v]=tangentBasis(plane.normal);
            std::sort(polygon.begin(),polygon.end(),[&](int a,int b){
                auto pa=vertices[a].position-center,pb=vertices[b].position-center;
                return std::atan2(dot(pa,v),dot(pa,u))<std::atan2(dot(pb,v),dot(pb,u));
            });
            for(size_t i=0;i<polygon.size();++i)loops[e].push_back(addEdge(polygon[i],polygon[(i+1)%polygon.size()]));
        }
    }
    faceEdges.resize(mesh.faces.size());
    for(size_t f=0;f<mesh.faces.size();++f) {
        auto face=mesh.faces[f];int a,b;
        a=face.tets[0]>=0?vertexId(mesh.tets[face.tets[0]].center):vertexId(face.center);
        b=face.tets[1]>=0?vertexId(mesh.tets[face.tets[1]].center):vertexId(face.center);
        auto it=edgeIds.find({std::min(a,b),std::max(a,b)});
        if(it==edgeIds.end())throw std::invalid_argument("Primal face has no matching clipped dual edge");
        faceEdges[f]=it->second;
        // Audit orientation: the loop boundary must implement d1^T exactly.
        int direction=a<b?1:-1;
        for(int k=0;k<3;++k) {
            int coefficient=0;for(auto edge:loops[face.edges[k]])if(edge.edge==it->second)coefficient+=edge.sign*direction;
            if(coefficient!=face.signs[k])throw std::runtime_error("Dual/primal incidence orientation mismatch");
        }
    }
}
std::vector<double> DualMesh::coordinates(int id,Point p) const {
    const auto& cell=cells.at(id);std::vector<double> weights(cell.corners.size(),0),distance(cell.planes.size());
    for(size_t v=0;v<weights.size();++v)if(length(vertices[cell.corners[v].vertex].position-p)<1e-12*size_) {weights[v]=1;return weights;}
    // Evaluate the boundary limit of the rational weights. This epsilon affects
    // only evaluation at a facet, not any dual measure or physical cochain.
    for(size_t f=0;f<distance.size();++f)distance[f]=std::max(1e-13*size_,cell.planes[f].offset-dot(cell.planes[f].normal,p));
    double sum=0;
    for(size_t v=0;v<weights.size();++v) {
        for(auto cone:cell.corners[v].cones) {
            Point a=cell.planes[cone[0]].normal,b=cell.planes[cone[1]].normal,c=cell.planes[cone[2]].normal;
            weights[v]+=std::abs(dot(a,cross(b,c)))/(distance[cone[0]]*distance[cone[1]]*distance[cone[2]]);
        }
        sum+=weights[v];
    }
    if(!std::isfinite(sum)||sum<=0)throw std::runtime_error("Invalid generalized barycentric weights");
    for(auto& w:weights)w/=sum;return weights;
}
Point DualMesh::interpolate(const Mesh& mesh,Point p,std::span<const Point> values,int tetHint) const {
    if(values.size()!=vertices.size())throw std::invalid_argument("Dual velocity size mismatch");
    if(!subdivisions_.empty()) {
        int ti=mesh.locate(p,tetHint);if(ti<0)return {};
        for(const auto& sub:subdivisions_[ti]) {
            Point delta=p-nodes_[sub.nodes[0]].position;
            std::array<double,4> w{0,dot(sub.gradients[0],delta),dot(sub.gradients[1],delta),dot(sub.gradients[2],delta)};
            w[0]=1-w[1]-w[2]-w[3];
            if(*std::min_element(w.begin(),w.end())>=-1e-8) {
                Point result{};for(int k=0;k<4;++k)result=result+samples_[sub.nodes[k]]*w[k];return result;
            }
        }
        throw std::runtime_error("Point outside circumcentric subdivision");
    }
    int site=0;double distance=std::numeric_limits<double>::max();
    for(int i=0;i<int(mesh.vertices.size());++i) {Point d=p-mesh.vertices[i];double squared=dot(d,d);if(squared<distance){site=i;distance=squared;}}
    auto weights=coordinates(site,p);Point value{};
    for(size_t i=0;i<weights.size();++i)value=value+values[cells[site].corners[i].vertex]*weights[i];
    return value;
}
std::vector<double> DualMesh::circulation(std::span<const double> integral) const {
    if(integral.size()!=edges.size())throw std::invalid_argument("Dual circulation size mismatch");
    std::vector<double> omega(loops.size(),0);
    for(size_t e=0;e<loops.size();++e)for(auto side:loops[e])omega[e]+=side.sign*integral[side.edge];
    return omega;
}
}

namespace simplicial3d {
Point tangentToNormals(Point p,std::span<const Point> normals) {
    std::vector<Point> basis;
    for(Point n:normals) {
        for(Point q:basis)n=n-q*dot(n,q);
        double l=length(n);if(l>1e-8)basis.push_back(n*(1/l));
    }
    for(Point q:basis)p=p-q*dot(p,q);
    return p;
}
void DualMesh::buildDomain(const Mesh& mesh) {
    const int nv=int(mesh.vertices.size()),ne=int(mesh.edges.size()),nf=int(mesh.faces.size()),nt=int(mesh.tets.size());
    loops.resize(ne);faceEdges.resize(nf);subdivisions_.resize(nt);
    nodes_.resize(nv+ne+nf+nt);
    std::vector<int> boundaryFace(nf,-1),boundaryEdge(ne,-1);
    for(int t=0;t<nt;++t) {
        auto w=mesh.barycentric(t,mesh.tets[t].center);
        if(*std::min_element(w.begin(),w.end())<=1e-8)
            throw std::invalid_argument("Domain interpolation requires strictly well-centered tetrahedra");
        vertices.push_back({mesh.tets[t].center,t,0,{}});
        nodes_[nv+ne+nf+t]={mesh.tets[t].center,{{t,1}}, {}};
        for(int v:mesh.tets[t].vertices)nodes_[v].weights.push_back({t,1});
    }
    for(int v=0;v<nv;++v)nodes_[v].position=mesh.vertices[v];
    for(int e=0;e<ne;++e) {
        auto ids=mesh.edges[e].vertices;
        nodes_[nv+e].position=(mesh.vertices[ids[0]]+mesh.vertices[ids[1]])*.5;
        if(mesh.edges[e].boundary) {
            boundaryEdge[e]=int(vertices.size());
            vertices.push_back({nodes_[nv+e].position,-1,1,{}});
            nodes_[nv+e].weights={{boundaryEdge[e],1}};
        }
    }
    for(int f=0;f<nf;++f) {
        const auto& face=mesh.faces[f];
        auto& node=nodes_[nv+ne+f];node.position=face.center;
        if(face.boundary()) {
            int id=boundaryFace[f]=int(vertices.size());Point n=face.areaNormal*(1/length(face.areaNormal));
            vertices.push_back({face.center,-1,1,{n}});node.weights={{id,1}};
            for(int e:face.edges)vertices[boundaryEdge[e]].normals.push_back(n);
            for(int v:face.vertices)nodes_[v].normals.push_back(n);
        } else {
            double a=length(face.center-mesh.tets[face.tets[0]].center),b=length(face.center-mesh.tets[face.tets[1]].center);
            node.weights={{face.tets[0],b/(a+b)},{face.tets[1],a/(a+b)}};
        }
        int a=face.tets[0]>=0?face.tets[0]:boundaryFace[f],b=face.tets[1]>=0?face.tets[1]:boundaryFace[f];
        faceEdges[f]=int(edges.size());edges.push_back({a,b,false});
        for(int k=0;k<3;++k)loops[face.edges[k]].push_back({faceEdges[f],face.signs[k]});
    }
    std::map<std::pair<int,int>,int> wallMidpoints;
    for(int e=0;e<ne;++e)if(mesh.edges[e].boundary) {
        std::map<int,int> balance;
        for(auto side:loops[e]){auto line=edges[side.edge];balance[line.a]+=side.sign;balance[line.b]-=side.sign;}
        int endpoints=0;
        for(auto [v,b]:balance)if(b) {
            if(std::abs(b)!=1||vertices[v].tet>=0)throw std::invalid_argument("Nonmanifold boundary dual face");
            Point midpoint=(vertices[boundaryEdge[e]].position+vertices[v].position)*.5;
            int mid=int(vertices.size());vertices.push_back({midpoint,-1,1,vertices[v].normals});
            int node=int(nodes_.size());nodes_.push_back({midpoint,{{mid,1}}, {}});
            wallMidpoints[{e,v}]=node;
            int id=int(edges.size());edges.push_back({boundaryEdge[e],mid,true});loops[e].push_back({id,b});
            id=int(edges.size());edges.push_back({mid,v,true});loops[e].push_back({id,b});++endpoints;
        }
        if(endpoints!=2)throw std::invalid_argument("Boundary edge must have two incident surface triangles");
    }
    for(int t=0;t<nt;++t) {
        const auto& tet=mesh.tets[t];
        for(int f:tet.faces)for(int e:mesh.faces[f].edges) {
            if(!mesh.edges[e].boundary)nodes_[nv+e].weights.push_back({t,1});
            for(int v:mesh.edges[e].vertices) {
                auto append=[&](std::array<int,4> ids) {
                    Subtet sub{ids,{}};
                    Point origin=nodes_[ids[0]].position,a=nodes_[ids[1]].position-origin,b=nodes_[ids[2]].position-origin,c=nodes_[ids[3]].position-origin;
                    double det=dot(a,cross(b,c));
                    if(std::abs(det)<1e-16*mesh.size*mesh.size*mesh.size)throw std::invalid_argument("Degenerate circumcentric subdivision");
                    sub.gradients={cross(b,c)*(1/det),cross(c,a)*(1/det),cross(a,b)*(1/det)};
                    subdivisions_[t].push_back(sub);
                };
                if(mesh.faces[f].boundary()) {
                    int mid=wallMidpoints.at({e,boundaryFace[f]});
                    append({v,nv+e,mid,nv+ne+nf+t});append({v,mid,nv+ne+f,nv+ne+nf+t});
                } else append({v,nv+e,nv+ne+f,nv+ne+nf+t});
            }
        }
    }
    for(auto& node:nodes_) {
        double total=0;for(auto [v,w]:node.weights)total+=w;
        if(total<=0)throw std::invalid_argument("Unconnected interpolation node");
        for(auto& entry:node.weights)entry.second/=total;
    }
    samples_.resize(nodes_.size());
    // Every circumcentric dual face must be a closed, oriented chain.
    for(const auto& loop:loops) {
        std::map<int,int> balance;
        for(auto side:loop){auto edge=edges[side.edge];balance[edge.a]+=side.sign;balance[edge.b]-=side.sign;}
        for(auto [v,b]:balance)if(b)throw std::invalid_argument("Open dual face");
    }
}
void DualMesh::prepareInterpolation(std::span<const Point> values) {
    if(values.size()!=vertices.size())throw std::invalid_argument("Dual velocity size mismatch");
    for(size_t i=0;i<nodes_.size();++i) {
        Point p{};for(auto [v,w]:nodes_[i].weights)p=p+values[v]*w;
        samples_[i]=tangentToNormals(p,nodes_[i].normals);
    }
}
}
