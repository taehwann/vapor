#include "SimplicialMesh3D.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>
#include <fstream>
#include <string>

namespace simplicial3d {
double length(Point a) { return std::sqrt(dot(a,a)); }
namespace {
Point solveRows(Point a, Point b, Point c, Point rhs) {
    const double det=dot(a,cross(b,c));
    if (std::abs(det)<1e-30) throw std::invalid_argument("Degenerate simplex");
    return (cross(b,c)*rhs.x+cross(c,a)*rhs.y+cross(a,b)*rhs.z)*(1/det);
}
void checkSize(size_t actual,size_t expected) {
    if(actual!=expected) throw std::invalid_argument("DEC cochain size mismatch");
}
}
Mesh::Mesh(std::vector<Point> v,std::vector<std::array<int,4>> t)
    :vertices(std::move(v)) { build(std::move(t)); }
Mesh Mesh::loadBox(const std::filesystem::path& path,double extent) {
    Mesh mesh=loadDomain(path,extent);mesh.boxDomain=true;
    double volume=0;
    for(const auto& t:mesh.tets) {
        volume+=t.volume;
        auto c=t.center;
        if(c.x< -1e-9*extent||c.y< -1e-9*extent||c.z< -1e-9*extent||c.x>extent*(1+1e-9)||c.y>extent*(1+1e-9)||c.z>extent*(1+1e-9))
            throw std::invalid_argument("Mesh circumcenter lies outside the box");
    }
    if(std::abs(volume-extent*extent*extent)>1e-8*extent*extent*extent)
        throw std::invalid_argument("Tetrahedra do not fill the box");
    for(const auto& face:mesh.faces)if(face.boundary()) {
        bool onBox=false;
        for(int axis=0;axis<3;++axis)for(double side:{0.,extent}) {
            bool onPlane=true;
            for(int v:face.vertices) {
                auto p=mesh.vertices[v];double coordinate=axis==0?p.x:axis==1?p.y:p.z;
                onPlane=onPlane&&std::abs(coordinate-side)<1e-10*extent;
            }
            if(onPlane)onBox=true;
        }
        if(!onBox)throw std::invalid_argument("Mesh has a boundary face away from the box boundary");
    }
    return mesh;
}
Mesh Mesh::loadDomain(const std::filesystem::path& path,double extent) {
    if(!std::isfinite(extent)||extent<=0) throw std::invalid_argument("Invalid mesh size");
    std::ifstream file(path); std::string magic;int version=0,nv=0,nt=0;
    if(!(file>>magic>>version>>nv>>nt)||magic!="VAPOR_TET"||version!=1||nv<4||nt<1||nv>1000000||nt>6000000)
        throw std::invalid_argument("Cannot read VAPOR_TET mesh: "+path.string());
    std::vector<Point> v(nv);std::vector<std::array<int,4>> cells(nt);
    for(auto& p:v) {
        if(!(file>>p.x>>p.y>>p.z)||p.x<0||p.x>1||p.y<0||p.y>1||p.z<0||p.z>1)
            throw std::invalid_argument("Mesh must contain finite unit-box coordinates");
        p=p*extent;
    }
    for(auto& t:cells) for(auto& i:t) if(!(file>>i)) throw std::invalid_argument("Truncated mesh connectivity");
    if(file>>magic) throw std::invalid_argument("Trailing mesh data");
    Mesh mesh(std::move(v),std::move(cells));mesh.size=extent;
    mesh.buildSpatialBins();return mesh;
}
Mesh Mesh::scaledBox(double extent) const {
    if(!std::isfinite(extent)||extent<=0)throw std::invalid_argument("Invalid box resize");
    auto v=vertices;for(auto& p:v)p=p*(extent/size);
    std::vector<std::array<int,4>> cells;for(const auto& t:tets)cells.push_back(t.vertices);
    Mesh mesh(std::move(v),std::move(cells));mesh.size=extent;mesh.boxDomain=boxDomain;mesh.buildSpatialBins();return mesh;
}
void Mesh::build(std::vector<std::array<int,4>> cells) {
    if(vertices.empty() || cells.empty()) throw std::invalid_argument("Empty tetrahedral mesh");
    for(auto p:vertices) if(!std::isfinite(p.x)||!std::isfinite(p.y)||!std::isfinite(p.z))
        throw std::invalid_argument("Non-finite mesh vertex");
    std::map<std::array<int,2>,int> edgeMap;
    std::map<std::array<int,3>,int> faceMap;
    std::map<std::array<int,4>,bool> unique;
    boundaryVertices.resize(vertices.size());
    minAltitude=std::numeric_limits<double>::max();
    for(auto ids:cells) {
        std::sort(ids.begin(),ids.end());
        if(ids[0]<0 || ids[3]>=int(vertices.size()) || std::adjacent_find(ids.begin(),ids.end())!=ids.end() || !unique.emplace(ids,true).second)
            throw std::invalid_argument("Invalid or duplicate tetrahedron");
        Tet t; t.vertices=ids;
        Point a=vertices[ids[0]], b=vertices[ids[1]]-a,c=vertices[ids[2]]-a,d=vertices[ids[3]]-a;
        double det=dot(b,cross(c,d));
        double scale=std::max({length(b),length(c),length(d)});
        if(std::abs(det)<1e-12*scale*scale*scale) throw std::invalid_argument("Degenerate tetrahedron");
        t.volume=std::abs(det)/6;
        t.gradients[1]=cross(c,d)*(1/det); t.gradients[2]=cross(d,b)*(1/det); t.gradients[3]=cross(b,c)*(1/det);
        t.gradients[0]=(t.gradients[1]+t.gradients[2]+t.gradients[3])*(-1);
        t.center=a+solveRows(b,c,d,{dot(b,b)/2,dot(c,c)/2,dot(d,d)/2});
        for(int k=0;k<4;++k) {
            std::array<int,3> fv; int j=0; for(int l=0;l<4;++l) if(l!=k) fv[j++]=ids[l];
            auto [it,inserted]=faceMap.emplace(fv,int(faces.size()));
            if(inserted) {
                Face f; f.vertices=fv;
                Point p=vertices[fv[0]],q=vertices[fv[1]],r=vertices[fv[2]];
                f.areaNormal=cross(q-p,r-p)*.5;
                f.center=p+solveRows(q-p,r-p,f.areaNormal,{dot(q-p,q-p)/2,dot(r-p,r-p)/2,0});
                std::array<std::array<int,2>,3> ev{{{fv[0],fv[1]},{fv[0],fv[2]},{fv[1],fv[2]}}};
                for(int l=0;l<3;++l) {
                    auto [ei,newEdge]=edgeMap.emplace(ev[l],int(edges.size()));
                    if(newEdge) edges.push_back({ev[l]});
                    f.edges[l]=ei->second;
                }
                faces.push_back(f);
            }
            int fi=it->second; Face& f=faces[fi];
            const int sign=dot(f.areaNormal,vertices[ids[k]]-vertices[fv[0]])<0 ? 1:-1;
            const int side=sign==1?0:1;
            if(f.tets[side]>=0) throw std::invalid_argument("Non-manifold or overlapping tetrahedra");
            f.tets[side]=int(tets.size()); t.faces[k]=fi; t.signs[k]=sign;
            minAltitude=std::min(minAltitude,3*t.volume/length(f.areaNormal));
        }
        tets.push_back(t);
    }
    star0.assign(vertices.size(),0); star1.assign(edges.size(),0); star2.assign(faces.size(),0); star3.resize(tets.size());
    for(size_t ti=0;ti<tets.size();++ti) {
        const auto& t=tets[ti]; star3[ti]=1/t.volume;
        for(int k=0;k<4;++k) {
            int fi=t.faces[k]; const auto& f=faces[fi];
            const double area=length(f.areaNormal);
            const double h=-t.signs[k]*dot(t.center-f.center,f.areaNormal)/area;
            star2[fi]+=h/area;
            for(int ei:f.edges) {
                auto e=edges[ei]; Point p=vertices[e.vertices[0]],q=vertices[e.vertices[1]],mid=(p+q)*.5;
                double dualArea=length(cross(f.center-mid,t.center-mid))*.5;
                double vol=std::abs(dot(mid-p,cross(f.center-p,t.center-p)))/6;
                int opposite=f.vertices[0]; for(int v:f.vertices) if(v!=e.vertices[0]&&v!=e.vertices[1]) opposite=v;
                double faceSign=dot(f.center-mid,vertices[opposite]-mid)>=0?1:-1;
                double tetSign=h>=0?1:-1;
                dualArea*=faceSign*tetSign; vol*=faceSign*tetSign;
                star1[ei]+=dualArea/length(q-p);
                star0[e.vertices[0]]+=vol; star0[e.vertices[1]]+=vol;
            }
        }
    }
    for(const auto& f:faces) if(f.boundary()) {
        for(int e:f.edges) edges[e].boundary=true;
        for(int v:f.vertices) boundaryVertices[v]=true;
    }
    int degree=0;
    for(const auto* star:{&star0,&star1,&star2,&star3}) {
        for(size_t i=0;i<star->size();++i)if(!std::isfinite((*star)[i])||(*star)[i]<=0)
            throw std::invalid_argument("Non-positive circumcentric star"+std::to_string(degree)+"["+std::to_string(i)+"]="+std::to_string((*star)[i])+"; regenerate the Delaunay mesh");
        ++degree;
    }
}
std::vector<double> Mesh::d0(std::span<const double> x) const {
    checkSize(x.size(),vertices.size()); std::vector<double> y(edges.size());
    for(size_t i=0;i<edges.size();++i) y[i]=x[edges[i].vertices[1]]-x[edges[i].vertices[0]]; return y;
}
std::vector<double> Mesh::d1(std::span<const double> x) const {
    checkSize(x.size(),edges.size()); std::vector<double> y(faces.size());
    for(size_t i=0;i<faces.size();++i) for(int k=0;k<3;++k) y[i]+=faces[i].signs[k]*x[faces[i].edges[k]]; return y;
}
std::vector<double> Mesh::d2(std::span<const double> x) const {
    checkSize(x.size(),faces.size()); std::vector<double> y(tets.size());
    for(size_t i=0;i<tets.size();++i) for(int k=0;k<4;++k) y[i]+=tets[i].signs[k]*x[tets[i].faces[k]]; return y;
}
std::vector<double> Mesh::d1Transpose(std::span<const double> x) const {
    checkSize(x.size(),faces.size()); std::vector<double> y(edges.size());
    for(size_t i=0;i<faces.size();++i) for(int k=0;k<3;++k) y[faces[i].edges[k]]+=faces[i].signs[k]*x[i]; return y;
}
std::array<double,4> Mesh::barycentric(int ti,Point p) const {
    auto& t=tets.at(ti); Point q=p-vertices[t.vertices[0]];
    std::array<double,4> w{0,dot(t.gradients[1],q),dot(t.gradients[2],q),dot(t.gradients[3],q)};
    w[0]=1-w[1]-w[2]-w[3]; return w;
}
int Mesh::locate(Point p,int hint) const {
    if(hint>=0&&hint<int(tets.size())) {
        auto w=barycentric(hint,p);
        if(*std::min_element(w.begin(),w.end())>=-1e-9)return hint;
    }
    if(binsPerAxis_>0) {
        if(p.x<0||p.y<0||p.z<0||p.x>size||p.y>size||p.z>size) return -1;
        int n=binsPerAxis_;
        int x=std::min(n-1,int(p.x/size*n)),y=std::min(n-1,int(p.y/size*n)),z=std::min(n-1,int(p.z/size*n));
        for(int i:bins_[x+n*(y+n*z)]) {auto w=barycentric(i,p); if(*std::min_element(w.begin(),w.end())>=-1e-9) return i;}
        return -1;
    }
    for(int k=0;k<int(tets.size());++k) {int i=(std::max(0,hint)+k)%int(tets.size()); auto w=barycentric(i,p); if(*std::min_element(w.begin(),w.end())>=-1e-9) return i;}
    return -1;
}
Point Mesh::clipSegment(Point start,Point end,int tetHint) const {
    if(tetHint>=0&&tetHint<int(tets.size())) {
        auto a=barycentric(tetHint,start),b=barycentric(tetHint,end);
        if(*std::min_element(a.begin(),a.end())>=-1e-10&&*std::min_element(b.begin(),b.end())>=-1e-10)return end;
    }
    Point delta=end-start;
    if(length(delta)<1e-14*size)return start;
    double first=1;
    auto inspect=[&](int fi) {
            const auto& face=faces[fi];
            int ti=face.tets[0]>=0?face.tets[0]:face.tets[1];
            Point normal=face.areaNormal*(face.tets[0]>=0?1.:-1.);
            double denominator=dot(normal,delta);
            if(denominator<=1e-15*length(normal)*length(delta))return;
            double hit=dot(normal,vertices[face.vertices[0]]-start)/denominator;
            if(hit< -1e-9||hit>first)return;
            auto weights=barycentric(ti,start+delta*hit);
            if(*std::min_element(weights.begin(),weights.end())< -1e-8)return;
            first=std::max(0.,hit);
    };
    if(binsPerAxis_>0) {
        int n=binsPerAxis_;auto bin=[&](double x){return std::clamp(int(x/size*n),0,n-1);};
        for(int z=bin(std::min(start.z,end.z));z<=bin(std::max(start.z,end.z));++z)
        for(int y=bin(std::min(start.y,end.y));y<=bin(std::max(start.y,end.y));++y)
        for(int x=bin(std::min(start.x,end.x));x<=bin(std::max(start.x,end.x));++x)
            for(int fi:boundaryBins_[x+n*(y+n*z)])inspect(fi);
    } else for(int fi=0;fi<int(faces.size());++fi)if(faces[fi].boundary())inspect(fi);
    return first<1 ? start+delta*std::max(0.,first-1e-8) : end;
}
void Mesh::buildSpatialBins() {
    int n=std::max(1,int(std::cbrt(tets.size()/6.)));binsPerAxis_=n;bins_.resize(n*n*n);
    for(int i=0;i<int(tets.size());++i) {
        Point lo{size,size,size},hi{};
        for(int v:tets[i].vertices) {auto p=vertices[v];lo={std::min(lo.x,p.x),std::min(lo.y,p.y),std::min(lo.z,p.z)};hi={std::max(hi.x,p.x),std::max(hi.y,p.y),std::max(hi.z,p.z)};}
        auto bin=[&](double x){return std::clamp(int(x/size*n),0,n-1);};
        for(int z=bin(lo.z);z<=bin(hi.z);++z)for(int y=bin(lo.y);y<=bin(hi.y);++y)for(int x=bin(lo.x);x<=bin(hi.x);++x)
            bins_[x+n*(y+n*z)].push_back(i);
    }
    boundaryBins_.resize(n*n*n);
    for(int i=0;i<int(faces.size());++i)if(faces[i].boundary()) {
        Point lo{size,size,size},hi{};
        for(int v:faces[i].vertices){auto p=vertices[v];lo={std::min(lo.x,p.x),std::min(lo.y,p.y),std::min(lo.z,p.z)};hi={std::max(hi.x,p.x),std::max(hi.y,p.y),std::max(hi.z,p.z)};}
        auto bin=[&](double x){return std::clamp(int(x/size*n),0,n-1);};
        for(int z=bin(lo.z);z<=bin(hi.z);++z)for(int y=bin(lo.y);y<=bin(hi.y);++y)for(int x=bin(lo.x);x<=bin(hi.x);++x)
            boundaryBins_[x+n*(y+n*z)].push_back(i);
    }
}
}
