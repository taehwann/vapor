#include "GpuSimplicial3D.hpp"
#include "solvers/simplicial3d/SimplicialFluidSolver3D.hpp"
#include "graphics/GlLoader.hpp"
#include "ShaderSources.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <stdexcept>
using namespace simplicial3d;
namespace {
GLuint compile(const char* source){
 GLuint shader=glCreateShader(GL_COMPUTE_SHADER);glShaderSource(shader,1,&source,nullptr);glCompileShader(shader);
 GLint ok=0;glGetShaderiv(shader,GL_COMPILE_STATUS,&ok);
 if(!ok){char log[8192];glGetShaderInfoLog(shader,sizeof(log),nullptr,log);glDeleteShader(shader);throw std::runtime_error(std::string("Simplicial GPU shader: ")+log);}
 GLuint program=glCreateProgram();glAttachShader(program,shader);glLinkProgram(program);glDeleteShader(shader);glGetProgramiv(program,GL_LINK_STATUS,&ok);
 if(!ok){char log[8192];glGetProgramInfoLog(program,sizeof(log),nullptr,log);glDeleteProgram(program);throw std::runtime_error(std::string("Simplicial GPU link: ")+log);}return program;
}
void setInt(GLuint program,const char* name,int value){glUniform1i(glGetUniformLocation(program,name),value);}
void setDouble(GLuint program,const char* name,double value){glUniform1d(glGetUniformLocation(program,name),value);}
template<class T>void upload(GLuint buffer,const std::vector<T>& data,GLenum usage=GL_STATIC_DRAW){glBindBuffer(GL_SHADER_STORAGE_BUFFER,buffer);glBufferData(GL_SHADER_STORAGE_BUFFER,GLsizeiptr(data.size()*sizeof(T)),data.empty()?nullptr:data.data(),usage);}
void allocate(GLuint buffer,size_t bytes){glBindBuffer(GL_SHADER_STORAGE_BUFFER,buffer);glBufferData(GL_SHADER_STORAGE_BUFFER,GLsizeiptr(bytes),nullptr,GL_DYNAMIC_DRAW);}
void write(GLuint buffer,size_t offset,size_t bytes,const void* data){glBindBuffer(GL_SHADER_STORAGE_BUFFER,buffer);glBufferSubData(GL_SHADER_STORAGE_BUFFER,GLintptr(offset),GLsizeiptr(bytes),data);}
void read(GLuint buffer,size_t offset,size_t bytes,void* data){glBindBuffer(GL_SHADER_STORAGE_BUFFER,buffer);glGetBufferSubData(GL_SHADER_STORAGE_BUFFER,GLintptr(offset),GLsizeiptr(bytes),data);}
}
struct GpuSimplicial3D::Impl {
 enum Buffer {Geometry,Topology,Input,Traces,Output,Status,Matrix,Indices,CG,Count};
 std::array<GLuint,Count> buffers{};
 GLuint advect=0,recovery=0;int advStage=-1,cgStage=-1,n=0,nv=0,nd=0,nn=0,groups=0;
 const Tet* meshKey=nullptr;double size=0;std::string renderer;
 Impl(){glGenBuffers(Count,buffers.data());renderer=reinterpret_cast<const char*>(glGetString(GL_RENDERER));}
 ~Impl(){if(advect)glDeleteProgram(advect);if(recovery)glDeleteProgram(recovery);glDeleteBuffers(Count,buffers.data());}
 void bindAdvection(){glUseProgram(advect);for(int i=0;i<=Status;++i)glBindBufferBase(GL_SHADER_STORAGE_BUFFER,i,buffers[i]);}
 void bindRecovery(){glUseProgram(recovery);glBindBufferBase(GL_SHADER_STORAGE_BUFFER,0,buffers[Matrix]);glBindBufferBase(GL_SHADER_STORAGE_BUFFER,1,buffers[Indices]);glBindBufferBase(GL_SHADER_STORAGE_BUFFER,2,buffers[CG]);}
 void dispatchAdvection(int stage,int count){glUniform1i(advStage,stage);glDispatchCompute((count+63)/64,1,1);glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT|GL_BUFFER_UPDATE_BARRIER_BIT);}
 void dispatchCG(int stage,bool reduction=false){glUniform1i(cgStage,stage);glDispatchCompute(reduction?1:groups,1,1);glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT|GL_BUFFER_UPDATE_BARRIER_BIT);}
 void uploadState(const SimplicialFluidSolver3D& solver){
  static_assert(sizeof(Point)==3*sizeof(double));
  write(buffers[Input],0,solver.dualVelocity_.size()*sizeof(Point),solver.dualVelocity_.data());
  write(buffers[Input],size_t(nd)*sizeof(Point),solver.dual_.samples_.size()*sizeof(Point),solver.dual_.samples_.data());
  std::vector<double> density(solver.vertexDensity_.begin(),solver.vertexDensity_.end());
  write(buffers[Input],size_t(nd+nn)*sizeof(Point),density.size()*sizeof(double),density.data());
  int error=0;write(buffers[Status],0,sizeof(error),&error);
 }
 void check(){int error=0;read(buffers[Status],0,sizeof(error),&error);if(error)throw std::runtime_error("GPU sample escaped circumcentric subdivision");if(glGetError()!=GL_NO_ERROR)throw std::runtime_error("Simplicial GPU OpenGL error");}
};
GpuSimplicial3D::GpuSimplicial3D():impl_(std::make_unique<Impl>()) {
 impl_->advect=compile(ShaderSources::simplicialAdvection);impl_->recovery=compile(ShaderSources::simplicialRecovery);
 impl_->advStage=glGetUniformLocation(impl_->advect,"stage");impl_->cgStage=glGetUniformLocation(impl_->recovery,"stage");
}
GpuSimplicial3D::~GpuSimplicial3D()=default;
const std::string& GpuSimplicial3D::device() const{return impl_->renderer;}
void GpuSimplicial3D::initialize(const SimplicialFluidSolver3D& solver){
 auto& gpu=*impl_;const auto& mesh=solver.mesh_;const auto& dual=solver.dual_;
 if(mesh.boxDomain)throw std::invalid_argument("GPU simplicial backend requires the domain subdivision (select Bunny)");
 if(gpu.meshKey==mesh.tets.data()&&gpu.size==mesh.size)return;
 gpu.bindAdvection();std::vector<double> geometry;std::vector<int> topology;
 auto gd=[&](const char* name,const std::vector<double>& data){setInt(gpu.advect,name,int(geometry.size()));geometry.insert(geometry.end(),data.begin(),data.end());};
 auto gi=[&](const char* name,const std::vector<int>& data){setInt(gpu.advect,name,int(topology.size()));topology.insert(topology.end(),data.begin(),data.end());};
 auto point=[](std::vector<double>& data,Point p){data.insert(data.end(),{p.x,p.y,p.z});};
 std::vector<double> data;for(const auto& v:dual.vertices)point(data,v.position);gd("dualPosition",data);
 data.clear();for(const auto& node:dual.nodes_)point(data,node.position);gd("nodePosition",data);
 data.clear();std::vector<int> ids;
 for(const auto& tet:mesh.tets){point(data,mesh.vertices[tet.vertices[0]]);for(int k=1;k<4;++k)point(data,tet.gradients[k]);ids.insert(ids.end(),tet.vertices.begin(),tet.vertices.end());}
 gd("tetGeometry",data);gi("tetVertices",ids);gi("dualHint",solver.dualTets_);gi("primalHint",solver.primalTets_);
 auto bins=[&](const auto& cells,const char* offsets,const char* items){std::vector<int> ptr{0},flat;for(const auto& cell:cells){flat.insert(flat.end(),cell.begin(),cell.end());ptr.push_back(int(flat.size()));}gi(offsets,ptr);gi(items,flat);};
 bins(mesh.bins_,"binOffsets","binItems");bins(mesh.boundaryBins_,"wallOffsets","wallItems");
 data.clear();ids.clear();for(const auto& f:mesh.faces){point(data,f.areaNormal*(f.tets[0]>=0?1.:-1.));point(data,mesh.vertices[f.vertices[0]]);ids.push_back(f.tets[0]>=0?f.tets[0]:f.tets[1]);}gd("wallGeometry",data);gi("wallOwners",ids);
 std::vector<int> offsets{0};ids.clear();data.clear();int count=0;
 for(const auto& subdivisions:dual.subdivisions_){for(const auto& sub:subdivisions){ids.insert(ids.end(),sub.nodes.begin(),sub.nodes.end());for(auto g:sub.gradients)point(data,g);++count;}offsets.push_back(count);}
 gi("subOffsets",offsets);gi("subNodes",ids);gd("subGradient",data);
 offsets={0};ids.clear();for(const auto& loop:dual.loops){for(auto side:loop)ids.push_back(side.sign*(side.edge+1));offsets.push_back(int(ids.size()));}gi("loopOffsets",offsets);gi("loopItems",ids);
 ids.clear();for(auto edge:dual.edges){ids.push_back(edge.a);ids.push_back(edge.b);}gi("dualEdges",ids);
 gpu.n=int(mesh.edges.size());gpu.nv=int(mesh.vertices.size());gpu.nd=int(dual.vertices.size());gpu.nn=int(dual.nodes_.size());gpu.groups=(gpu.n+127)/128;
 setInt(gpu.advect,"nDual",gpu.nd);setInt(gpu.advect,"nPrimal",gpu.nv);setInt(gpu.advect,"nEdges",gpu.n);setInt(gpu.advect,"nNodes",gpu.nn);setInt(gpu.advect,"bins",mesh.binsPerAxis_);setDouble(gpu.advect,"size",mesh.size);
 upload(gpu.buffers[Impl::Geometry],geometry);upload(gpu.buffers[Impl::Topology],topology);
 allocate(gpu.buffers[Impl::Input],sizeof(double)*(3*size_t(gpu.nd+gpu.nn)+gpu.nv));allocate(gpu.buffers[Impl::Traces],sizeof(double)*6*gpu.nd);
 allocate(gpu.buffers[Impl::Output],sizeof(double)*std::max(gpu.n,gpu.nv));allocate(gpu.buffers[Impl::Status],sizeof(int));
 // Assemble the same Dirichlet potential operator as the CPU's matrix-free path.
 std::vector<std::map<int,double>> rows(gpu.n);
 auto add=[&](int i,int j,double value){if(!mesh.edges[i].boundary&&!mesh.edges[j].boundary)rows[i][j]+=value;};
 for(size_t f=0;f<mesh.faces.size();++f){const auto& face=mesh.faces[f];for(int a=0;a<3;++a)for(int b=0;b<3;++b)add(face.edges[a],face.edges[b],face.signs[a]*face.signs[b]*mesh.star2[f]);}
 std::vector<std::vector<std::pair<int,double>>> incident(mesh.vertices.size());
 for(int e=0;e<gpu.n;++e)if(!mesh.edges[e].boundary){auto v=mesh.edges[e].vertices;incident[v[0]].push_back({e,-mesh.star1[e]});incident[v[1]].push_back({e,mesh.star1[e]});}
 for(size_t v=0;v<incident.size();++v)if(!mesh.boundaryVertices[v])for(auto [i,a]:incident[v])for(auto [j,b]:incident[v])add(i,j,a*b/mesh.star0[v]);
 offsets={0};ids.clear();data.clear();std::vector<double> diagonal;
 for(int i=0;i<gpu.n;++i){if(mesh.edges[i].boundary)rows[i][i]=1;diagonal.push_back(rows[i].at(i));for(auto [j,value]:rows[i]){ids.push_back(j);data.push_back(value);}offsets.push_back(int(ids.size()));}
 int columns=int(offsets.size());offsets.insert(offsets.end(),ids.begin(),ids.end());data.insert(data.end(),diagonal.begin(),diagonal.end());
 upload(gpu.buffers[Impl::Matrix],data);upload(gpu.buffers[Impl::Indices],offsets);
 allocate(gpu.buffers[Impl::CG],sizeof(double)*(6*size_t(gpu.n)+2*gpu.groups+8));
 gpu.bindRecovery();setInt(gpu.recovery,"n",gpu.n);setInt(gpu.recovery,"groups",gpu.groups);setInt(gpu.recovery,"columns",columns);
 if(glGetError()!=GL_NO_ERROR)throw std::runtime_error("Could not allocate simplicial GPU buffers");
 gpu.meshKey=mesh.tets.data();gpu.size=mesh.size;
}
void GpuSimplicial3D::advectCirculation(const SimplicialFluidSolver3D& solver,double dt,std::vector<double>& omega){
 initialize(solver);auto& gpu=*impl_;gpu.bindAdvection();gpu.uploadState(solver);setDouble(gpu.advect,"dt",dt);
 gpu.dispatchAdvection(0,gpu.nd);gpu.dispatchAdvection(1,gpu.n);omega.resize(gpu.n);read(gpu.buffers[Impl::Output],0,omega.size()*sizeof(double),omega.data());gpu.check();
}
void GpuSimplicial3D::advectDensity(const SimplicialFluidSolver3D& solver,double dt,std::vector<float>& density){
 initialize(solver);auto& gpu=*impl_;gpu.bindAdvection();gpu.uploadState(solver);setDouble(gpu.advect,"dt",dt);setDouble(gpu.advect,"decay",std::exp(-solver.parameters_.smokeDecay*dt));
 gpu.dispatchAdvection(2,gpu.nv);std::vector<double> values(gpu.nv);read(gpu.buffers[Impl::Output],0,values.size()*sizeof(double),values.data());density.assign(values.begin(),values.end());gpu.check();
}
LinearSolveResult GpuSimplicial3D::recover(const SimplicialFluidSolver3D& solver,std::span<const double> rhs,std::vector<double>& x){
 initialize(solver);auto& gpu=*impl_;gpu.bindRecovery();
 const auto& parameters=solver.parameters_;if(parameters.cgIterations<0||!std::isfinite(parameters.cgAbsoluteTolerance)||parameters.cgAbsoluteTolerance<0||!std::isfinite(parameters.cgRelativeTolerance)||parameters.cgRelativeTolerance<0)throw std::invalid_argument("Invalid GPU CG parameters");
 std::vector<double> b(rhs.begin(),rhs.end());double norm=0;
 for(int i=0;i<gpu.n;++i){if(!std::isfinite(b[i]))throw std::runtime_error("Non-finite GPU RHS");if(solver.mesh_.edges[i].boundary)b[i]=0;norm=std::max(norm,std::abs(b[i]));}
 double target=std::max(parameters.cgAbsoluteTolerance,parameters.cgRelativeTolerance*norm);setDouble(gpu.recovery,"target",target);
 write(gpu.buffers[Impl::CG],0,gpu.n*sizeof(double),solver.phi_.data());write(gpu.buffers[Impl::CG],size_t(5)*gpu.n*sizeof(double),gpu.n*sizeof(double),b.data());
 gpu.dispatchCG(0);gpu.dispatchCG(1,true);std::array<double,8> status{};
 auto getStatus=[&]{read(gpu.buffers[Impl::CG],sizeof(double)*(6*size_t(gpu.n)+2*gpu.groups),sizeof(status),status.data());};getStatus();
 for(int done=0;done<parameters.cgIterations&&status[5]==0;){
  int batch=std::min(8,parameters.cgIterations-done);
  for(int i=0;i<batch;++i){gpu.dispatchCG(2);gpu.dispatchCG(3,true);gpu.dispatchCG(4);gpu.dispatchCG(5,true);gpu.dispatchCG(6);}done+=batch;getStatus();
 }
 if(status[7]!=0)throw std::runtime_error("GPU CG lost positive definiteness");
 LinearSolveResult result;result.initialResidual=status[3];result.iterations=int(status[6]);
 gpu.dispatchCG(7);gpu.dispatchCG(8,true);getStatus();result.finalResidual=status[4];result.converged=std::isfinite(status[4])&&status[5]==1;
 x.resize(gpu.n);read(gpu.buffers[Impl::CG],0,x.size()*sizeof(double),x.data());
 if(glGetError()!=GL_NO_ERROR)throw std::runtime_error("GPU CG OpenGL error");return result;
}
