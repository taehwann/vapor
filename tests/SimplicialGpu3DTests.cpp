#include "graphics/GlLoader.hpp"
#include "solvers/simplicial3d/gpu/GpuSimplicial3D.hpp"
#include "solvers/simplicial3d/SimplicialFluidSolver3D.hpp"
#include <chrono>
#include <iostream>
#include <cmath>
namespace {
std::pair<double,double> smokeCenter(const SimplicialFluidSolver3D& solver) {
 double mass=0,height=0;
 for(size_t i=0;i<solver.domain().tets.size();++i) {
  const auto& tet=solver.domain().tets[i];
  const double amount=solver.density()[i]*tet.volume;
  mass+=amount;
  for(int v:tet.vertices)height+=.25*amount*solver.domain().vertices[v].z;
 }
 return {mass>0?height/mass:0,mass};
}
}
int main(int argc,char** argv) {
 if(!glfwInit())return 77;
 glfwWindowHint(GLFW_VISIBLE,GLFW_FALSE);
 glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR,4);glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR,3);
 auto* window=glfwCreateWindow(64,64,"Simplicial compute test",nullptr,nullptr);
 if(!window){glfwTerminate();return 77;}
 glfwMakeContextCurrent(window);if(!loadOpenGLFunctions())return 1;
 int status=0;
 try {
  if(argc>1&&std::string(argv[1])=="motion") {
   const char* path=argc>2?argv[2]:"examples/simplicial3d/box-dense.tet";
   const int frames=argc>3?std::stoi(argv[3]):90;
   SimplicialFluidSolver3D solver(simplicial3d::Mesh::loadBox(path,3));
   GpuSimplicial3D backend;backend.initialize(solver);solver.setComputeBackend(&backend);
   solver.parameters().smokeDecay=0;
   double firstHeight=0,firstMass=0,elapsedMs=0;
   for(int i=0;i<frames;++i) {
    if(i==30){auto center=smokeCenter(solver);firstHeight=center.first;firstMass=center.second;solver.parameters().emitterEnabled=false;}
    const auto start=std::chrono::steady_clock::now();
    solver.advance(1.f/60);
    elapsedMs+=std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-start).count();
   }
   const auto [height,mass]=smokeCenter(solver);
   std::cout<<"GPU motion tets="<<solver.domain().tets.size()<<" frames="<<frames
            <<" center z="<<firstHeight<<" -> "<<height<<" mass="<<firstMass<<" -> "<<mass
            <<" speed="<<solver.diagnostics().maxSpeed<<" ms/frame="<<elapsedMs/frames<<std::endl;
   if(frames>=90 && !(height>firstHeight+.03))throw std::runtime_error("GPU smoke did not rise");
  } else {
  auto mesh=argc>2?simplicial3d::Mesh::loadBox(argv[2],3):simplicial3d::Mesh::loadDomain("examples/simplicial3d-bunny/bunny-fluid.tet",2);
  SimplicialFluidSolver3D cpu(mesh),gpu(mesh);GpuSimplicial3D backend;
  backend.initialize(gpu);gpu.setComputeBackend(&backend);
  std::cout<<backend.device()<<std::endl;
  if(mesh.boxDomain) {
   std::vector<double> potential(mesh.edges.size());
   for(size_t e=0;e<potential.size();++e) {
    auto ids=mesh.edges[e].vertices;auto a=mesh.vertices[ids[0]],b=mesh.vertices[ids[1]],p=(a+b)*.5;
    potential[e]=std::sin(3.141592653589793*p.x/mesh.size)*std::sin(3.141592653589793*p.y/mesh.size)*std::sin(3.141592653589793*p.z/mesh.size)*(b.z-a.z);
   }
   cpu.setPotential(potential);gpu.setPotential(potential);
   cpu.advectVorticity(.017);gpu.advectVorticity(.017);
   double error=0;for(size_t i=0;i<cpu.vorticity().size();++i)error=std::max(error,std::abs(cpu.vorticity()[i]-gpu.vorticity()[i]));
   if(error>1e-8)throw std::runtime_error("GPU box seeded circulation differs from CPU");
   std::cout<<"seeded circulation error="<<error<<std::endl;
   cpu.resetState();gpu.resetState();
  }
  int steps=argc>1?std::stoi(argv[1]):60;double cpuMs=0,gpuMs=0;
  for(int i=0;i<steps;++i){
   auto a=std::chrono::steady_clock::now();cpu.advance(1.f/30);
   auto b=std::chrono::steady_clock::now();gpu.advance(1.f/30);
   auto c=std::chrono::steady_clock::now();
   cpuMs+=std::chrono::duration<double,std::milli>(b-a).count();gpuMs+=std::chrono::duration<double,std::milli>(c-b).count();
   double fluxError=0,densityError=0;
   for(size_t j=0;j<cpu.flux().size();++j)fluxError=std::max(fluxError,std::abs(cpu.flux()[j]-gpu.flux()[j]));
   for(size_t j=0;j<cpu.vertexDensity().size();++j)densityError=std::max(densityError,double(std::abs(cpu.vertexDensity()[j]-gpu.vertexDensity()[j])));
   auto d=gpu.diagnostics();
   if(i%10==0 || i+1==steps)std::cout<<i<<" fluxError="<<fluxError<<" densityError="<<densityError<<" cg="<<gpu.lastRecovery().iterations<<" substeps="<<gpu.timings().substeps<<" gpuMs="<<gpuMs/(i+1)<<" cpuMs="<<cpuMs/(i+1)<<std::endl;
   if(!std::isfinite(d.maxSpeed)||d.maxBoundaryFlux>1e-10||d.maxDivergence>1e-8||fluxError>1e-5||densityError>1e-3)throw std::runtime_error("GPU parity or boundary check failed");
  }
  const auto& timing=gpu.timings();
  std::cout<<"last GPU-path step ms: circulation="<<timing.advectionMs<<" forces="<<timing.forcesMs<<" recovery="<<timing.recoveryMs<<" smoke="<<timing.smokeMs<<std::endl;
  cpu.setBoxSize(2.5f);gpu.setBoxSize(2.5f);
  cpu.advance(.01f);gpu.advance(.01f);
  double resizeError=0;for(size_t i=0;i<cpu.flux().size();++i)resizeError=std::max(resizeError,std::abs(cpu.flux()[i]-gpu.flux()[i]));
  if(resizeError>1e-5||!gpu.lastRecovery().converged)throw std::runtime_error("GPU resize failed parity/recovery");
  auto replay=gpu.flux();gpu.resetState();gpu.advance(.01f);
  double replayError=0;for(size_t i=0;i<replay.size();++i)replayError=std::max(replayError,std::abs(replay[i]-gpu.flux()[i]));
  if(replayError>1e-10)throw std::runtime_error("GPU reset replay mismatch");
  std::cout<<"PASS resize and reset; resize flux error="<<resizeError<<" replay error="<<replayError<<std::endl;
  }
 }catch(const std::exception& e){std::cerr<<e.what()<<std::endl;status=1;}
 glfwDestroyWindow(window);glfwTerminate();return status;
}
