#include "graphics/GlLoader.hpp"
#include "solvers/simplicial3d/gpu/GpuSimplicial3D.hpp"
#include "solvers/simplicial3d/SimplicialFluidSolver3D.hpp"
#include <chrono>
#include <iostream>
#include <cmath>
int main(int argc,char** argv) {
 if(!glfwInit())return 77;
 glfwWindowHint(GLFW_VISIBLE,GLFW_FALSE);
 glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR,4);glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR,3);
 auto* window=glfwCreateWindow(64,64,"Simplicial compute test",nullptr,nullptr);
 if(!window){glfwTerminate();return 77;}
 glfwMakeContextCurrent(window);if(!loadOpenGLFunctions())return 1;
 int status=0;
 try {
  auto mesh=simplicial3d::Mesh::loadDomain("examples/simplicial3d-bunny/bunny-fluid.tet",2);
  SimplicialFluidSolver3D cpu(mesh),gpu(mesh);GpuSimplicial3D backend;
  backend.initialize(gpu);gpu.setComputeBackend(&backend);
  std::cout<<backend.device()<<std::endl;
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
   if(i%10==0 || i+1==steps)std::cout<<i<<" fluxError="<<fluxError<<" densityError="<<densityError<<" cg="<<gpu.lastRecovery().iterations<<" gpuMs="<<gpuMs/(i+1)<<" cpuMs="<<cpuMs/(i+1)<<std::endl;
   if(!std::isfinite(d.maxSpeed)||d.maxBoundaryFlux>1e-10||d.maxDivergence>1e-8||fluxError>1e-5||densityError>1e-3)throw std::runtime_error("GPU parity or boundary check failed");
  }
 }catch(const std::exception& e){std::cerr<<e.what()<<std::endl;status=1;}
 glfwDestroyWindow(window);glfwTerminate();return status;
}
