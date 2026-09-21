#include "solvers/simplicial3d/SimplicialFluidSolver3D.hpp"
#include "renderer/SimplicialRenderData3D.hpp"
#include "renderer/VolumeRenderer.hpp"
#include "renderer/OrbitCamera.hpp"
#include "renderer/TetrahedralWireRenderer.hpp"
#include "graphics/GlLoader.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>

int main(int argc,char** argv) {
    if(!glfwInit()) return 77;
    glfwWindowHint(GLFW_VISIBLE,GLFW_FALSE);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR,4);glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR,3);
    auto* window=glfwCreateWindow(512,512,"Simplicial 3D validation",nullptr,nullptr);
    if(!window){glfwTerminate();return 77;}
    glfwMakeContextCurrent(window);
    if(!loadOpenGLFunctions()){glfwDestroyWindow(window);glfwTerminate();return 1;}
    int status=0;
    {
        VolumeRenderer renderer;
        try {
            SimplicialFluidSolver3D solver(simplicial3d::Mesh::loadBox("examples/simplicial3d/box.tet",2));solver.parameters().emitterRadius=.15f;
            // Rendering is tested with a manufactured density field. Sustained
            // inviscid dynamics remain covered by the separate stability tests.
            solver.setDensity([](simplicial3d::Point p) {
                auto d=p-simplicial3d::Point{1,1,1};return std::exp(-4*simplicial3d::dot(d,d));
            });
            SimplicialRenderData3D raster(48);raster.update(solver);renderer.init(raster.renderData());
            glViewport(0,0,512,512);glClearColor(.03f,.04f,.07f,1);glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
            float mvp[16],cx,cy,cz;buildMVP(mvp,cx,cy,cz,-1.2f,.3f,5.f,1.f,2.f);
            renderer.render(raster.renderData(),512,512,mvp,cx,cy,cz,.5f,30,.2683f,.3578f,.8944f,.3f,.1f,nullptr);
            std::vector<unsigned char> pixels(512*512*3);glReadBuffer(GL_BACK);
            glReadPixels(0,0,512,512,GL_RGB,GL_UNSIGNED_BYTE,pixels.data());
            if(glGetError()!=GL_NO_ERROR)throw std::runtime_error("OpenGL error");
            size_t lit=0;for(size_t i=0;i<pixels.size();i+=3)if(pixels[i]>32)++lit;
            if(lit<512)throw std::runtime_error("No substantial tetrahedral smoke in volume renderer");
            if(argc>1){std::ofstream out(argv[1],std::ios::binary);out<<"P6\n512 512\n255\n";
                for(int y=511;y>=0;--y)out.write(reinterpret_cast<char*>(pixels.data()+y*512*3),512*3);
                if(!out)throw std::runtime_error("Capture write failed");}
            auto d=solver.diagnostics();
            std::cout<<"PASS tetrahedral solver -> volume: lit pixels="<<lit<<" divergence="<<d.maxDivergence<<" energy="<<d.kineticEnergy<<'\n';
            TetrahedralWireRenderer wires;
            for (float size : {2.f, 3.f}) {
                if (size != solver.boxSize()) solver.setBoxSize(size);
                TetrahedralWireData data;
                for (auto e : solver.domain().edges)
                    data.primal.push_back({solver.domain().vertices[e.vertices[0]],solver.domain().vertices[e.vertices[1]]});
                for (auto e : solver.dual().edges)
                    data.circumcentric.push_back({solver.dual().vertices[e.a].position,solver.dual().vertices[e.b].position});
                wires.init(data);
                for(int mode=0;mode<3;++mode) {
                const bool primal=mode!=2,dual=mode!=1;
                glDepthMask(GL_TRUE);glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
                buildMVP(mvp,cx,cy,cz,-1.2f,.3f,size*2.5f,1.f,size);
                wires.render(mvp,512,512,primal,dual,false,size*.6f,true);
                glReadPixels(0,0,512,512,GL_RGB,GL_UNSIGNED_BYTE,pixels.data());
                size_t white=0,green=0;
                for(size_t i=0;i<pixels.size();i+=3){
                    if(pixels[i]>200&&pixels[i+1]>200&&pixels[i+2]>200)++white;
                    if(pixels[i]<80&&pixels[i+1]>200&&pixels[i+2]<100)++green;
                }
                if(glGetError()!=GL_NO_ERROR||(primal?white<1000:white!=0)||(dual?green<1000:green!=0))
                    throw std::runtime_error("Solver primal/dual wire rendering failed");
                std::cout<<"PASS solver wires at box size "<<size<<": white="<<white<<" green="<<green<<'\n';
                }
            }
        } catch(const std::exception& e){std::cerr<<"FAIL "<<e.what()<<'\n';status=1;}
    }
    glfwDestroyWindow(window);glfwTerminate();return status;
}
