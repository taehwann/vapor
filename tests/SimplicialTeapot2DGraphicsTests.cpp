#include "solvers/simplicial2d/SimplicialFluidSolver2D.hpp"
#include "solvers/simplicial2d/SimplicialRasterizer2D.hpp"
#include "renderer/Renderer2D.hpp"
#include "renderer/SimplicialWireData2D.hpp"
#include "renderer/TetrahedralWireRenderer.hpp"
#include "graphics/GlLoader.hpp"
#include <iostream>

int main(int argc,char** argv) {
    if(!glfwInit())return 77;
    glfwWindowHint(GLFW_VISIBLE,GLFW_FALSE);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR,4);glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR,3);
    auto* window=glfwCreateWindow(1000,800,"Teapot validation",nullptr,nullptr);
    if(!window){glfwTerminate();return 77;}
    glfwMakeContextCurrent(window);
    if(!loadOpenGLFunctions()){glfwDestroyWindow(window);glfwTerminate();return 1;}
    int status=0;
    try {
        SimplicialFluidSolver2D solver{SimplicialFluidState2D(SimplicialMesh2D::teapot())};
        for(int i=0;i<180;++i)solver.advance(1.f/30);
        SimplicialRasterizer2D raster(solver.state());
        Renderer2D smoke;smoke.init(raster.renderData());
        TetrahedralWireRenderer wire,outline;
        wire.init(simplicialWires2D(solver.domain()));outline.init(simplicialWires2D(solver.domain(),true));
        float matrix[16];simplicialMvp2D(matrix,1000,800,4);
        for(int mode=0;mode<5;++mode) {
            const bool showSmoke=mode==0||mode==2;
            const bool showPrimal=mode!=0&&mode!=4;
            const bool showDual=mode!=0&&mode!=3;
            glClearColor(.025f,.025f,.035f,1);glClear(GL_COLOR_BUFFER_BIT);
            if(showSmoke)smoke.render(raster.renderData(),1000,800,3);
            if(mode!=0)wire.render(matrix,1000,800,showPrimal,showDual,false,5,false);
            if(mode<3)outline.render(matrix,1000,800,true,false,false,5,false);
            std::vector<unsigned char> pixels(1000*800*3);
            glReadPixels(0,0,1000,800,GL_RGB,GL_UNSIGNED_BYTE,pixels.data());
            size_t white=0,green=0,smokePixels=0;
            for(size_t i=0;i<pixels.size();i+=3) {
                if(pixels[i]>230&&pixels[i+1]>230&&pixels[i+2]>230)++white;
                if(pixels[i]<80&&pixels[i+1]>200&&pixels[i+2]<100)++green;
                if(pixels[i]>40&&pixels[i]<225&&pixels[i+2]>=pixels[i])++smokePixels;
            }
            if(glGetError()!=GL_NO_ERROR||(mode!=4&&white<500)||(showDual&&green<1000)||(showSmoke&&smokePixels<500)
               ||(mode==3&&green!=0)||(mode==4&&white!=0))
                throw std::runtime_error("Teapot smoke/wire render failed");
            // The hole's center must stay dark even though the mesh bounds contain it.
            int x=100+int(.20*800),y=int(.47*800),offset=3*(x+1000*y);
            if(pixels[offset]>15||pixels[offset+1]>15||pixels[offset+2]>25)
                throw std::runtime_error("Rendering filled the handle hole");
            if(argc>1)saveWireframeImage(std::string(argv[1])+"-"+std::to_string(mode)+".ppm",1000,800);
            std::cout<<"PASS teapot render "<<mode<<": white="<<white<<" green="<<green<<" smoke="<<smokePixels<<'\n';
        }
    }catch(const std::exception& e){std::cerr<<e.what()<<'\n';status=1;}
    glfwDestroyWindow(window);glfwTerminate();return status;
}
