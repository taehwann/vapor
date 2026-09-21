#include "renderer/TetrahedralWireRenderer.hpp"
#include "renderer/OrbitCamera.hpp"
#include "graphics/GlLoader.hpp"
#include <iostream>

int main(int argc,char** argv) {
    if(!glfwInit())return 77;
    glfwWindowHint(GLFW_VISIBLE,GLFW_FALSE);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR,4);glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR,3);
    auto* window=glfwCreateWindow(800,800,"Bunny wire validation",nullptr,nullptr);
    if(!window){glfwTerminate();return 77;}
    glfwMakeContextCurrent(window);
    if(!loadOpenGLFunctions()){glfwDestroyWindow(window);glfwTerminate();return 1;}
    int status=0;
    try {
        TetrahedralWireRenderer renderer;
        renderer.init(TetrahedralWireData::load("examples/simplicial3d-bunny/bunny.tet",2));
        float mvp[16],cx,cy,cz;buildMVP(mvp,cx,cy,cz,-1.2f,.3f,4.f,1.f,2.f);
        for(int mode=0;mode<6;++mode) {
            glClearColor(.018f,.024f,.035f,1);glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
            bool primal=mode!=1,dual=mode!=0,bary=mode==3;
            renderer.render(mvp,800,800,primal,dual,bary,mode==4?1.f:mode==5?-100.f:1e6f,true);
            std::vector<unsigned char> pixels(800*800*3);
            glReadPixels(0,0,800,800,GL_RGB,GL_UNSIGNED_BYTE,pixels.data());
            size_t white=0,green=0;
            for(size_t i=0;i<pixels.size();i+=3){if(pixels[i]>200&&pixels[i+1]>200&&pixels[i+2]>200)++white;
                if(pixels[i]<80&&pixels[i+1]>200&&pixels[i+2]<100)++green;}
            if(glGetError()!=GL_NO_ERROR)throw std::runtime_error("Wireframe OpenGL error");
            if((mode==0&&(white<1000||green))||(mode==1&&(green<1000||white))||
               (mode>=2&&mode<=4&&(white<1000||green<1000))||(mode==5&&(white||green)))
                throw std::runtime_error("Wireframe color/toggle/cutaway check failed");
            if(argc>1&&mode==2)saveWireframeImage(argv[1],800,800);
            std::cout<<"PASS mode "<<mode<<": white="<<white<<" green="<<green<<'\n';
        }
    }catch(const std::exception& e){std::cerr<<e.what()<<'\n';status=1;}
    glfwDestroyWindow(window);glfwTerminate();return status;
}
