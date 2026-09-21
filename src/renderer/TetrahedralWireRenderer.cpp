#include "renderer/TetrahedralWireRenderer.hpp"
#include "graphics/GlLoader.hpp"
#include <limits>

TetrahedralWireRenderer::~TetrahedralWireRenderer() {
    if (buffer_) glDeleteBuffers(1,&buffer_);
    if (vao_) glDeleteVertexArrays(1,&vao_);
    if (program_) glDeleteProgram(program_);
}
void TetrahedralWireRenderer::init(const TetrahedralWireData& data) {
    if (!program_) {
        const char* vertex=R"(#version 430 core
layout(location=0) in vec3 position;
uniform mat4 mvp;
out float worldX;
void main(){worldX=position.x;gl_Position=mvp*vec4(position,1);})";
        const char* fragment=R"(#version 430 core
in float worldX;
uniform float cutX;
uniform vec3 color;
out vec4 frag;
void main(){if(worldX>cutX)discard;frag=vec4(color,1);})";
        auto compile=[](unsigned int type,const char* source) {
            auto shader=glCreateShader(type);
            glShaderSource(shader,1,&source,nullptr);glCompileShader(shader);
            int ok;glGetShaderiv(shader,GL_COMPILE_STATUS,&ok);
            if(!ok){char log[2048];glGetShaderInfoLog(shader,2048,nullptr,log);glDeleteShader(shader);throw std::runtime_error(log);}
            return shader;
        };
        auto vs=compile(GL_VERTEX_SHADER,vertex);
        unsigned int fs=0;
        try { fs=compile(GL_FRAGMENT_SHADER,fragment); }
        catch (...) { glDeleteShader(vs); throw; }
        program_=glCreateProgram();glAttachShader(program_,vs);glAttachShader(program_,fs);glLinkProgram(program_);
        glDeleteShader(vs);glDeleteShader(fs);
        int ok;glGetProgramiv(program_,GL_LINK_STATUS,&ok);
        if(!ok)throw std::runtime_error("Wireframe shader link failed");
        glGenVertexArrays(1,&vao_);glGenBuffers(1,&buffer_);
    }
    std::vector<float> vertices;
    auto append=[&](const auto& segments) {
        for(auto s:segments) for(auto p:{s.a,s.b})
            for(double x:{p.x,p.y,p.z}) {
                if(!std::isfinite(x)||std::abs(x)>std::numeric_limits<float>::max())
                    throw std::invalid_argument("Non-finite wireframe position");
                vertices.push_back(float(x));
            }
    };
    append(data.primal);append(data.circumcentric);append(data.barycentric);
    primalCount_=int(data.primal.size()*2);circumCount_=int(data.circumcentric.size()*2);baryCount_=int(data.barycentric.size()*2);
    if (data.surfaceTriangles.size() % 3) throw std::invalid_argument("Incomplete surface triangle");
    surfaceCount_=int(data.surfaceTriangles.size());
    for (auto p : data.surfaceTriangles)
        for (double x : {p.x,p.y,p.z}) {
            if (!std::isfinite(x) || std::abs(x)>std::numeric_limits<float>::max())
                throw std::invalid_argument("Non-finite surface position");
            vertices.push_back(float(x));
        }
    glBindVertexArray(vao_);glBindBuffer(GL_ARRAY_BUFFER,buffer_);
    glBufferData(GL_ARRAY_BUFFER,vertices.size()*sizeof(float),vertices.data(),GL_STATIC_DRAW);
    glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,3*sizeof(float),nullptr);glEnableVertexAttribArray(0);
    glBindVertexArray(0);
}
void TetrahedralWireRenderer::render(const float* mvp,int width,int height,bool primal,bool dual,
                                   bool barycentric,float cutX,bool depthTest) {
    glViewport(0,0,width,height);glDisable(GL_CULL_FACE);glDisable(GL_BLEND);
    if(depthTest)glEnable(GL_DEPTH_TEST);else glDisable(GL_DEPTH_TEST);
    glDepthMask(GL_TRUE);glLineWidth(1);
    glUseProgram(program_);glBindVertexArray(vao_);
    glUniformMatrix4fv(glGetUniformLocation(program_,"mvp"),1,GL_FALSE,mvp);
    glUniform1f(glGetUniformLocation(program_,"cutX"),cutX);
    // Hide the far side of a surface wireframe without painting over smoke.
    if (depthTest && surfaceCount_) {
        glColorMask(GL_FALSE,GL_FALSE,GL_FALSE,GL_FALSE);
        glEnable(GL_POLYGON_OFFSET_FILL);glPolygonOffset(1,1);
        glDrawArrays(GL_TRIANGLES,primalCount_+circumCount_+baryCount_,surfaceCount_);
        glDisable(GL_POLYGON_OFFSET_FILL);
        glColorMask(GL_TRUE,GL_TRUE,GL_TRUE,GL_TRUE);
    }
    if(primal){glUniform3f(glGetUniformLocation(program_,"color"),1,1,1);glDrawArrays(GL_LINES,0,primalCount_);}
    if(dual){glUniform3f(glGetUniformLocation(program_,"color"),.1f,1,.25f);
        glDrawArrays(GL_LINES,primalCount_+(barycentric?circumCount_:0),barycentric?baryCount_:circumCount_);}
    glBindVertexArray(0);glUseProgram(0);glDisable(GL_DEPTH_TEST);
}
void saveWireframeImage(const std::filesystem::path& path,int width,int height) {
    if(width<=0||height<=0)throw std::runtime_error("Cannot capture an empty framebuffer");
    std::vector<unsigned char> pixels(size_t(width)*height*3);
    int alignment;glGetIntegerv(GL_PACK_ALIGNMENT,&alignment);glPixelStorei(GL_PACK_ALIGNMENT,1);
    glReadPixels(0,0,width,height,GL_RGB,GL_UNSIGNED_BYTE,pixels.data());glPixelStorei(GL_PACK_ALIGNMENT,alignment);
    std::ofstream out(path,std::ios::binary);out<<"P6\n"<<width<<' '<<height<<"\n255\n";
    for(int y=height-1;y>=0;--y)out.write(reinterpret_cast<const char*>(pixels.data()+size_t(y)*width*3),width*3);
    if(!out)throw std::runtime_error("Cannot save render: "+path.string());
}
