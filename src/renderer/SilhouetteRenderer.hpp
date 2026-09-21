#pragma once
#include "graphics/GlLoader.hpp"
#include "solvers/simplicial3d/SimplicialMesh3D.hpp"
#include <stdexcept>
#include <vector>

// Project the exterior into a separate mask; composite only its screen-space
// contour. No surface color or depth is ever written over the smoke.
class SilhouetteRenderer {
    GLuint maskProgram_=0, outlineProgram_=0, vao_=0, buffer_=0, texture_=0, fbo_=0;
    int count_=0, width_=0, height_=0;
    static GLuint program(const char* vertex, const char* fragment) {
        auto compile=[](GLenum type,const char* source) {
            GLuint shader=glCreateShader(type);
            glShaderSource(shader,1,&source,nullptr);glCompileShader(shader);
            GLint ok;glGetShaderiv(shader,GL_COMPILE_STATUS,&ok);
            if(!ok){char log[2048];glGetShaderInfoLog(shader,2048,nullptr,log);glDeleteShader(shader);throw std::runtime_error(log);}
            return shader;
        };
        GLuint vs=compile(GL_VERTEX_SHADER,vertex),fs=0;
        try {fs=compile(GL_FRAGMENT_SHADER,fragment);} catch(...){glDeleteShader(vs);throw;}
        GLuint p=glCreateProgram();glAttachShader(p,vs);glAttachShader(p,fs);glLinkProgram(p);
        glDeleteShader(vs);glDeleteShader(fs);
        GLint ok;glGetProgramiv(p,GL_LINK_STATUS,&ok);
        if(!ok){glDeleteProgram(p);throw std::runtime_error("Silhouette shader link failed");}
        return p;
    }
public:
    SilhouetteRenderer()=default;
    SilhouetteRenderer(const SilhouetteRenderer&)=delete;
    SilhouetteRenderer& operator=(const SilhouetteRenderer&)=delete;
    ~SilhouetteRenderer(){
        if(maskProgram_)glDeleteProgram(maskProgram_);
        if(outlineProgram_)glDeleteProgram(outlineProgram_);
        if(buffer_)glDeleteBuffers(1,&buffer_);
        if(vao_)glDeleteVertexArrays(1,&vao_);
        if(texture_)glDeleteTextures(1,&texture_);
        if(fbo_)glDeleteFramebuffers(1,&fbo_);
    }
    void init(const simplicial3d::Mesh& mesh) {
        if(!maskProgram_) {
            maskProgram_=program(R"(#version 430 core
layout(location=0) in vec3 position;
uniform mat4 mvp;
void main(){gl_Position=mvp*vec4(position,1);})",R"(#version 430 core
out vec4 color;
void main(){color=vec4(1);})");
            outlineProgram_=program(R"(#version 430 core
void main(){vec2 p=vec2((gl_VertexID<<1)&2,gl_VertexID&2);gl_Position=vec4(p*2-1,0,1);})",R"(#version 430 core
uniform sampler2D maskTexture;
out vec4 color;
float maskAt(ivec2 p){ivec2 s=textureSize(maskTexture,0);if(any(lessThan(p,ivec2(0)))||any(greaterThanEqual(p,s)))return 0;return texelFetch(maskTexture,p,0).r;}
void main(){
 ivec2 p=ivec2(gl_FragCoord.xy);
 if(maskAt(p)>.5)discard;
 float neighbor=0;
 for(int y=-1;y<=1;++y)for(int x=-1;x<=1;++x)neighbor=max(neighbor,maskAt(p+ivec2(x,y)));
 if(neighbor<.5)discard;
 color=vec4(1);
})");
            glGenVertexArrays(1,&vao_);glGenBuffers(1,&buffer_);
            glGenTextures(1,&texture_);glGenFramebuffers(1,&fbo_);
        }
        std::vector<float> vertices;
        for(const auto& f:mesh.faces)if(f.boundary())for(int v:f.vertices) {
            auto p=mesh.vertices[v];vertices.insert(vertices.end(),{float(p.x),float(p.y),float(p.z)});
        }
        count_=int(vertices.size()/3);
        glBindVertexArray(vao_);glBindBuffer(GL_ARRAY_BUFFER,buffer_);
        glBufferData(GL_ARRAY_BUFFER,vertices.size()*sizeof(float),vertices.data(),GL_STATIC_DRAW);
        glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,3*sizeof(float),nullptr);glEnableVertexAttribArray(0);
        glBindVertexArray(0);
    }
    void render(const float* mvp,int width,int height) {
        if(width<=0||height<=0)return;
        constexpr GLenum drawFramebufferBinding=0x8CA6, red8=0x8229;
        GLint destination;glGetIntegerv(drawFramebufferBinding,&destination);
        glActiveTexture(GL_TEXTURE0);glBindTexture(GL_TEXTURE_2D,texture_);
        if(width!=width_||height!=height_) {
            glTexImage2D(GL_TEXTURE_2D,0,red8,width,height,0,GL_RED,GL_UNSIGNED_BYTE,nullptr);
            glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_NEAREST);
            glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_NEAREST);
            width_=width;height_=height;
        }
        glBindFramebuffer(GL_FRAMEBUFFER,fbo_);
        glFramebufferTexture2D(GL_FRAMEBUFFER,GL_COLOR_ATTACHMENT0,GL_TEXTURE_2D,texture_,0);
        glViewport(0,0,width,height);glDisable(GL_DEPTH_TEST);glDisable(GL_CULL_FACE);glDisable(GL_BLEND);
        glDepthMask(GL_FALSE);glClearColor(0,0,0,0);glClear(GL_COLOR_BUFFER_BIT);
        glUseProgram(maskProgram_);glUniformMatrix4fv(glGetUniformLocation(maskProgram_,"mvp"),1,GL_FALSE,mvp);
        glBindVertexArray(vao_);glDrawArrays(GL_TRIANGLES,0,count_);
        glBindFramebuffer(GL_FRAMEBUFFER,GLuint(destination));
        glUseProgram(outlineProgram_);glUniform1i(glGetUniformLocation(outlineProgram_,"maskTexture"),0);
        glDrawArrays(GL_TRIANGLES,0,3);
        glBindVertexArray(0);glUseProgram(0);glDepthMask(GL_TRUE);
    }
};
