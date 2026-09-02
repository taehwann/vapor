#include "renderer/Renderer.hpp"
#include "graphics/GlLoader.hpp"
#include "ShaderSources.hpp"

void Renderer::init(const RenderData& data) {
    n = data.width;
    const float s = data.boxSize;
    float cube[] = {
        // -Z face
        0,0,0,  0,s,0,  s,s,0,  0,0,0,  s,s,0,  s,0,0,
        // +Z face
        0,0,s,  s,0,s,  s,s,s,  0,0,s,  s,s,s,  0,s,s,
        // -Y face
        0,0,0,  s,0,0,  s,0,s,  0,0,0,  s,0,s,  0,0,s,
        // +Y face
        0,s,0,  0,s,s,  s,s,s,  0,s,0,  s,s,s,  s,s,0,
        // -X face
        0,0,0,  0,0,s,  0,s,s,  0,0,0,  0,s,s,  0,s,0,
        // +X face
        s,0,0,  s,s,0,  s,s,s,  s,0,0,  s,s,s,  s,0,s,
    };
    glGenVertexArrays(1, &vao);
    glGenBuffers(1, &vbo);
    glBindVertexArray(vao);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER, sizeof(cube), cube, GL_STATIC_DRAW);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);

    glGenTextures(1, &volumeTex);
    glBindTexture(GL_TEXTURE_3D, volumeTex);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);

    glTexImage3D(GL_TEXTURE_3D, 0, GL_R16F, n, n, n, 0, GL_RED, GL_FLOAT, nullptr);

    const char* vert = ShaderSources::vert;
    const char* frag = ShaderSources::frag;

    auto compile = [](GLenum type, const char* s) {
        GLuint sh = glCreateShader(type);
        glShaderSource(sh, 1, &s, nullptr); glCompileShader(sh);
        int ok = 0; glGetShaderiv(sh, GL_COMPILE_STATUS, &ok);
        if (!ok) { char log[1024]; glGetShaderInfoLog(sh, 1024, nullptr, log); std::fprintf(stderr, "GLSL %s compile: %s\n", type == GL_VERTEX_SHADER ? "VS" : "FS", log); }
        return sh;
        };
    GLuint vs = compile(GL_VERTEX_SHADER, vert), fs = compile(GL_FRAGMENT_SHADER, frag);
    prog = glCreateProgram();
    glAttachShader(prog, vs); glAttachShader(prog, fs); glLinkProgram(prog);
    int linkOk = 0; glGetProgramiv(prog, GL_LINK_STATUS, &linkOk);
    if (!linkOk) { char log[1024]; glGetProgramInfoLog(prog, 1024, nullptr, log); std::fprintf(stderr, "Link: %s\n", log); }
    glDeleteShader(vs); glDeleteShader(fs);
}

void Renderer::upload(const RenderData& data) {
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_3D, volumeTex);
    glTexSubImage3D(GL_TEXTURE_3D, 0, 0, 0, 0,
        data.width, data.height, data.depth, GL_RED, GL_FLOAT, data.density.data());
}

void Renderer::updateDomain(const RenderData& data) {
    const float s = data.boxSize;
    float cube[] = {
        // -Z face
        0,0,0,  0,s,0,  s,s,0,  0,0,0,  s,s,0,  s,0,0,
        // +Z face
        0,0,s,  s,0,s,  s,s,s,  0,0,s,  s,s,s,  0,s,s,
        // -Y face
        0,0,0,  s,0,0,  s,0,s,  0,0,0,  s,0,s,  0,0,s,
        // +Y face
        0,s,0,  0,s,s,  s,s,s,  0,s,0,  s,s,s,  s,s,0,
        // -X face
        0,0,0,  0,0,s,  0,s,s,  0,0,0,  0,s,s,  0,s,0,
        // +X face
        s,0,0,  s,s,0,  s,s,s,  s,0,0,  s,s,s,  s,0,s,
    };
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER, sizeof(cube), cube, GL_STATIC_DRAW);
}

void Renderer::render(const RenderData& data, int w, int h, const float* mvp, float cx, float cy, float cz, float stepScale, float alphaMul, float lightX, float lightY, float lightZ, float shadowStr, float shadowStep, FILE* logFile) {
    upload(data);
    glViewport(0, 0, w, h);
    glDisable(GL_DEPTH_TEST);
    glEnable(GL_CULL_FACE);
    glCullFace(GL_FRONT);
    glUseProgram(prog);
    glUniformMatrix4fv(glGetUniformLocation(prog, "u_MVP"), 1, GL_FALSE, mvp);
    glUniform3f(glGetUniformLocation(prog, "u_CamPos"), cx, cy, cz);
    glUniform1f(glGetUniformLocation(prog, "u_StepScale"), stepScale);
    glUniform1f(glGetUniformLocation(prog, "u_BoxSize"), data.boxSize);
    glUniform1f(glGetUniformLocation(prog, "u_GridRes"), float(data.width));
    glUniform1f(glGetUniformLocation(prog, "u_AlphaMul"), alphaMul);
    glUniform3f(glGetUniformLocation(prog, "u_LightDir"), lightX, lightY, lightZ);
    glUniform1f(glGetUniformLocation(prog, "u_ShadowStr"), shadowStr);
    glUniform1f(glGetUniformLocation(prog, "u_ShadowStep"), shadowStep);
    glUniform1i(glGetUniformLocation(prog, "u_Volume"), 0);
    glBindTexture(GL_TEXTURE_3D, volumeTex);
    glBindVertexArray(vao);
    glDrawArrays(GL_TRIANGLES, 0, 36);
    glBindVertexArray(0);
    glUseProgram(0);
    glDisable(GL_CULL_FACE);
    glCullFace(GL_BACK);
}

void Renderer::shutdown() {
    if (prog) glDeleteProgram(prog);
    if (vao) glDeleteVertexArrays(1, &vao);
    if (vbo) glDeleteBuffers(1, &vbo);
    if (volumeTex) glDeleteTextures(1, &volumeTex);
}
