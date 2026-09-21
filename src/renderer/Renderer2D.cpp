#include "renderer/Renderer2D.hpp"
#include "graphics/GlLoader.hpp"
#include "ShaderSources.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

namespace {
void validate(const ImageRenderData& data) {
    if (data.width <= 0 || data.height <= 0 || data.density.size() != std::size_t(data.width) * data.height ||
        !std::isfinite(data.worldWidth) || data.worldWidth <= 0.f ||
        !std::isfinite(data.worldHeight) || data.worldHeight <= 0.f)
        throw std::invalid_argument("Invalid 2D image dimensions or density extent");
}
GLuint compile(GLenum type, const char* source) {
    GLuint shader = glCreateShader(type);
    glShaderSource(shader, 1, &source, nullptr);
    glCompileShader(shader);
    GLint ok = 0;
    glGetShaderiv(shader, GL_COMPILE_STATUS, &ok);
    if (!ok) {
        char log[2048]{};
        glGetShaderInfoLog(shader, sizeof(log), nullptr, log);
        glDeleteShader(shader);
        throw std::runtime_error(std::string("2D shader compilation failed: ") + log);
    }
    return shader;
}
}
void Renderer2D::init(const ImageRenderData& data) {
    validate(data);
    shutdown();
    GLuint vertex = 0, fragment = 0;
    try {
        vertex = compile(GL_VERTEX_SHADER, ShaderSources::imageVert);
        fragment = compile(GL_FRAGMENT_SHADER, ShaderSources::imageFrag);
        program_ = glCreateProgram();
        glAttachShader(program_, vertex);
        glAttachShader(program_, fragment);
        glLinkProgram(program_);
        GLint ok = 0;
        glGetProgramiv(program_, GL_LINK_STATUS, &ok);
        if (!ok) {
            char log[2048]{};
            glGetProgramInfoLog(program_, sizeof(log), nullptr, log);
            throw std::runtime_error(std::string("2D shader linking failed: ") + log);
        }
        glGenVertexArrays(1, &vao_);
        glGenTextures(1, &texture_);
        glBindTexture(GL_TEXTURE_2D, texture_);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        upload(data);
    } catch (...) {
        if (vertex) glDeleteShader(vertex);
        if (fragment) glDeleteShader(fragment);
        shutdown();
        throw;
    }
    glDeleteShader(vertex);
    glDeleteShader(fragment);
}
void Renderer2D::upload(const ImageRenderData& data) {
    validate(data);
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, texture_);
    if (width_ != data.width || height_ != data.height) {
        GLint limit = 0;
        glGetIntegerv(GL_MAX_TEXTURE_SIZE, &limit);
        if (data.width > limit || data.height > limit) throw std::invalid_argument("2D image exceeds GPU texture limit");
        glTexImage2D(GL_TEXTURE_2D, 0, GL_R16F, data.width, data.height, 0, GL_RED, GL_FLOAT, data.density.data());
        width_ = data.width;
        height_ = data.height;
    } else {
        glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, data.width, data.height, GL_RED, GL_FLOAT, data.density.data());
    }
}
void Renderer2D::render(const ImageRenderData& data, int w, int h, float exposure) {
    if (!program_) throw std::logic_error("2D renderer is not initialized");
    if (!std::isfinite(exposure) || exposure < 0.f) throw std::invalid_argument("Invalid 2D exposure");
    if (w <= 0 || h <= 0) return;
    upload(data);
    // Fit physical bounds without stretching when the window aspect changes.
    const double aspect = double(data.worldWidth) / data.worldHeight;
    int vw = w, vh = h;
    if (double(w) / h > aspect) vw = std::max(1, int(h * aspect));
    else vh = std::max(1, int(w / aspect));
    glViewport((w - vw) / 2, (h - vh) / 2, vw, vh);
    glDisable(GL_DEPTH_TEST);
    glDisable(GL_CULL_FACE);
    glDisable(GL_BLEND);
    glDisable(GL_SCISSOR_TEST);
    glUseProgram(program_);
    glUniform1i(glGetUniformLocation(program_, "u_Density"), 0);
    glUniform1f(glGetUniformLocation(program_, "u_Exposure"), exposure);
    glUniform1i(glGetUniformLocation(program_, "u_HistoricalPalette"), data.historicalPalette);
    glBindVertexArray(vao_);
    glDrawArrays(GL_TRIANGLES, 0, 3);
    glBindVertexArray(0);
    glUseProgram(0);
}
void Renderer2D::shutdown() noexcept {
    if (program_) glDeleteProgram(program_);
    if (vao_) glDeleteVertexArrays(1, &vao_);
    if (texture_) glDeleteTextures(1, &texture_);
    program_ = vao_ = texture_ = 0;
    width_ = height_ = 0;
}
