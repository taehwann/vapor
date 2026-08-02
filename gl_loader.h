#pragma once

// GLFW includes the platform OpenGL 1.1 header on Windows.  Functions added
// by later OpenGL versions must be resolved from the current context.
#include <GLFW/glfw3.h>
#include <cstddef>

// These typedefs first appeared in newer OpenGL headers than the legacy
// header supplied by Windows.
using GLchar = char;
using GLsizeiptr = std::ptrdiff_t;

#ifndef GL_ARRAY_BUFFER
#define GL_ARRAY_BUFFER 0x8892
#define GL_STATIC_DRAW 0x88E4
#define GL_TEXTURE0 0x84C0
#define GL_TEXTURE_3D 0x806F
#define GL_TEXTURE_WRAP_R 0x8072
#define GL_CLAMP_TO_EDGE 0x812F
#define GL_RED 0x1903
#define GL_R16F 0x822D
#define GL_RGBA8 0x8058
#define GL_FRAMEBUFFER 0x8D40
#define GL_COLOR_ATTACHMENT0 0x8CE0
#define GL_VERTEX_SHADER 0x8B31
#define GL_FRAGMENT_SHADER 0x8B30
#define GL_COMPILE_STATUS 0x8B81
#define GL_LINK_STATUS 0x8B82
#endif

#define VAPOR_GL_FUNCTIONS(X) \
    X(GLuint, glCreateShader, (GLenum)) \
    X(void, glShaderSource, (GLuint, GLsizei, const GLchar* const*, const GLint*)) \
    X(void, glCompileShader, (GLuint)) \
    X(void, glGetShaderiv, (GLuint, GLenum, GLint*)) \
    X(void, glGetShaderInfoLog, (GLuint, GLsizei, GLsizei*, GLchar*)) \
    X(void, glGetProgramiv, (GLuint, GLenum, GLint*)) \
    X(void, glGetProgramInfoLog, (GLuint, GLsizei, GLsizei*, GLchar*)) \
    X(void, glDeleteProgram, (GLuint)) \
    X(void, glGenVertexArrays, (GLsizei, GLuint*)) \
    X(void, glDeleteVertexArrays, (GLsizei, GLuint*)) \
    X(void, glBindVertexArray, (GLuint)) \
    X(void, glGenBuffers, (GLsizei, GLuint*)) \
    X(void, glDeleteBuffers, (GLsizei, GLuint*)) \
    X(void, glBindBuffer, (GLenum, GLuint)) \
    X(void, glBufferData, (GLenum, GLsizeiptr, const void*, GLenum)) \
    X(void, glVertexAttribPointer, (GLuint, GLint, GLenum, GLboolean, GLsizei, const void*)) \
    X(void, glEnableVertexAttribArray, (GLuint)) \
    X(GLuint, glCreateProgram, ()) \
    X(void, glAttachShader, (GLuint, GLuint)) \
    X(void, glLinkProgram, (GLuint)) \
    X(void, glDeleteShader, (GLuint)) \
    X(void, glTexImage3D, (GLenum, GLint, GLint, GLsizei, GLsizei, GLsizei, GLint, GLenum, GLenum, const void*)) \
    X(void, glGenFramebuffers, (GLsizei, GLuint*)) \
    X(void, glDeleteFramebuffers, (GLsizei, GLuint*)) \
    X(void, glBindFramebuffer, (GLenum, GLuint)) \
    X(void, glFramebufferTexture2D, (GLenum, GLenum, GLenum, GLuint, GLint)) \
    X(void, glUseProgram, (GLuint)) \
    X(GLint, glGetUniformLocation, (GLuint, const GLchar*)) \
    X(void, glUniform1i, (GLint, GLint)) \
    X(void, glUniform1f, (GLint, GLfloat)) \
    X(void, glUniform3f, (GLint, GLfloat, GLfloat, GLfloat)) \
    X(void, glActiveTexture, (GLenum))

#define VAPOR_DECLARE_GL_FUNCTION(return_type, name, arguments) \
    using name##Proc = return_type(APIENTRY*) arguments; \
    inline name##Proc name = nullptr;
VAPOR_GL_FUNCTIONS(VAPOR_DECLARE_GL_FUNCTION)
#undef VAPOR_DECLARE_GL_FUNCTION

inline bool loadOpenGLFunctions() {
    bool loaded = true;
#define VAPOR_LOAD_GL_FUNCTION(return_type, name, arguments) \
    name = reinterpret_cast<name##Proc>(glfwGetProcAddress(#name)); \
    loaded = loaded && name != nullptr;
    VAPOR_GL_FUNCTIONS(VAPOR_LOAD_GL_FUNCTION)
#undef VAPOR_LOAD_GL_FUNCTION
    return loaded;
}

#undef VAPOR_GL_FUNCTIONS
