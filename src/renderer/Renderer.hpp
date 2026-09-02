#pragma once

#include "renderer/RenderData.hpp"

#include <cstdio>

// The rendering API depends on RenderData, never on a fluid solver.
class Renderer {
public:
    Renderer() = default;
    Renderer(const Renderer&) = delete;
    Renderer& operator=(const Renderer&) = delete;

    void init(const RenderData& data);

    void upload(const RenderData& data);

    void updateDomain(const RenderData& data);

    void render(const RenderData& data, int w, int h, const float* mvp, float cx, float cy, float cz, float stepScale, float alphaMul, float lightX, float lightY, float lightZ, float shadowStr, float shadowStep, FILE* logFile);

    void shutdown();

private:
    unsigned int prog = 0, vao = 0, vbo = 0, volumeTex = 0;
    int n = 0;
};
