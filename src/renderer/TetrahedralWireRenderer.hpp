#pragma once
#include "renderer/TetrahedralWireData.hpp"

class TetrahedralWireRenderer {
public:
    TetrahedralWireRenderer() = default;
    ~TetrahedralWireRenderer();
    TetrahedralWireRenderer(const TetrahedralWireRenderer&) = delete;
    TetrahedralWireRenderer& operator=(const TetrahedralWireRenderer&) = delete;
    void init(const TetrahedralWireData& data);
    void render(const float* mvp, int width, int height, bool primal, bool dual,
                bool barycentric, float cutX, bool depthTest);
private:
    unsigned int program_=0, vao_=0, buffer_=0;
    int primalCount_=0, circumCount_=0, baryCount_=0, surfaceCount_=0;
};

// Saves the current framebuffer before drawing UI; PPM needs no image library.
void saveWireframeImage(const std::filesystem::path& path, int width, int height);
