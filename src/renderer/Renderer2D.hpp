#pragma once
#include "renderer/ImageRenderData.hpp"

// Orthographic scalar-image display. GL resources must be shut down before the context.
class Renderer2D {
public:
    Renderer2D() = default;
    ~Renderer2D() { shutdown(); }
    Renderer2D(const Renderer2D&) = delete;
    Renderer2D& operator=(const Renderer2D&) = delete;
    void init(const ImageRenderData& data);
    void render(const ImageRenderData& data, int framebufferWidth, int framebufferHeight, float exposure = 3.f);
    void shutdown() noexcept;
private:
    void upload(const ImageRenderData& data);
    unsigned int program_ = 0, vao_ = 0, texture_ = 0;
    int width_ = 0, height_ = 0;
};
