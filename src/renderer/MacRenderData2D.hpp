#pragma once
#include "renderer/ImageRenderData.hpp"
class MacGridState2D;
// Borrows density until the next state mutation or destruction.
ImageRenderData imageData(const MacGridState2D& state) noexcept;
