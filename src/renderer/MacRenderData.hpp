#pragma once
#include "renderer/RenderData.hpp"
class MacGridState3D;
// Borrows density until the next state mutation or destruction.
RenderData volumeData(const MacGridState3D& state) noexcept;
