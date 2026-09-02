#pragma once

#include "domain/IDomain.hpp"
#include "renderer/ImageRenderData.hpp"
#include "renderer/RenderData.hpp"

#include <variant>

using PresentationData = std::variant<ImageRenderData, RenderData>;

// Type-erased view of simulation state for the application. Concrete adapters
// remain responsible for translating each domain layout into renderer input.
class IPresentationSource {
public:
    virtual ~IPresentationSource() = default;
    [[nodiscard]] virtual Dimension dimension() const noexcept = 0;
    [[nodiscard]] virtual PresentationData presentationData() const noexcept = 0;
};
