#pragma once
#include "presentation/IPresentationSource.hpp"
#include "renderer/ImageRenderData.hpp"

class MacGridState2D;

class MacGridPresentationAdapter2D final : public IPresentationSource {
public:
    // State must outlive the adapter and must not be moved from during use.
    explicit MacGridPresentationAdapter2D(const MacGridState2D& state) noexcept : state_(state) {}
    MacGridPresentationAdapter2D(const MacGridState2D&&) = delete;
    [[nodiscard]] ImageRenderData renderData() const noexcept;
    [[nodiscard]] Dimension dimension() const noexcept override { return Dimension::D2; }
    [[nodiscard]] PresentationData presentationData() const noexcept override { return renderData(); }
private:
    const MacGridState2D& state_;
};
