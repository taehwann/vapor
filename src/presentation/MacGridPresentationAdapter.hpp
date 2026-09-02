#pragma once

#include "presentation/IPresentationSource.hpp"
#include "renderer/RenderData.hpp"

class MacGridState;

// Converts MAC physical state to the existing volume renderer's input.
// The state must outlive this adapter and must not be moved from while in use.
class MacGridPresentationAdapter final : public IPresentationSource {
public:
    explicit MacGridPresentationAdapter(const MacGridState& state) noexcept;
    MacGridPresentationAdapter(const MacGridState&&) = delete;

    // Fresh metadata and a borrowed density view, without copying field data.
    // Consume the view before the next state mutation or destruction.
    [[nodiscard]] RenderData renderData() const noexcept;
    [[nodiscard]] Dimension dimension() const noexcept override { return Dimension::D3; }
    [[nodiscard]] PresentationData presentationData() const noexcept override { return renderData(); }

private:
    const MacGridState& state_;
};
