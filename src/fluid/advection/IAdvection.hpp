#pragma once

#include <string_view>

// Operation describes fields, geometry, tracing data and workspace for one
// discretization. A future mesh operation may carry circulation/cochain data;
// the generic interface does not require Cartesian face fields.
template<class Operation>
class IAdvection {
public:
    virtual ~IAdvection() = default;
    [[nodiscard]] virtual std::string_view name() const noexcept = 0;
    virtual void advect(const Operation& operation) const = 0;
};
