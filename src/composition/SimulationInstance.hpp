#pragma once

#include "config/SimulationConfig.hpp"
#include "fluid/IFluidSolver.hpp"
#include "presentation/IPresentationSource.hpp"

#include <memory>

// One application-ready simulation. The solver is the numerical abstraction;
// presentation remains a separate adapter with the same lifetime.
class SimulationInstance {
public:
    SimulationInstance(SimulationConfig configuration,
                       std::unique_ptr<IFluidSolver> solver,
                       std::unique_ptr<IPresentationSource> presentation);
    SimulationInstance(SimulationInstance&&) noexcept = default;
    SimulationInstance& operator=(SimulationInstance&&) noexcept = default;
    SimulationInstance(const SimulationInstance&) = delete;
    SimulationInstance& operator=(const SimulationInstance&) = delete;

    [[nodiscard]] const SimulationConfig& configuration() const noexcept { return configuration_; }
    [[nodiscard]] IFluidSolver& solver() noexcept { return *solver_; }
    [[nodiscard]] const IFluidSolver& solver() const noexcept { return *solver_; }
    [[nodiscard]] IPresentationSource& presentation() noexcept { return *presentation_; }
    [[nodiscard]] const IPresentationSource& presentation() const noexcept { return *presentation_; }

private:
    SimulationConfig configuration_;
    // Presentation is destroyed first because it may borrow solver state.
    std::unique_ptr<IFluidSolver> solver_;
    std::unique_ptr<IPresentationSource> presentation_;
};
