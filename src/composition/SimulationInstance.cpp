#include "composition/SimulationInstance.hpp"

#include <stdexcept>
#include <utility>

SimulationInstance::SimulationInstance(
    SimulationConfig configuration,
    std::unique_ptr<IFluidSolver> solver,
    std::unique_ptr<IPresentationSource> presentation)
    : configuration_(configuration),
      solver_(std::move(solver)),
      presentation_(std::move(presentation)) {
    if (!solver_ || !presentation_)
        throw std::invalid_argument("A simulation instance requires a solver and presentation source");
    if (solver_->domain().dimension() != configuration_.dimension ||
        presentation_->dimension() != configuration_.dimension)
        throw std::invalid_argument("Solver, presentation, and configuration dimensions do not match");
}
