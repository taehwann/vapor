#include "solvers/mac2d/MacGridFluidSolver2D.hpp"
#include "solvers/mac2d/ReflectionMacFluidSolver2D.hpp"
#include "solvers/mac3d/MacGridFluidSolver3D.hpp"
#include "solvers/simplicial2d/SimplicialFluidSolver2D.hpp"
#include "solvers/simplicial3d/SimplicialFluidSolver3D.hpp"
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {
void require(bool ok, const char* message) {
    if (!ok)
        throw std::runtime_error(message);
}
} // namespace
int main() {
    try {
        MacGridFluidSolver2D simple;
        require(simple.state().resolution() == 128 && simple.boxSize() == 1.f, "Simple scene defaults");
        require(simple.parameters().projectIterations == 120 && simple.parameters().sourceStrength == 1.f,
                "Simple numerical defaults");
        simple.parameters().emitterKickMultiplier = 0;
        simple.parameters().buoyancy = 0;
        simple.advance(1.f / 60.f);
        require(simple.diagnostics().kineticEnergy == 0, "Disabled simple kick still moves fluid");

        for (auto method : {AdvectionMethod::SemiLagrangian, AdvectionMethod::MacCormack}) {
            ReflectionMacFluidSolver2D reflection(method);
            require(reflection.n == 128 && reflection.projectIterations == 120, "Reflection scene defaults");
            require(reflection.macCormackSmoke == (method == AdvectionMethod::MacCormack),
                    "Reflection smoke selection");
            reflection.advance(1.f / 60.f);
            require(reflection.diagnostics().kineticEnergy > 0, "Reflection default plume");
            reflection.resetState();
            reflection.emitterKickMultiplier = 0;
            reflection.buoyancy = 0;
            reflection.advance(1.f / 60.f);
            require(reflection.diagnostics().kineticEnergy == 0,
                    "Disabled reflection kick still moves fluid");
        }

        for (auto algorithm : {MacAlgorithm::SimpleSemiLagrangian, MacAlgorithm::SimpleMacCormack,
                               MacAlgorithm::ReflectionSemiLagrangian, MacAlgorithm::ReflectionMacCormack}) {
            MacGridFluidSolver3D solver(algorithm);
            const auto& p = solver.parameters();
            require(solver.resolution() == 32 && solver.boxSize() == 4.5f, "CPU 3D scene defaults");
            require(p.projectIterations == 600 && p.sorOmega == 1.9f, "CPU pressure defaults");
            const bool reflection = algorithm == MacAlgorithm::ReflectionSemiLagrangian ||
                                    algorithm == MacAlgorithm::ReflectionMacCormack;
            const bool corrected = algorithm == MacAlgorithm::SimpleMacCormack ||
                                   algorithm == MacAlgorithm::ReflectionMacCormack;
            require(p.reflection == reflection && p.macCormackSmoke == corrected &&
                        p.macCormackVel == corrected,
                    "CPU algorithm construction");
            solver.parameters().stirStrength = 0;
            solver.parameters().emitterKickMultiplier = 0;
            solver.parameters().buoyancy = 0;
            solver.advance(1.f / 60.f);
            require(solver.kineticEnergy() == 0, "Default CPU stages changed resting flow");
        }
        bool rejected = false;
        try {
            MacGridFluidSolver3D invalid(static_cast<MacAlgorithm>(-1));
        } catch (const std::invalid_argument&) {
            rejected = true;
        }
        require(rejected, "Invalid algorithm accepted");

        SimplicialFluidSolver2D triangle;
        require(triangle.boxSize() == 4.f && triangle.parameters().cgIterations == 600, "Triangle defaults");
        triangle.advance(1.f / 30.f);
        require(std::isfinite(triangle.diagnostics().kineticEnergy), "Triangle default step");

        SimplicialFluidSolver3D tetra(std::filesystem::path("examples/simplicial3d/box.tet"));
        require(tetra.boxSize() == 2.f &&
                    tetra.parameters().cgIterations == 1000,
                "Tetrahedral inviscid defaults");
        tetra.advance(1.f / 60.f);
        require(tetra.lastRecovery().converged, "Default tetrahedral recovery");
        std::cout << "PASS concrete solver construction and built-in numerical stages\n";
    } catch (const std::exception& error) {
        std::cerr << "FAIL construction: " << error.what() << '\n';
        return 1;
    }
}
