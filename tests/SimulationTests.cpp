#include "fluid/MacGridFluidSolver.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {
struct Baseline {
    int n, reflect, flags;
    double energy, maxVelocity, divergence, densitySum, weightedDensity;
};
// Captured from the pre-split main.cpp CPU implementation (MSVC Release).
// Six 1/60-second frames, 16 projection iterations, radius=.18, centerZ=.3,
// buoyancy=2.5, source=1, stirring disabled; flags: smoke MC=1, velocity MC=2.
const Baseline baselines[] = {
    {8, 0, 0, 5.32227563858, 1.79323649406, 0.00224227714352, 12.5633865102, 132.161955374},
    {8, 0, 1, 5.31826114655, 1.79280233383, 0.00225179060362, 12.3054398399, 130.187893328},
    {8, 0, 2, 5.63376235962, 1.80886173248, 0.00317340437323, 12.5620898544, 132.073204243},
    {8, 0, 3, 5.62872838974, 1.8083230257, 0.0031795501709, 12.3062464306, 130.151977868},
    {8, 1, 0, 5.3048157692, 1.79005801678, 0.0247040055692, 12.560423065, 132.109861564},
    {8, 1, 1, 5.30073595047, 1.78961789608, 0.0246352888644, 12.3041853862, 130.164586014},
    {8, 1, 2, 5.63356494904, 1.80607974529, 0.0130067132413, 12.5588216419, 132.010692464},
    {8, 1, 3, 5.62848520279, 1.80551791191, 0.0129552418366, 12.3049406262, 130.123926043},
    {12, 0, 0, 6.23518705368, 1.70936644077, 0.0195950269699, 47.985539313, 721.707724965},
    {12, 0, 1, 6.22765159607, 1.70904874802, 0.0196563154459, 46.408280425, 697.377404871},
    {12, 0, 2, 6.55779933929, 1.77157914639, 0.0131081836298, 47.9547437223, 721.24064451},
    {12, 0, 3, 6.55442667007, 1.76948368549, 0.0132885873318, 46.4151888065, 697.523262066},
    {12, 1, 0, 6.25715303421, 1.71856784821, 0.115855634212, 48.0296325346, 722.386514655},
    {12, 1, 1, 6.24975967407, 1.71815788746, 0.11610905081, 46.4359233264, 697.809387388},
    {12, 1, 2, 6.60593128204, 1.801410079, 0.0750236585736, 47.997528489, 721.888705373},
    {12, 1, 3, 6.60220718384, 1.79932284355, 0.0749910548329, 46.446071623, 698.00428824},
};
void near(double actual, double expected, const char* name) {
    const double tolerance = 2e-5 * std::max(1.0, std::abs(expected));
    if (!std::isfinite(actual) || std::abs(actual-expected)>tolerance) {
        throw std::runtime_error(std::string(name)+": actual="+std::to_string(actual)+
                                 ", baseline="+std::to_string(expected));
    }
}
}

int main() {
    try {
        for (const auto& baseline : baselines) {
            MacGridFluidSolver solver(baseline.n);
            solver.parameters().stirStrength=0.f;
            solver.parameters().emitterRadius=.18f;
            solver.parameters().emitterCenterZ=.3f;
            solver.parameters().macCormackSmoke=(baseline.flags & 1)!=0;
            solver.parameters().macCormackVel=(baseline.flags & 2)!=0;
            solver.parameters().reflection = baseline.reflect != 0;
            solver.parameters().projectIterations = 16;
            IFluidSolver& interface = solver;
            for (int frame=0; frame<6; ++frame) {
                interface.advance(1.f/60.f);
            }
            double sum=0, weighted=0;
            const auto density=solver.state().density();
            for (std::size_t i=0; i<density.size(); ++i) {
                sum+=density[i];
                weighted+=density[i]*double(1+i%31);
            }
            near(solver.kineticEnergy(),baseline.energy,"Kinetic energy");
            near(solver.maxVelocity(),baseline.maxVelocity,"Maximum velocity");
            near(solver.computeDivergenceNorm(),baseline.divergence,"Divergence");
            near(sum,baseline.densitySum,"Density sum");
            near(weighted,baseline.weightedDensity,"Weighted density");
            solver.resetState();
            near(solver.kineticEnergy(),0,"Reset energy");
            for (float value : solver.state().density()) near(value,0,"Reset density");
            solver.setBoxSize(3.f);
            near(interface.domain().bounds().max.x,3,"Physical domain update");
            std::cout << "PASS n=" << baseline.n << " reflect=" << baseline.reflect
                      << " MC flags=" << baseline.flags << '\n';
        }
    } catch (const std::exception& error) {
        std::cerr << "FAIL simulation regression: " << error.what() << '\n';
        return 1;
    }
    return 0;
}
