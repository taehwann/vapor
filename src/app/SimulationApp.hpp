#pragma once
#include "config/SimulationConfig.hpp"

class SimulationApp {
public:
    int run(const SimulationConfig& config);
private:
    int run2D(const SimulationConfig& config);
    int run3D(const SimulationConfig& config);
};
