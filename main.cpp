#include "app/SimulationApp.hpp"
#include <exception>
#include <iostream>

int main() {
    try {
        const SimulationConfig config = loadSimulationConfig("vapor.cfg");
        SimulationApp app;
        return app.run(config);
    } catch (const std::exception& error) {
        std::cerr << "Vapor: " << error.what() << '\n';
        return 1;
    }
}
