#include "app/SimulationApp.hpp"
#include <exception>
#include <filesystem>
#include <iostream>

int main(int argc, char** argv) {
    try {
        const auto assets = std::filesystem::absolute(argv[0]).parent_path();
        SimulationApp app;
        return app.run(assets);
    } catch (const std::exception& error) {
        std::cerr << "Vapor: " << error.what() << '\n';
        return 1;
    }
}
