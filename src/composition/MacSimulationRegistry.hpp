#pragma once

#include "config/SimulationConfig.hpp"
#include "domain/MacGridDomain.hpp"
#include "domain/MacGridDomain2D.hpp"
#include "fluid/MacGridFluidSolver.hpp"
#include "fluid/MacGridFluidSolver2D.hpp"

#include <cstddef>
#include <functional>
#include <memory>
#include <unordered_map>

struct MacAdvectionKey {
    DomainType domain;
    AdvectionType method;
    ComputeBackend backend;
    bool operator==(const MacAdvectionKey&) const noexcept = default;
};

struct MacProjectionKey {
    DomainType domain;
    ProjectionSolverType method;
    bool operator==(const MacProjectionKey&) const noexcept = default;
};

struct MacFluidKey {
    DomainType domain;
    FluidSolverType method;
    bool operator==(const MacFluidKey&) const noexcept = default;
};

// Modules are destroyed in reverse order: solver first, then the borrowed
// pressure/advection implementations.
struct MacSimulation2D {
    std::unique_ptr<IAdvection<MacAdvectionOperation2D>> advection;
    std::unique_ptr<IPressureSolver<MacGridState2D>> pressure;
    std::unique_ptr<MacGridFluidSolver2D> solver;
};

struct MacSimulation3D {
    SimulationConfig configuration;
    std::unique_ptr<IAdvection<MacAdvectionOperation>> advection;
    std::unique_ptr<IPressureSolver<MacGridState>> pressure;
    std::unique_ptr<MacGridFluidSolver> solver;
};

// Typed composition root for the MAC family. Each registry entry creates one
// module. A configuration is resolved domain -> advection -> projection ->
// fluid algorithm, then the modules are wired together.
class MacSimulationRegistry {
public:
    using Domain2DFactory = std::function<MacGridDomain2D(const SimulationConfig&)>;
    using Advection2DFactory = std::function<std::unique_ptr<IAdvection<MacAdvectionOperation2D>>()>;
    using Pressure2DFactory = std::function<std::unique_ptr<IPressureSolver<MacGridState2D>>()>;
    using Fluid2DFactory =
        std::function<std::unique_ptr<MacGridFluidSolver2D>(MacGridState2D, const SimulationConfig&)>;

    using Domain3DFactory = std::function<MacGridDomain(const SimulationConfig&)>;
    using Advection3DFactory = std::function<std::unique_ptr<IAdvection<MacAdvectionOperation>>()>;
    using Pressure3DFactory = std::function<std::unique_ptr<IPressureSolver<MacGridState>>()>;
    using Fluid3DFactory =
        std::function<std::unique_ptr<MacGridFluidSolver>(MacGridState, const SimulationConfig&)>;
    using Prepare3D = std::function<bool(const SimulationConfig&, const MacGridState&)>;

    MacSimulationRegistry();

    void registerDomain2D(DomainType type, Domain2DFactory factory);
    void registerAdvection2D(MacAdvectionKey key, Advection2DFactory factory);
    void registerPressure2D(MacProjectionKey key, Pressure2DFactory factory);
    void registerFluid2D(MacFluidKey key, Fluid2DFactory factory);

    void registerDomain3D(DomainType type, Domain3DFactory factory);
    void registerAdvection3D(MacAdvectionKey key, Advection3DFactory factory);
    void registerPressure3D(MacProjectionKey key, Pressure3DFactory factory);
    void registerFluid3D(MacFluidKey key, Fluid3DFactory factory);
    void registerPrepare3D(Prepare3D prepare);

    [[nodiscard]] MacSimulation2D create2D(const SimulationConfig& config) const;
    [[nodiscard]] MacSimulation3D create3D(const SimulationConfig& config) const;

private:
    struct EnumHash {
        template<class Enum>
        std::size_t operator()(Enum value) const noexcept {
            return static_cast<std::size_t>(value);
        }
    };
    struct AdvectionKeyHash {
        std::size_t operator()(const MacAdvectionKey& key) const noexcept;
    };
    struct ProjectionKeyHash {
        std::size_t operator()(const MacProjectionKey& key) const noexcept;
    };
    struct FluidKeyHash {
        std::size_t operator()(const MacFluidKey& key) const noexcept;
    };

    std::unordered_map<DomainType, Domain2DFactory, EnumHash> domains2D_;
    std::unordered_map<MacAdvectionKey, Advection2DFactory, AdvectionKeyHash> advections2D_;
    std::unordered_map<MacProjectionKey, Pressure2DFactory, ProjectionKeyHash> pressures2D_;
    std::unordered_map<MacFluidKey, Fluid2DFactory, FluidKeyHash> fluids2D_;

    std::unordered_map<DomainType, Domain3DFactory, EnumHash> domains3D_;
    std::unordered_map<MacAdvectionKey, Advection3DFactory, AdvectionKeyHash> advections3D_;
    std::unordered_map<MacProjectionKey, Pressure3DFactory, ProjectionKeyHash> pressures3D_;
    std::unordered_map<MacFluidKey, Fluid3DFactory, FluidKeyHash> fluids3D_;
    Prepare3D prepare3D_;
};
