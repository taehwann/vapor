#pragma once
// Four supported timestep/advection combinations. CPU and GPU share these choices.
enum class MacAlgorithm {
    SimpleSemiLagrangian,
    SimpleMacCormack,
    ReflectionSemiLagrangian,
    ReflectionMacCormack
};
enum class AdvectionMethod { SemiLagrangian, MacCormack };
