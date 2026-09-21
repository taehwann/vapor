#pragma once
#include "common/LinearSolveOptions.hpp"
#include <span>
#include <vector>
class SimplicialFluidSolver3D;
// Optional compute backend; the CPU implementation remains the reference.
class SimplicialCompute3D {
public:
    virtual ~SimplicialCompute3D()=default;
    virtual void advectCirculation(const SimplicialFluidSolver3D&,double,std::vector<double>&)=0;
    virtual void advectDensity(const SimplicialFluidSolver3D&,double,std::vector<float>&)=0;
    virtual LinearSolveResult recover(const SimplicialFluidSolver3D&,std::span<const double>,std::vector<double>&)=0;
};
