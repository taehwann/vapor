#pragma once
#include "solvers/simplicial3d/SimplicialCompute3D.hpp"
#include <memory>
#include <string>
class GpuSimplicial3D final:public SimplicialCompute3D {
public:
    GpuSimplicial3D();
    ~GpuSimplicial3D();
    GpuSimplicial3D(const GpuSimplicial3D&)=delete;
    GpuSimplicial3D& operator=(const GpuSimplicial3D&)=delete;
    void initialize(const SimplicialFluidSolver3D&);
    void advectCirculation(const SimplicialFluidSolver3D&,double,std::vector<double>&) override;
    void advectDensity(const SimplicialFluidSolver3D&,double,std::vector<float>&) override;
    LinearSolveResult recover(const SimplicialFluidSolver3D&,std::span<const double>,std::vector<double>&) override;
    const std::string& device() const;
private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};
