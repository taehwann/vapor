#include "solvers/mac3d/MacGridSampler3D.hpp"
#include "solvers/mac3d/CpuAdvection.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string_view>
#include <vector>

namespace {
void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}
void near(float actual, float expected, float tolerance = 1e-5f) {
    require(std::isfinite(actual) && std::abs(actual - expected) <= tolerance, "Unexpected sampled/advected value");
}
void same(std::span<const float> a, std::span<const float> b) {
    require(a.size() == b.size(), "Field sizes differ");
    for (std::size_t i = 0; i < a.size(); ++i) near(a[i], b[i]);
}
template<class F> void rejects(F function) {
    try { function(); }
    catch (const std::invalid_argument&) { return; }
    throw std::runtime_error("Expected rejection of invalid advection input");
}

void sampling() {
    MacGridState3D grid(4, 4.f);
    const int n = grid.resolution();
    auto value = [](float x, float y, float z) { return x + 2*y + 3*z; };
    for (int z=0; z<n; ++z) for (int y=0; y<n; ++y) for (int x=0; x<n; ++x)
        grid.density()[grid.idC(x,y,z)] = value(x+.5f,y+.5f,z+.5f);
    for (int z=0; z<n; ++z) for (int y=0; y<n; ++y) for (int x=0; x<=n; ++x)
        grid.velocityX()[grid.idX(x,y,z)] = value(float(x),y+.5f,z+.5f);
    for (int z=0; z<n; ++z) for (int y=0; y<=n; ++y) for (int x=0; x<n; ++x)
        grid.velocityY()[grid.idY(x,y,z)] = value(x+.5f,float(y),z+.5f);
    for (int z=0; z<=n; ++z) for (int y=0; y<n; ++y) for (int x=0; x<n; ++x)
        grid.velocityZ()[grid.idZ(x,y,z)] = value(x+.5f,y+.5f,float(z));
    MacGridSampler3D sampler(grid);
    const Vec3 p{1.2f,1.6f,2.1f};
    const float expected = value(p.x,p.y,p.z);
    near(sampler.sampleCell(grid.density(),p), expected);
    const auto velocity = sampler.velocity(p,velocityView(grid));
    near(velocity.x, expected); near(velocity.y, expected); near(velocity.z, expected);
    near(sampler.sampleCell(grid.density(),{-1,-1,-1}), 3.f);
    near(sampler.sampleVxField(grid.velocityX(),{-1,-1,-1}), 2.5f);
    near(sampler.sampleCell(grid.density(),{9,9,9}), 21.f);
    near(sampler.sampleVxField(grid.velocityX(),{9,9,9}), 21.5f);
    near(sampler.sampleVyField(grid.velocityY(),{9,9,9}), 22.f);
    near(sampler.sampleVzField(grid.velocityZ(),{9,9,9}), 22.5f);
    MacGridState3D single(1);
    single.density()[0] = 2.f;
    near(MacGridSampler3D(single).sampleCell(single.density(),{1,1,1}), 2.f);
}

void stationary() {
    CpuAdvection sl;
    CpuAdvection mc;
    const AdvectionMethod methods[] = {AdvectionMethod::SemiLagrangian, AdvectionMethod::MacCormack};
    CpuAdvection cpu;
    for (int n : {1,4,5}) for (auto method : methods) {
        MacGridState3D grid(n,float(n)); // Zero tracer velocity.
        AdvectionWorkspace workspace;
        for (std::size_t i=0; i<grid.cellCount(); ++i) grid.density()[i] = float(i%7)/7.f;
        const std::vector<float> source(grid.density().begin(),grid.density().end());
        cpu.advectScalar({grid,velocityView(grid),.2f},grid.density(),workspace, method);
        same(grid.density(),source);
        std::vector<float> face(grid.velocityX().size()), output(face.size());
        for (std::size_t i=0; i<face.size(); ++i) face[i] = float(i%13)/13.f;
        for (auto axis : {FaceAxis::X,FaceAxis::Y,FaceAxis::Z}) {
            cpu.advectFace({grid,velocityView(grid),.2f},axis,face,output,workspace, method);
            same(output,face);
        }
        workspace.reset();
        require(std::all_of(workspace.forward.begin(),workspace.forward.end(),[](float v){return v==0.f;}),
                "Workspace reset must clear forward values");
        require(std::all_of(workspace.backward.begin(),workspace.backward.end(),[](float v){return v==0.f;}),
                "Workspace reset must clear backward values");
    }
}

void transport() {
    MacGridState3D grid(6,6.f);
    std::fill(grid.velocityX().begin(),grid.velocityX().end(),.5f);
    CpuAdvection sl;
    CpuAdvection mc;
    const AdvectionMethod methods[] = {AdvectionMethod::SemiLagrangian, AdvectionMethod::MacCormack};
    CpuAdvection cpu;
    AdvectionWorkspace workspace;
    for (auto method : methods) {
        for (int z=0; z<6; ++z) for (int y=0; y<6; ++y) for (int x=0; x<6; ++x)
            grid.density()[grid.idC(x,y,z)] = float(x+2*y+3*z);
        cpu.advectScalar({grid,velocityView(grid),.5f},grid.density(),workspace, method);
        near(grid.density()[grid.idC(2,2,2)], 12.f-.25f);
        for (auto axis : {FaceAxis::X,FaceAxis::Y,FaceAxis::Z}) {
            std::vector<float> source(grid.velocityX().size()), output(source.size());
            auto index = [&](int x,int y,int z) {
                return axis==FaceAxis::X ? grid.idX(x,y,z) : axis==FaceAxis::Y ? grid.idY(x,y,z) : grid.idZ(x,y,z);
            };
            for (int z=0; z<6+(axis==FaceAxis::Z); ++z)
                for (int y=0; y<6+(axis==FaceAxis::Y); ++y)
                    for (int x=0; x<6+(axis==FaceAxis::X); ++x)
                        source[index(x,y,z)] = float(x+2*y+3*z);
            cpu.advectFace({grid,velocityView(grid),.5f},axis,source,output,workspace, method);
            near(output[index(2,2,2)],12.f-.25f);
        }
    }
}

void limiter() {
    // Reflection must not change the answer by reversing the CPU traversal.
    {
        MacGridState3D right(10, 10.f), left(10, 10.f);
        std::fill(right.velocityX().begin(), right.velocityX().end(), .7f);
        std::fill(left.velocityX().begin(), left.velocityX().end(), -.7f);
        for (int z = 0; z < 10; ++z) for (int y = 0; y < 10; ++y) for (int x = 0; x < 10; ++x) {
            const float value = float((x * 17 + y * 7 + z * 3) % 13) / 12.f;
            right.density()[right.idC(x,y,z)] = value;
            left.density()[left.idC(9-x,y,z)] = value;
        }
        CpuAdvection method;
        AdvectionWorkspace scratch;
        auto fallback = std::vector<float>(right.density().begin(), right.density().end());
        auto reference = fallback;
        method.advectScalar({right,velocityView(right),.6f,.001f},fallback,scratch,AdvectionMethod::MacCormack);
        method.advectScalar({right,velocityView(right),.6f},reference,scratch,AdvectionMethod::SemiLagrangian);
        same(fallback, reference);
        method.advectScalar({right,velocityView(right),.6f},right.density(),scratch,AdvectionMethod::MacCormack);
        method.advectScalar({left,velocityView(left),.6f},left.density(),scratch,AdvectionMethod::MacCormack);
        for (int z = 2; z < 8; ++z) for (int y = 2; y < 8; ++y) for (int x = 2; x < 8; ++x)
            near(right.density()[right.idC(x,y,z)], left.density()[left.idC(9-x,y,z)], 2e-6f);
    }
    MacGridState3D grid(6,6.f);
    std::fill(grid.velocityX().begin(),grid.velocityX().end(),.5f);
    CpuAdvection sl;
    CpuAdvection mc;
    AdvectionWorkspace workspace;
    grid.density()[grid.idC(3,3,3)] = 1.f;
    std::vector<float> semi(grid.density().begin(),grid.density().end());
    sl.advectScalar({grid,velocityView(grid),.4f},semi,workspace, AdvectionMethod::SemiLagrangian);
    mc.advectScalar({grid,velocityView(grid),.4f},grid.density(),workspace, AdvectionMethod::MacCormack);
    require(grid.density()[grid.idC(3,3,3)] > semi[grid.idC(3,3,3)]+.01f,
            "MacCormack must execute a correction, not just select a label");
    for (float value : grid.density()) require(value>=0.f && value<=1.f, "Scalar limiter bounds");
    for (auto axis : {FaceAxis::X,FaceAxis::Y,FaceAxis::Z}) {
        std::vector<float> source(grid.velocityX().size(),0.f), output(source.size()), reference(source.size());
        const auto i = axis==FaceAxis::X ? grid.idX(3,3,3) : axis==FaceAxis::Y ? grid.idY(3,3,3) : grid.idZ(3,3,3);
        source[i] = 1.f;
        mc.advectFace({grid,velocityView(grid),.4f},axis,source,output,workspace, AdvectionMethod::MacCormack);
        for (float value : output) require(value>=0.f && value<=1.f, "Velocity limiter bounds");
        sl.advectFace({grid,velocityView(grid),.4f},axis,source,reference,workspace, AdvectionMethod::SemiLagrangian);
        mc.advectFace({grid,velocityView(grid),.4f,.001f},axis,source,output,workspace, AdvectionMethod::MacCormack);
        same(output,reference); // Every nonzero trace exceeds the deliberately tiny CFL threshold.
    }
}

void contracts() {
    MacGridState3D grid(3);
    CpuAdvection method;
    AdvectionWorkspace workspace;
    AdvectionContext context{grid,velocityView(grid),.1f};
    std::vector<float> source(grid.velocityX().size()), output(source.size());
    rejects([&]{method.advectScalar(context,{},workspace, AdvectionMethod::SemiLagrangian);});
    rejects([&]{method.advectFace(context,FaceAxis::X,source,source,workspace, AdvectionMethod::SemiLagrangian);});
    rejects([&]{method.advectFace(context,FaceAxis::X,source,grid.velocityX(),workspace, AdvectionMethod::SemiLagrangian);});
    rejects([&]{method.advectFace(context,static_cast<FaceAxis>(99),source,output,workspace, AdvectionMethod::SemiLagrangian);});
    rejects([&]{method.advectFace(context,FaceAxis::X,{},output,workspace, AdvectionMethod::SemiLagrangian);});
    context.dt = -1.f;
    rejects([&]{method.advectScalar(context,grid.density(),workspace, AdvectionMethod::SemiLagrangian);});
    context.dt = std::numeric_limits<float>::quiet_NaN();
    rejects([&]{method.advectScalar(context,grid.density(),workspace, AdvectionMethod::SemiLagrangian);});
    context.dt = .1f;
    context.flow.x = {};
    rejects([&]{method.advectScalar(context,grid.density(),workspace, AdvectionMethod::SemiLagrangian);});
}

struct Test { const char* name; void (*run)(); };
const Test tests[] = {{"sampling",sampling},{"advection_stationary",stationary},
    {"advection_transport",transport},{"advection_limiter",limiter},{"advection_contracts",contracts}};
}

int main(int argc,char** argv) {
    int passed=0, failed=0;
    for (const auto& test : tests) {
        if (argc>1 && std::string_view(argv[1])!=test.name) continue;
        try { test.run(); std::cout << "PASS " << test.name << '\n'; ++passed; }
        catch (const std::exception& error) {
            std::cerr << "FAIL " << test.name << ": " << error.what() << '\n'; ++failed;
        }
    }
    return failed==0 && passed>0 ? 0 : 1;
}
