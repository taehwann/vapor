#include "fluid/MacGridSampler.hpp"
#include "fluid/advection/MacCormack.hpp"
#include "fluid/advection/SemiLagrangian.hpp"

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
    MacGridState grid(4, 4.f);
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
    MacGridSampler sampler(grid);
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
    MacGridState single(1);
    single.density()[0] = 2.f;
    near(MacGridSampler(single).sampleCell(single.density(),{1,1,1}), 2.f);
}

void stationary() {
    SemiLagrangian sl;
    MacCormack mc;
    const MacGridAdvection* methods[] = {&sl, &mc};
    for (int n : {1,4,5}) for (const auto* method : methods) {
        MacGridState grid(n,float(n)); // Zero tracer velocity.
        AdvectionWorkspace workspace;
        for (std::size_t i=0; i<grid.cellCount(); ++i) grid.density()[i] = float(i%7)/7.f;
        const std::vector<float> source(grid.density().begin(),grid.density().end());
        method->advectScalar({grid,velocityView(grid),.2f},grid.density(),workspace);
        same(grid.density(),source);
        std::vector<float> face(grid.velocityX().size()), output(face.size());
        for (std::size_t i=0; i<face.size(); ++i) face[i] = float(i%13)/13.f;
        for (auto axis : {FaceAxis::X,FaceAxis::Y,FaceAxis::Z}) {
            method->advectFace({grid,velocityView(grid),.2f},axis,face,output,workspace);
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
    MacGridState grid(6,6.f);
    std::fill(grid.velocityX().begin(),grid.velocityX().end(),.5f);
    SemiLagrangian sl;
    MacCormack mc;
    const MacGridAdvection* methods[] = {&sl,&mc};
    AdvectionWorkspace workspace;
    for (const auto* method : methods) {
        for (int z=0; z<6; ++z) for (int y=0; y<6; ++y) for (int x=0; x<6; ++x)
            grid.density()[grid.idC(x,y,z)] = float(x+2*y+3*z);
        method->advectScalar({grid,velocityView(grid),.5f},grid.density(),workspace);
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
            method->advectFace({grid,velocityView(grid),.5f},axis,source,output,workspace);
            near(output[index(2,2,2)],12.f-.25f);
        }
    }
}

void limiter() {
    MacGridState grid(6,6.f);
    std::fill(grid.velocityX().begin(),grid.velocityX().end(),.5f);
    SemiLagrangian sl;
    MacCormack mc;
    AdvectionWorkspace workspace;
    grid.density()[grid.idC(3,3,3)] = 1.f;
    std::vector<float> semi(grid.density().begin(),grid.density().end());
    sl.advectScalar({grid,velocityView(grid),.4f},semi,workspace);
    mc.advectScalar({grid,velocityView(grid),.4f},grid.density(),workspace);
    require(grid.density()[grid.idC(3,3,3)] > semi[grid.idC(3,3,3)]+.01f,
            "MacCormack must execute a correction, not just select a label");
    for (float value : grid.density()) require(value>=0.f && value<=1.f, "Scalar limiter bounds");
    for (auto axis : {FaceAxis::X,FaceAxis::Y,FaceAxis::Z}) {
        std::vector<float> source(grid.velocityX().size(),0.f), output(source.size()), reference(source.size());
        const auto i = axis==FaceAxis::X ? grid.idX(3,3,3) : axis==FaceAxis::Y ? grid.idY(3,3,3) : grid.idZ(3,3,3);
        source[i] = 1.f;
        mc.advectFace({grid,velocityView(grid),.4f},axis,source,output,workspace);
        for (float value : output) require(value>=0.f && value<=1.f, "Velocity limiter bounds");
        sl.advectFace({grid,velocityView(grid),.4f},axis,source,reference,workspace);
        mc.advectFace({grid,velocityView(grid),.4f,.001f},axis,source,output,workspace);
        same(output,reference); // Every nonzero trace exceeds the deliberately tiny CFL threshold.
    }
}

void contracts() {
    MacGridState grid(3);
    SemiLagrangian method;
    AdvectionWorkspace workspace;
    AdvectionContext context{grid,velocityView(grid),.1f};
    std::vector<float> source(grid.velocityX().size()), output(source.size());
    rejects([&]{method.advectScalar(context,{},workspace);});
    rejects([&]{method.advectFace(context,FaceAxis::X,source,source,workspace);});
    rejects([&]{method.advectFace(context,FaceAxis::X,source,grid.velocityX(),workspace);});
    rejects([&]{method.advectFace(context,static_cast<FaceAxis>(99),source,output,workspace);});
    rejects([&]{method.advectFace(context,FaceAxis::X,{},output,workspace);});
    context.dt = -1.f;
    rejects([&]{method.advectScalar(context,grid.density(),workspace);});
    context.dt = std::numeric_limits<float>::quiet_NaN();
    rejects([&]{method.advectScalar(context,grid.density(),workspace);});
    context.dt = .1f;
    context.flow.x = {};
    rejects([&]{method.advectScalar(context,grid.density(),workspace);});
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
