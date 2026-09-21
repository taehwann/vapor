#include "graphics/GlLoader.hpp"
#include "renderer/MacRenderData.hpp"
#include "renderer/OrbitCamera.hpp"
#include "renderer/VolumeRenderer.hpp"
#include "solvers/mac3d/MacGridFluidSolver3D.hpp"
#include "solvers/mac3d/gpu/GpuAdvection.hpp"
#include "solvers/mac3d/gpu/GpuMacSimulation3D.hpp"
#include "solvers/mac3d/gpu/GpuPressureSolver.hpp"
#include <algorithm>
#include <chrono>
#include <iostream>
#include "ReflectionComparison.hpp"
using Clock = std::chrono::steady_clock;
template <class F> double measure(F &&f) {
    glFinish();
    auto start = Clock::now();
    f();
    glFinish();
    return std::chrono::duration<double, std::milli>(Clock::now() - start).count();
}
struct TimedAdvection final : Advection3D {
    GpuAdvection impl;
    double &ms;
    TimedAdvection(std::shared_ptr<GpuMacBackend> b, double &m) : impl(b), ms(m) {}
    void advectScalar(const AdvectionContext &c, std::span<float> f, AdvectionWorkspace &w,
                      AdvectionMethod m) const override {
        ms += measure([&] { impl.advectScalar(c, f, w, m); });
    }
    void advectFace(const AdvectionContext &c, FaceAxis a, std::span<const float> s, std::span<float> o,
                    AdvectionWorkspace &w, AdvectionMethod m) const override {
        ms += measure([&] { impl.advectFace(c, a, s, o, w, m); });
    }
};
struct TimedPressure final : Projection3D {
    GpuPressureSolver impl;
    double &ms;
    TimedPressure(std::shared_ptr<GpuMacBackend> b, double &m) : impl(b), ms(m) {}
    ProjectionResult project(MacGridState3D &s, const ProjectionOptions &o) override {
        ProjectionResult r;
        ms += measure([&] { r = impl.project(s, o); });
        return r;
    }
    void reset() noexcept override { impl.reset(); }
};
int main(int argc, char **argv) {
    if (!glfwInit())
        return 77;
    glfwWindowHint(GLFW_VISIBLE, GLFW_FALSE);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    auto window = glfwCreateWindow(960, 960, "GPU benchmark", nullptr, nullptr);
    if (!window) {
        glfwTerminate();
        return 77;
    }
    glfwMakeContextCurrent(window);
    loadOpenGLFunctions();
    int status = 0;
    try {
        std::cout << "GPU: " << glGetString(GL_RENDERER) << std::endl;
        const int n = argc > 1 ? std::stoi(argv[1]) : 128, frames = argc > 2 ? std::stoi(argv[2]) : 6;
        if (argc > 3 && std::string(argv[3]) == "compare") {
            compareReflection(n, frames, argc > 4 ? std::stoi(argv[4]) : 120);
        } else if (argc > 3 && std::string(argv[3]) == "resident") {
            MacGridParameters3D p;
            p.reflection = false;
            p.macCormackSmoke = p.macCormackVel = false;
            p.sourceStrength = 1;
            p.buoyancy = .8f;
            p.smokeDecay = .3007525f;
            p.emitterRadius = .04f;
            p.projectIterations = 60;
            p.sorOmega = 1.9f;
            GpuMacSimulation3D s(n, 1, p);
            VolumeRenderer renderer;
            renderer.init(s.renderData());
            float mvp[16], cx, cy, cz;
            buildMVP(mvp, cx, cy, cz, -1.2f, .3f, 2.5f, 1, 1);
            double step = 0, draw = 0;
            std::uint64_t uploads = 0;
            for (int f = -2; f < frames; ++f) {
                if (f == 0)
                    uploads = s.uploadedBytes();
                double t = measure([&] { s.advance(1.f / 60); });
                double r = measure([&] {
                    renderer.render(s.renderData(), 960, 960, mvp, cx, cy, cz, .5f, 30, .3f, .4f, 1, .3f, .1f,
                                    nullptr);
                });
                if (f >= 0) {
                    step += t;
                    draw += r;
                }
            }
            std::cout << "resident n=" << n << " frames=" << frames << " step_ms=" << step / frames
                      << " render_ms=" << draw / frames
                      << " upload_bytes_per_step=" << (s.uploadedBytes() - uploads) / frames
                      << " total_download_bytes=" << s.downloadedBytes() << std::endl;
        } else {
            MacGridFluidSolver3D s(MacGridState3D(n, 1));
            auto &p = s.parameters();
            p.reflection = false;
            p.macCormackSmoke = p.macCormackVel = false;
            p.sourceStrength = 1;
            p.buoyancy = .8f;
            p.smokeDecay = .3007525f;
            p.emitterRadius = .04f;
            p.projectIterations = 60;
            p.sorOmega = 1.9f;
            auto backend = std::make_shared<GpuMacBackend>();
            backend->init(s.state());
            double advection = 0, pressure = 0;
            s.setAdvection(std::make_unique<TimedAdvection>(backend, advection));
            s.setPressureSolver(std::make_unique<TimedPressure>(backend, pressure));
            VolumeRenderer renderer;
            renderer.init(volumeData(s.state()));
            float mvp[16], cx, cy, cz;
            buildMVP(mvp, cx, cy, cz, -1.2f, .3f, 2.5f, 1, 1);
            double step = 0, draw = 0, stats = 0;
            for (int f = -2; f < frames; ++f) {
                if (f == 0)
                    advection = pressure = 0;
                double t = measure([&] { s.advance(1.f / 60); });
                double r = measure([&] {
                    renderer.render(volumeData(s.state()), 960, 960, mvp, cx, cy, cz, .5f, 30, .3f, .4f, 1,
                                    .3f, .1f, nullptr);
                });
                double d = measure([&] {
                    auto a = s.diagnostics();
                    volatile float keep = a.kineticEnergy + *std::max_element(s.state().density().begin(),
                                                                              s.state().density().end());
                    (void)keep;
                });
                if (f >= 0) {
                    step += t;
                    draw += r;
                    stats += d;
                }
            }
            std::cout << "baseline n=" << n << " frames=" << frames << " step_ms=" << step / frames
                      << " advection_including_transfers_ms=" << advection / frames
                      << " pressure_including_transfers_and_cpu_divergence_ms=" << pressure / frames
                      << " other_cpu_ms=" << (step - advection - pressure) / frames
                      << " render_upload_ms=" << draw / frames << " diagnostics_ms=" << stats / frames
                      << std::endl;
        }
    } catch (const std::exception &e) {
        std::cerr << e.what() << std::endl;
        status = 1;
    }
    glfwDestroyWindow(window);
    glfwTerminate();
    return status;
}
