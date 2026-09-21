#pragma once
#include "solvers/mac3d/MacGridOperators3D.hpp"
#include "solvers/mac3d/MacGridSampler3D.hpp"
#include <cmath>
#include <filesystem>
#include <fstream>

// Diagnostic experiment, deliberately separate from pass/fail regression tests.
inline void compareReflection(int n, int frames, int iterations) {
    std::filesystem::create_directories("build/reflection-comparison");
    for (int variant = 0; variant < 4; ++variant) {
        MacGridParameters3D p;
        p.reflection = true;
        p.macCormackVel = (variant & 1) != 0;
        p.macCormackSmoke = (variant & 2) != 0;
        p.emitterRadius = .04f;
        p.sourceStrength = 1;
        p.buoyancy = .8f;
        p.smokeDecay = .3007525f;
        p.projectIterations = iterations;
        p.sorOmega = 1.9f;
        GpuMacSimulation3D simulation(n, 1, p);
        VolumeRenderer renderer;
        renderer.init(simulation.renderData());
        for (int frame = 1; frame <= frames; ++frame) {
            simulation.advance(1.f / 60);
            if (frame % 60 != 0 && frame != frames)
                continue;
            MacGridState3D s(n, 1);
            simulation.download(s);
            MacGridSampler3D sampler(s);
            double mass = 0, grad = 0, square = 0, divSquare = 0, maxCfl = 0;
            size_t smokeCells = 0, cflFallback = 0, boundaryFallback = 0;
            std::vector<float> div(s.cellCount());
            MacGridOperators3D::computeDivergence(s, div);
            for (float d : div)
                divSquare += double(d) * d;
            for (int z = 0; z < n; ++z)
                for (int y = 0; y < n; ++y)
                    for (int x = 0; x < n; ++x) {
                        float d = s.density()[s.idC(x, y, z)];
                        mass += d;
                        square += double(d) * d;
                        if (x + 1 < n)
                            grad += std::abs(d - s.density()[s.idC(x + 1, y, z)]);
                        if (y + 1 < n)
                            grad += std::abs(d - s.density()[s.idC(x, y + 1, z)]);
                        if (z + 1 < n)
                            grad += std::abs(d - s.density()[s.idC(x, y, z + 1)]);
                        if (d <= .001f)
                            continue;
                        ++smokeCells;
                        Vec3 pos{(x + .5f) / n, (y + .5f) / n, (z + .5f) / n};
                        Vec3 v = sampler.velocity(pos, velocityView(s));
                        double cfl = std::sqrt(v.x * v.x + v.y * v.y + v.z * v.z) * n / 60;
                        maxCfl = std::max(maxCfl, cfl);
                        if (cfl > p.cflMc)
                            ++cflFallback;
                        Vec3 a = pos - v * (1.f / 60), b = pos + v * (1.f / 60);
                        if (a.x < 0 || a.x > 1 || a.y < 0 || a.y > 1 || a.z < 0 || a.z > 1 || b.x < 0 ||
                            b.x > 1 || b.y < 0 || b.y > 1 || b.z < 0 || b.z > 1)
                            ++boundaryFallback;
                    }
            std::cout << "compare variant=" << variant << " mcV=" << p.macCormackVel
                      << " mcD=" << p.macCormackSmoke << " frame=" << frame << " mass=" << mass / (n * n * n)
                      << " grad_per_mass=" << grad / mass << " square_per_mass=" << square / mass
                      << " rms_div=" << std::sqrt(divSquare / div.size())
                      << " max_div=" << MacGridOperators3D::maxDivergence(s) << " smoke_cells=" << smokeCells
                      << " smoke_cfl_fallback=" << cflFallback
                      << " smoke_boundary_fallback=" << boundaryFallback << " smoke_max_cfl=" << maxCfl
                      << std::endl;
            std::string base =
                "build/reflection-comparison/v" + std::to_string(variant) + "_f" + std::to_string(frame);
            std::ofstream raw(base + ".density", std::ios::binary);
            raw.write(reinterpret_cast<const char *>(s.density().data()), s.density().size_bytes());
            float mvp[16], cx, cy, cz;
            buildMVP(mvp, cx, cy, cz, -1.2f, .3f, 2.5f, 1, 1);
            glClearColor(.03f, .04f, .07f, 1);
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
            const float norm = std::sqrt(1.25f);
            renderer.render(simulation.renderData(), 960, 960, mvp, cx, cy, cz, .5f, 30, .3f / norm,
                            .4f / norm, 1 / norm, .3f, .1f, nullptr);
            std::vector<unsigned char> pixels(960 * 960 * 3);
            glReadPixels(0, 0, 960, 960, GL_RGB, GL_UNSIGNED_BYTE, pixels.data());
            std::ofstream image(base + ".ppm", std::ios::binary);
            image << "P6\n960 960\n255\n";
            for (int row = 959; row >= 0; --row)
                image.write(reinterpret_cast<char *>(pixels.data() + row * 960 * 3), 960 * 3);
        }
    }
}
