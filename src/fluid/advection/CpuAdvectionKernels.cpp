#include "fluid/advection/CpuAdvectionKernels.hpp"
#include "fluid/MacGridSampler.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <stdexcept>

namespace {
// These are the extracted numerical loops, including the existing in-place
// scalar limiter. Keep the loops serial until that limiter is made alias-safe.
class Kernel {
public:
    Kernel(const AdvectionContext& context, AdvectionWorkspace& workspace, bool correction)
        : state_(context.grid), sampler_(context.grid), flow_(context.flow),
          correction_(correction), fwdScratch(workspace.forward),
          backScratch(workspace.backward), cflMc(context.cflLimit) {}

    void advectScalar(std::span<float> field, float dt) {
        const float dx = h();
        const int N = state_.resolution() * state_.resolution() * state_.resolution();

        for (int z = 0; z < state_.resolution(); ++z) {
            for (int y = 0; y < state_.resolution(); ++y) {
                for (int x = 0; x < state_.resolution(); ++x) {
                    Vec3 p{ (x + .5f) * dx, (y + .5f) * dx, (z + .5f) * dx };
                    fwdScratch[idC(x, y, z)] = sampleCell(field, p - velocity(p) * dt);
                }
            }
        }

        if (!correction_) {
            std::copy(fwdScratch.begin(), fwdScratch.begin() + N, field.begin());
            return;
        }

        for (int z = 0; z < state_.resolution(); ++z) {
            for (int y = 0; y < state_.resolution(); ++y) {
                for (int x = 0; x < state_.resolution(); ++x) {
                    Vec3 p{ (x + .5f) * dx, (y + .5f) * dx, (z + .5f) * dx };
                    backScratch[idC(x, y, z)] = sampleCell(fwdScratch, p + velocity(p) * dt);
                }
            }
        }

        for (int z = 0; z < state_.resolution(); ++z) {
            for (int y = 0; y < state_.resolution(); ++y) {
                for (int x = 0; x < state_.resolution(); ++x) {
                    const int i = idC(x, y, z);
                    Vec3 p{ (x + .5f) * dx, (y + .5f) * dx, (z + .5f) * dx };
                    Vec3 v = velocity(p);
                    Vec3 posFwd = p - v * dt;
                    Vec3 posBwd = p + v * dt;
                    if (posFwd.x < 0 || posFwd.x > state_.boxSize() || posFwd.y < 0 || posFwd.y > state_.boxSize() || posFwd.z < 0 || posFwd.z > state_.boxSize() ||
                        posBwd.x < 0 || posBwd.x > state_.boxSize() || posBwd.y < 0 || posBwd.y > state_.boxSize() || posBwd.z < 0 || posBwd.z > state_.boxSize()) {
                        field[i] = fwdScratch[i];
                        continue;
                    }
                    float corrected = fwdScratch[i] + 0.5f * (field[i] - backScratch[i]);
                    int gx = std::clamp(int(posFwd.x / dx - 0.5f), 0, state_.resolution() - 1);
                    int gy = std::clamp(int(posFwd.y / dx - 0.5f), 0, state_.resolution() - 1);
                    int gz = std::clamp(int(posFwd.z / dx - 0.5f), 0, state_.resolution() - 1);
                    float fMin = field[idC(gx, gy, gz)], fMax = fMin;
                    for (int dz = 0; dz <= 1; ++dz)
                        for (int dy = 0; dy <= 1; ++dy)
                            for (int dx2 = 0; dx2 <= 1; ++dx2) {
                                float val = field[idC(std::clamp(gx + dx2, 0, state_.resolution() - 1),
                                    std::clamp(gy + dy, 0, state_.resolution() - 1),
                                    std::clamp(gz + dz, 0, state_.resolution() - 1))];
                                fMin = std::min(fMin, val);
                                fMax = std::max(fMax, val);
                            }
                    field[i] = std::clamp(corrected, fMin, fMax);
                }
            }
        }
    }

    void advectVx(float dt, std::span<const float> vxSrc,
        std::span<const float> vxVel, std::span<const float> vyVel, std::span<const float> vzVel,
        std::span<float> out, bool useMC = true) {
        const float dx = h();
        const int N = (int)out.size();

        for (int z = 0; z < state_.resolution(); ++z) {
            for (int y = 0; y < state_.resolution(); ++y) {
                for (int x = 0; x <= state_.resolution(); ++x) {
                    Vec3 face{ x * dx, (y + .5f) * dx, (z + .5f) * dx };
                    fwdScratch[idX(x, y, z)] = sampleVxField(vxSrc, face - velocity(face, vxVel, vyVel, vzVel) * dt);
                }
            }
        }

        if (!useMC) {
            std::copy(fwdScratch.begin(), fwdScratch.begin() + N, out.begin());
            return;
        }

        for (int z = 0; z < state_.resolution(); ++z) {
            for (int y = 0; y < state_.resolution(); ++y) {
                for (int x = 0; x <= state_.resolution(); ++x) {
                    Vec3 face{ x * dx, (y + .5f) * dx, (z + .5f) * dx };
                    backScratch[idX(x, y, z)] = sampleVxField(fwdScratch, face + velocity(face, vxVel, vyVel, vzVel) * dt);
                }
            }
        }

        for (int z = 0; z < state_.resolution(); ++z) {
            for (int y = 0; y < state_.resolution(); ++y) {
                for (int x = 0; x <= state_.resolution(); ++x) {
                    const int i = idX(x, y, z);
                    Vec3 face{ x * dx, (y + .5f) * dx, (z + .5f) * dx };
                    Vec3 v = velocity(face, vxVel, vyVel, vzVel);
                    Vec3 posFwd = face - v * dt, posBwd = face + v * dt;
                    float cfl = std::sqrt(v.x * v.x + v.y * v.y + v.z * v.z) * dt / dx;
                    if (posFwd.x < 0 || posFwd.x > state_.boxSize() || posFwd.y < 0 || posFwd.y > state_.boxSize() || posFwd.z < 0 || posFwd.z > state_.boxSize() ||
                        posBwd.x < 0 || posBwd.x > state_.boxSize() || posBwd.y < 0 || posBwd.y > state_.boxSize() || posBwd.z < 0 || posBwd.z > state_.boxSize() ||
                        cfl > cflMc) {
                        out[i] = fwdScratch[i];
                        continue;
                    }
                    float corrected = fwdScratch[i] + 0.5f * (vxSrc[i] - backScratch[i]);
                    int gx = std::clamp(int(posFwd.x / dx), 0, state_.resolution());
                    int gy = std::clamp(int(posFwd.y / dx - 0.5f), 0, state_.resolution() - 1);
                    int gz = std::clamp(int(posFwd.z / dx - 0.5f), 0, state_.resolution() - 1);
                    float fMin = vxSrc[idX(gx, gy, gz)], fMax = fMin;
                    for (int dz = 0; dz <= 1; ++dz)
                        for (int dy = 0; dy <= 1; ++dy)
                            for (int dx2 = 0; dx2 <= 1; ++dx2) {
                                float val = vxSrc[idX(std::clamp(gx + dx2, 0, state_.resolution()),
                                    std::clamp(gy + dy, 0, state_.resolution() - 1),
                                    std::clamp(gz + dz, 0, state_.resolution() - 1))];
                                fMin = std::min(fMin, val);
                                fMax = std::max(fMax, val);
                            }
                    out[i] = std::clamp(corrected, fMin, fMax);
                }
            }
        }
    }

    void advectVy(float dt, std::span<const float> vySrc,
        std::span<const float> vxVel, std::span<const float> vyVel, std::span<const float> vzVel,
        std::span<float> out, bool useMC = true) {
        const float dx = h();
        const int N = (int)out.size();

        for (int z = 0; z < state_.resolution(); ++z) {
            for (int y = 0; y <= state_.resolution(); ++y) {
                for (int x = 0; x < state_.resolution(); ++x) {
                    Vec3 face{ (x + .5f) * dx, y * dx, (z + .5f) * dx };
                    fwdScratch[idY(x, y, z)] = sampleVyField(vySrc, face - velocity(face, vxVel, vyVel, vzVel) * dt);
                }
            }
        }

        if (!useMC) {
            std::copy(fwdScratch.begin(), fwdScratch.begin() + N, out.begin());
            return;
        }

        for (int z = 0; z < state_.resolution(); ++z) {
            for (int y = 0; y <= state_.resolution(); ++y) {
                for (int x = 0; x < state_.resolution(); ++x) {
                    Vec3 face{ (x + .5f) * dx, y * dx, (z + .5f) * dx };
                    backScratch[idY(x, y, z)] = sampleVyField(fwdScratch, face + velocity(face, vxVel, vyVel, vzVel) * dt);
                }
            }
        }

        for (int z = 0; z < state_.resolution(); ++z) {
            for (int y = 0; y <= state_.resolution(); ++y) {
                for (int x = 0; x < state_.resolution(); ++x) {
                    const int i = idY(x, y, z);
                    Vec3 face{ (x + .5f) * dx, y * dx, (z + .5f) * dx };
                    Vec3 v = velocity(face, vxVel, vyVel, vzVel);
                    Vec3 posFwd = face - v * dt, posBwd = face + v * dt;
                    float cfl = std::sqrt(v.x * v.x + v.y * v.y + v.z * v.z) * dt / dx;
                    if (posFwd.x < 0 || posFwd.x > state_.boxSize() || posFwd.y < 0 || posFwd.y > state_.boxSize() || posFwd.z < 0 || posFwd.z > state_.boxSize() ||
                        posBwd.x < 0 || posBwd.x > state_.boxSize() || posBwd.y < 0 || posBwd.y > state_.boxSize() || posBwd.z < 0 || posBwd.z > state_.boxSize() ||
                        cfl > cflMc) {
                        out[i] = fwdScratch[i];
                        continue;
                    }
                    float corrected = fwdScratch[i] + 0.5f * (vySrc[i] - backScratch[i]);
                    int gx = std::clamp(int(posFwd.x / dx - 0.5f), 0, state_.resolution() - 1);
                    int gy = std::clamp(int(posFwd.y / dx), 0, state_.resolution());
                    int gz = std::clamp(int(posFwd.z / dx - 0.5f), 0, state_.resolution() - 1);
                    float fMin = vySrc[idY(gx, gy, gz)], fMax = fMin;
                    for (int dz = 0; dz <= 1; ++dz)
                        for (int dy = 0; dy <= 1; ++dy)
                            for (int dx2 = 0; dx2 <= 1; ++dx2) {
                                float val = vySrc[idY(std::clamp(gx + dx2, 0, state_.resolution() - 1),
                                    std::clamp(gy + dy, 0, state_.resolution()),
                                    std::clamp(gz + dz, 0, state_.resolution() - 1))];
                                fMin = std::min(fMin, val);
                                fMax = std::max(fMax, val);
                            }
                    out[i] = std::clamp(corrected, fMin, fMax);
                }
            }
        }
    }

    void advectVz(float dt, std::span<const float> vzSrc,
        std::span<const float> vxVel, std::span<const float> vyVel, std::span<const float> vzVel,
        std::span<float> out, bool useMC = true) {
        const float dx = h();
        const int N = (int)out.size();

        for (int z = 0; z <= state_.resolution(); ++z) {
            for (int y = 0; y < state_.resolution(); ++y) {
                for (int x = 0; x < state_.resolution(); ++x) {
                    Vec3 face{ (x + .5f) * dx, (y + .5f) * dx, z * dx };
                    fwdScratch[idZ(x, y, z)] = sampleVzField(vzSrc, face - velocity(face, vxVel, vyVel, vzVel) * dt);
                }
            }
        }

        if (!useMC) {
            std::copy(fwdScratch.begin(), fwdScratch.begin() + N, out.begin());
            return;
        }

        for (int z = 0; z <= state_.resolution(); ++z) {
            for (int y = 0; y < state_.resolution(); ++y) {
                for (int x = 0; x < state_.resolution(); ++x) {
                    Vec3 face{ (x + .5f) * dx, (y + .5f) * dx, z * dx };
                    backScratch[idZ(x, y, z)] = sampleVzField(fwdScratch, face + velocity(face, vxVel, vyVel, vzVel) * dt);
                }
            }
        }

        for (int z = 0; z <= state_.resolution(); ++z) {
            for (int y = 0; y < state_.resolution(); ++y) {
                for (int x = 0; x < state_.resolution(); ++x) {
                    const int i = idZ(x, y, z);
                    Vec3 face{ (x + .5f) * dx, (y + .5f) * dx, z * dx };
                    Vec3 v = velocity(face, vxVel, vyVel, vzVel);
                    Vec3 posFwd = face - v * dt, posBwd = face + v * dt;
                    float cfl = std::sqrt(v.x * v.x + v.y * v.y + v.z * v.z) * dt / dx;
                    if (posFwd.x < 0 || posFwd.x > state_.boxSize() || posFwd.y < 0 || posFwd.y > state_.boxSize() || posFwd.z < 0 || posFwd.z > state_.boxSize() ||
                        posBwd.x < 0 || posBwd.x > state_.boxSize() || posBwd.y < 0 || posBwd.y > state_.boxSize() || posBwd.z < 0 || posBwd.z > state_.boxSize() ||
                        cfl > cflMc) {
                        out[i] = fwdScratch[i];
                        continue;
                    }
                    float corrected = fwdScratch[i] + 0.5f * (vzSrc[i] - backScratch[i]);
                    int gx = std::clamp(int(posFwd.x / dx - 0.5f), 0, state_.resolution() - 1);
                    int gy = std::clamp(int(posFwd.y / dx - 0.5f), 0, state_.resolution() - 1);
                    int gz = std::clamp(int(posFwd.z / dx), 0, state_.resolution());
                    float fMin = vzSrc[idZ(gx, gy, gz)], fMax = fMin;
                    for (int dz = 0; dz <= 1; ++dz)
                        for (int dy = 0; dy <= 1; ++dy)
                            for (int dx2 = 0; dx2 <= 1; ++dx2) {
                                float val = vzSrc[idZ(std::clamp(gx + dx2, 0, state_.resolution() - 1),
                                    std::clamp(gy + dy, 0, state_.resolution() - 1),
                                    std::clamp(gz + dz, 0, state_.resolution()))];
                                fMin = std::min(fMin, val);
                                fMax = std::max(fMax, val);
                            }
                    out[i] = std::clamp(corrected, fMin, fMax);
                }
            }
        }
    }

private:
    const MacGridState& state_;
    MacGridSampler sampler_;
    VelocityFieldView flow_;
    bool correction_;
    std::vector<float>& fwdScratch;
    std::vector<float>& backScratch;
    float cflMc;

    inline int idC(int x, int y, int z) const { return static_cast<int>(state_.idC(x, y, z)); }
    inline int idX(int x, int y, int z) const { return static_cast<int>(state_.idX(x, y, z)); }
    inline int idY(int x, int y, int z) const { return static_cast<int>(state_.idY(x, y, z)); }
    inline int idZ(int x, int y, int z) const { return static_cast<int>(state_.idZ(x, y, z)); }
    inline float h() const { return state_.cellSize(); }
    float sampleCell(std::span<const float> field, Vec3 p) const { return sampler_.sampleCell(field, p); }
    float sampleVxField(std::span<const float> field, Vec3 p) const { return sampler_.sampleVxField(field, p); }
    float sampleVyField(std::span<const float> field, Vec3 p) const { return sampler_.sampleVyField(field, p); }
    float sampleVzField(std::span<const float> field, Vec3 p) const { return sampler_.sampleVzField(field, p); }
    Vec3 velocity(Vec3 p) const { return sampler_.velocity(p, flow_); }
    Vec3 velocity(Vec3 p, std::span<const float> x, std::span<const float> y, std::span<const float> z) const {
        return sampler_.velocity(p, {x, y, z});
    }
};
}

void CpuAdvectionKernels::advectScalar(const AdvectionContext& context, std::span<float> field,
                                      AdvectionWorkspace& workspace, bool correction) {
    validateScalarAdvection(context, field);
    workspace.resize(context.grid);
    Kernel(context, workspace, correction).advectScalar(field, context.dt);
}

void CpuAdvectionKernels::advectFace(const AdvectionContext& context, FaceAxis axis,
                                    std::span<const float> source, std::span<float> output,
                                    AdvectionWorkspace& workspace, bool correction) {
    validateFaceAdvection(context, axis, source, output);
    workspace.resize(context.grid);
    Kernel kernel(context, workspace, correction);
    const auto flow = context.flow;
    switch (axis) {
        case FaceAxis::X: kernel.advectVx(context.dt, source, flow.x, flow.y, flow.z, output, correction); break;
        case FaceAxis::Y: kernel.advectVy(context.dt, source, flow.x, flow.y, flow.z, output, correction); break;
        case FaceAxis::Z: kernel.advectVz(context.dt, source, flow.x, flow.y, flow.z, output, correction); break;
    }
}
