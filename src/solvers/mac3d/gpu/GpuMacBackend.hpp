#pragma once

class MacGridState3D;

// Shared OpenGL resource/dispatch layer used by GPU pressure and advection.
// init(), use, and destruction of live resources require a current context.
class GpuMacBackend {
public:
    GpuMacBackend() = default;
    ~GpuMacBackend();
    GpuMacBackend(const GpuMacBackend&) = delete;
    GpuMacBackend& operator=(const GpuMacBackend&) = delete;

    void init(const MacGridState3D& state);
    [[nodiscard]] bool initialized() const noexcept { return initialized_; }
    void validateGrid(const MacGridState3D& state) const;

    void project(float dt, int iterations, float omega, float hInv, float h,
        const float* vxData, const float* vyData, const float* vzData,
        float* vxOut, float* vyOut, float* vzOut);

    void advectSmokeGPU(float dt, float* smokeData, const float* vxVel, const float* vyVel, const float* vzVel, float boxSize, float cflMc, bool useMC);

    void advectVxGPU(float dt, const float* vxSrc, int srcSz, const float* vxVel, const float* vyVel, const float* vzVel, float* vxOut, float boxSize, float cflMc, bool useMC);

    void advectVyGPU(float dt, const float* vySrc, int srcSz, const float* vxVel, const float* vyVel, const float* vzVel, float* vyOut, float boxSize, float cflMc, bool useMC);

    void advectVzGPU(float dt, const float* vzSrc, int srcSz, const float* vxVel, const float* vyVel, const float* vzVel, float* vzOut, float boxSize, float cflMc, bool useMC);

    void shutdown();

private:
    friend class GpuMacSimulation3D;
    void projectDevice(float dt, int iterations, float omega, float hInv, float h);
    int vxSz = 0, vySz = 0, vzSz = 0;
    bool initialized_ = false;
    unsigned int progDiv = 0, progRBGS = 0, progCorrect = 0;
    unsigned int progAdvectSL = 0, progAdvectMC = 0, progCopy = 0;
    unsigned int ssboP = 0, ssboDiv = 0, ssboVx = 0, ssboVy = 0, ssboVz = 0;
    unsigned int ssboSmoke = 0, ssboFwd = 0, ssboAdvectSrc = 0, ssboAdvectOut = 0;
    int n = 0, maxSz = 0;

    void uploadVelocity(const float* vxData, const float* vyData, const float* vzData);
    void dispatchSL(int mode, int count, float dx, float dt,
                    unsigned int source = 0, unsigned int output = 0);
    void dispatchMC(int mode, float dx, float dt, float boxSize, float cflMc, unsigned int srcBo, unsigned int outBo);
    void copySSBO(unsigned int dstBo, unsigned int srcBo, int count);
};
