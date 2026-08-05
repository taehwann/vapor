# Vapor

An interactive **3D smoke-plume** simulation written in C++, rendered with OpenGL
volume ray-marching and Dear ImGui. The solver integrates the incompressible
Navier–Stokes equations on a 32×32×32 staggered (MAC) grid using semi-Lagrangian
advection, a midpoint-reflect (second-order) velocity update, and a red-black SOR
pressure solve.

---

## Grid layout

The domain is `[0, boxSize]³` (`boxSize` adjustable 0.5–5.0, default 3.0) divided
into `n×n×n` cells (`n = 32`, cell size `h = boxSize / n`). Quantities live on a
**staggered (MAC) grid**:

| Field       | Location           | Size               |
| ----------- | ------------------ | ------------------ |
| `smoke`     | cell centers       | `n³`               |
| `pressure`  | cell centers       | `n³`               |
| `divergence`| cell centers       | `n³`               |
| `vx`        | x-faces            | `(n+1) × n × n`    |
| `vy`        | y-faces            | `n × (n+1) × n`    |
| `vz`        | z-faces            | `n × n × (n+1)`    |

## Boundaries

Each of the 6 faces can be open or closed. Defaults: all open except bottom (-Z)
is a **solid wall**. Closed walls zero the face velocity. Open boundaries clamp
face velocities to outflow-only after pressure projection, preventing inflow.

The SOR pressure solve uses a one-sided stencil at boundaries — missing neighbors
are skipped, equivalent to a homogeneous Neumann condition (`∂p/∂n = 0`).

## Emitter

A spherical source at the **center of the box** (`emitterCenter = (0.5, 0.5, 0.5)`
in normalized coordinates) with adjustable radius. Emits constant density `1.0`
inside the sphere, with an upward velocity kick and random stirring jitter. All
emitter coordinates are in local [0,1] space and scale proportionally with
`boxSize`.

## Advection

**Semi-Lagrangian backtrace** with optional **MacCormack correction** (toggleable
separately for smoke and velocity). The MacCormack limiter uses a 3×3×3
neighborhood clamp to prevent overshoot/undershoot. Rays that leave the domain
fall back to the pure semi-Lagrangian result. A CFL-based guard (`cflMc`) also
disables correction for cells that would require excessive extrapolation.

## Time integration

Two-stage **midpoint-reflect** (Richardson extrapolation, second-order) or a
simpler first-order single-pass method (toggleable via "Reflection" checkbox).

1. Advect smoke by `halfDt`.
2. Predict velocity by `halfDt`, add buoyancy + emitter, project → midpoint `u₁/₂`.
3. Reflect: `û = 2·u₁/₂ − ũ₁/₂`.
4. Advect smoke and reflected velocity by `halfDt` with midpoint field, project.

Buoyancy is split across the two `vz` faces of each cell.

## Box resizing

When `boxSize` changes via the UI slider, all simulation fields (smoke, velocity,
pressure, divergence) are **cleared to zero** — a fresh start at the new domain
size. The emitter stays at its relative [0,1] position and all other parameters
are preserved.

## Volume rendering

3D texture (`GL_R16F`) uploaded to GPU, rendered via ray-marching through a unit
cube proxy. The fragment shader casts rays from the camera through each pixel,
intersects the domain box `[0, boxSize]³`, and accumulates smoke density by
alpha-compositing. Non-smoke rays are discarded, leaving the clear-color
background. Mouse drag orbits the camera; scroll wheel zooms.

## Build

Requires CMake ≥ 3.20, a C++20 compiler, and OpenGL 3.3. GLFW and Dear ImGui are
pulled in as git submodules.

```powershell
git submodule update --init --recursive
cmake -S . -B build
cmake --build build --config Release
.\build\Release\vapor.exe
```

## UI controls

| Control         | Description                                    |
| --------------- | ---------------------------------------------- |
| Pause / Resume  | Freeze or advance the simulation               |
| Reset           | Clear all fields, keep parameters              |
| Reflection      | Toggle midpoint-reflect vs single-pass         |
| MC smoke / vel  | Toggle MacCormack correction per field         |
| Buoyancy        | Upward force strength (0–10)                   |
| Source          | Emitter injection strength                     |
| Box size        | Domain extent [0, boxSize]³ (0.5–5.0)         |
| Emitt radius    | Emitter sphere radius (0.05–0.5)               |
| Stir            | Random turbulence in emitter (0–2)             |
| Alpha mul       | Volume rendering opacity multiplier            |
| SOR its         | Pressure solve iterations (10–500)             |
| MC CFL          | MacCormack CFL guard threshold                 |
| Step            | Ray-marching step scale                        |
| FPS             | Current frame rate                             |

## Files

| File          | Purpose                                              |
| ------------- | ---------------------------------------------------- |
| `main.cpp`    | Simulation (`SmokeSim3D`), volume renderer, UI loop  |
| `gl_loader.h` | Loads OpenGL 3.3 entry points at runtime             |
| `dependency/` | GLFW and Dear ImGui submodules                       |
