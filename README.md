# Vapor

An interactive **3D smoke-plume** simulation written in C++20, rendered with OpenGL
4.3 volume ray-marching and Dear ImGui. The solver integrates the incompressible
Navier–Stokes equations on a **closed-box** staggered (MAC) grid using semi-Lagrangian
advection, an optional midpoint-reflect (Richardson extrapolation) velocity update,
and a red-black SOR pressure solve. Both the SOR solve and the semi-Lagrangian
advection can be toggled between CPU (OpenMP-parallelised) and GPU (OpenGL 4.3
compute shader) independently at runtime.

---

## Grid layout

The domain is `[0, boxSize]³` (adjustable 0.5–8.0, default 4.5) divided into
`n×n×n` cells (`n = 32` at runtime, cell size `h = boxSize / n`). Quantities live
on a **staggered (MAC) grid**:

| Field         | Location       | Size               | Index function                    |
| ------------- | -------------- | ------------------ | --------------------------------- |
| `smoke`       | cell centers   | `n³`               | `x + n*(y + n*z)`                 |
| `pressure`    | cell centers   | `n³`               | same as above                     |
| `divergence`  | cell centers   | `n³`               | same as above                     |
| `vx`          | x-faces        | `(n+1) × n × n`    | `x + (n+1)*(y + n*z)`            |
| `vy`          | y-faces        | `n × (n+1) × n`    | `x + n*(y + (n+1)*z)`            |
| `vz`          | z-faces        | `n × n × (n+1)`    | `x + n*(y + n*z)`                 |

`x` ranges `[0, n-1]` for cell/scalar fields, `[0, n]` for face-velocity fields
along their staggered axis. The extra `+1` in one dimension per velocity component
is baked into the index stride: `idX` uses `(n+1)` as the x-stride, `idY` uses
`(n+1)` as the y-stride, and `idZ` shares `idC`'s formula because z is the
outermost dimension — `n*n` advance per z-layer naturally covers `n+1` layers.

---

## Boundaries — fully closed solid walls

All 6 faces are **no-slip solid walls**: normal velocity is zeroed at every outer
face. The method `enforceVelocityBoundaries()` unconditionally sets:

| Face   | Assignment                     |
| ------ | ------------------------------ |
| -X, +X | `vx[idX(0,*,*)] = 0`, `vx[idX(n,*,*)] = 0` |
| -Y, +Y | `vy[idY(*,0,*)] = 0`, `vy[idY(*,n,*)] = 0` |
| -Z, +Z | `vz[idZ(*,*,0)] = 0`, `vz[idZ(*,*,n)] = 0` |

It is called:
- At the **start** of `project()` (before divergence computation)
- At the **end** of `project()` (after pressure correction)
- At the **start** of `step()` / `stepNoReflect()` (after emit + buoyancy)
- During the pressure correction loop itself, velocity updates are restricted to
  interior faces only: `x ∈ [1, n-1]`, `y ∈ [1, n-1]`, `z ∈ [1, n-1]`, which
  preserves the zero-wall condition without needing a separate boundary pass.

There are no open-boundary / outflow / clamping options in the current build.

---

## Sampling and interpolation

All field sampling uses **clamped** coordinate lookups before trilinear
interpolation. This ensures that traces which leave the domain return the nearest
boundary value rather than zero, which is correct for a closed box.

| Function           | Valid range per axis          |
| ------------------ | ----------------------------- |
| `sampleCell`       | `g ∈ [0, n-1]` for all three |
| `sampleVxField`    | `gx ∈ [0, n]`, `gy,gz ∈ [0, n-1]` |
| `sampleVyField`    | `gy ∈ [0, n]`, `gx,gz ∈ [0, n-1]` |
| `sampleVzField`    | `gz ∈ [0, n]`, `gx,gy ∈ [0, n-1]` |

For cell-centered fields, `x0 = int(g)`, `x1 = min(x0+1, n-1)`. For face-velocity
fields on the staggered axis, `x0 = min(int(g), n-1)`, `x1 = x0+1`.

---

## Emitter

A **spherical source** near the bottom of the domain:

| Property            | Default               | Description                          |
| ------------------- | --------------------- | ------------------------------------ |
| `emitterCenterX/Y`  | `0.5`                 | horizontal centre (normalised [0,1]) |
| `emitterCenterZ`    | `0.1`                 | low elevation (normalised [0,1])     |
| `emitterRadius`     | `0.035`               | sphere radius (normalised [0,1])     |
| `emitterKickBase`   | `0.8`                 | base upward velocity                 |
| `emitterKickScale`  | `0.8`                 | velocity multiplier × sourceStrength |

All positions scale proportionally with `boxSize`. Inside the sphere:
- `smoke[i] = 1.0f`
- Both `vz` faces of each emitter cell get `max(vz, kick)`
- Random stirring jitter is added to `vx` and `vy` faces, scaled by
  `stirStrength * dt * 60.0` and `(frand() - 0.5)`

---

## Buoyancy

For each cell `(x,y,z)`:

```
f = buoyancyVal * smoke[i] * dt
```

The force is split evenly across the two `vz` faces of the cell:

```
vz[idZ(x,y,z)]   += buoyancySplit * f    (z > 0, skip bottom wall)
vz[idZ(x,y,z+1)] += buoyancySplit * f    (z < n-1, skip top wall)
```

`buoyancySplit = 0.5` by default.

---

## Smoke dissipation

A per-frame density decay, applied AFTER advection:

```
sFactor = clamp(1.0 - smokeDecay * dt, 0.0, 1.0)
smoke[i] *= sFactor
```

`smokeDecay = 0.06` by default (tuned so the plume fades before reaching the top
wall in a closed box). Adjustable 0–2 via UI slider.

---

## Advection — semi-Lagrangian with optional MacCormack

All advection uses a **semi-Lagrangian backtrace**:

```
fwdScratch[i] = sampleField(field, position - velocity(position) * dt)
```

If MacCormack is disabled (`macCormackSmoke = false` or `macCormackVel = false`),
`fwdScratch` is copied directly to the output field.

### MacCormack correction

When enabled, a second pass computes a forward trace from `fwdScratch`:

```
backScratch[i] = sampleField(fwdScratch, position + velocity(position) * dt)
```

The corrected value is:

```
corrected = fwdScratch[i] + 0.5 * (field[i] - backScratch[i])
```

A **3×3×3 neighborhood limiter** prevents overshoot: the forward trace position
`posFwd = p - v*dt` is mapped to a grid cell `(gx,gy,gz)`, and the 8 surrounding
cells of `field` give `[fMin, fMax]`. The corrected value is clamped to this range.

### Boundary fallback

If either the forward or backward trace position leaves the domain
`[0, boxSize]³`, the cell uses `fwdScratch[i]` (the pure semi-Lagrangian result)
without MacCormack correction.

### CFL guard

For velocity advection, a CFL-based guard (`cfl > cflMc`, default `cflMc = 3.0`)
also disables MacCormack correction for cells with very high velocity.

---

## Pressure projection — Neumann boundary conditions

The pressure Poisson equation `Δp = ∇·u / dt` is solved with a **Red-Black
Gauss-Seidel SOR** method. Two backends are available, toggleable via the
"GPU SOR" checkbox:

| Backend | Description |
|---------|-------------|
| **CPU** (default) | OpenMP-parallelised red-black SOR, runs on host. Same algorithm as GPU. |
| **GPU** | OpenGL 4.3 compute shaders with SSBOs. Divergence, SOR, and velocity correction run entirely on the GPU via `glDispatchCompute`. Data is uploaded before and downloaded after each `project()` call. |

### GPU advection

Semi-Lagrangian advection with optional MacCormack correction can also be offloaded
to the GPU, toggleable via the "GPU Advection" checkbox:

| Backend | Description |
|---------|-------------|
| **CPU** (default) | OpenMP-parallelised advection on host. All fields (`smoke`, `vx`, `vy`, `vz`) are advected in-place with scratch vectors. |
| **GPU** | OpenGL 4.3 compute shaders with SSBOs. Two shader programs (`progAdvectSL` and `progAdvectMC`) handle all four field types (scalar smoke, x/y/z velocity faces) via a mode uniform. Each advection call uploads source and tracer velocity fields, dispatches one or two compute passes (semi-Lagrangian backtrace + optional MacCormack correction), and downloads the result. |

The GPU advection reuses the same velocity SSBOs created for the GPU SOR solver
and adds four additional buffers: `ssboSmoke`, `ssboFwd` (scratch for the
semi-Lagrangian backtrace), `ssboAdvectSrc`, and `ssboAdvectOut`. A simple
`progCopy` shader handles buffer-to-buffer copies when MacCormack is disabled.

The GLSL implementations mirror the CPU advection logic exactly:
- **Trilinear interpolation** for velocity sampling and field lookups, with the
  same clamped coordinate ranges and staggered-grid index functions (`idC`, `idX`,
  `idY`, `idZ`).
- **MacCormack correction** with boundary fallback (out-of-domain traces revert to
  pure semi-Lagrangian), CFL guard, and 2×2×2 neighbourhood clamping.
- Workgroup size 8×8×4 per dispatch, matching the divergence shader layout.

GPU advection and GPU SOR can be enabled independently. When both are active, the
entire advection–projection pipeline runs on the GPU with upload/download per
operation (future work: persistent GPU-resident fields to eliminate the round-trip).

### Divergence computation

```
div[i] = (vx[x+1] - vx[x] + vy[y+1] - vy[y] + vz[z+1] - vz[z]) / h
```

### SOR stencil with Neumann BCs

For each cell, all 6 neighbours are considered. If a neighbour lies outside the
domain, the homogeneous Neumann condition `∂p/∂n = 0` is enforced by using the
current cell's pressure as the ghost value:

```
if (!inside(a,b,c)) {
    sum += pressure[i];    // p_ghost = p_current (Neumann)
    ++terms;
    continue;
}
++terms;
sum += pressure[idC(a,b,c)];
```

The SOR update is:

```
p[i] = (1 - ω) * p[i] + ω * (sum - div[i] * h² / dt) / max(terms, 1)
```

where `ω = sorOmega = 1.95` by default.

### Velocity correction

Pressure gradient is subtracted from interior face velocities only (boundary
faces are skipped to maintain `v=0` at solid walls):

```
vx[x] -= dt * (p[x,y,z] - p[x-1,y,z]) / h    for x ∈ [1, n-1]
vy[y] -= dt * (p[x,y,z] - p[x,y-1,z]) / h    for y ∈ [1, n-1]
vz[z] -= dt * (p[x,y,z] - p[x,y,z-1]) / h    for z ∈ [1, n-1]
```

### Iteration count

`projectIterations = n * 2` by default (64 for n=32). Adjustable 10–500 via UI.

---

## Time integration — midpoint-reflect

Two modes are available, toggleable via the "Reflection" checkbox:

### Reflection mode (`step`, second-order)

1. Emit smoke, apply buoyancy, enforce wall boundaries
2. Save current velocity: `v0 = v`
3. Advect each velocity component by `halfDt` using `v0` as both source and tracer → `vTilde`
4. Project with `halfDt / halfIters`
5. Save clean projected velocity: `vDiv0 = v`
6. Compute reflected velocity: `vHat = 2*v - vTilde`
7. Advect `vHat` by `halfDt` using `vDiv0` as tracer
8. Project with `halfDt / halfIters` again
9. Advect smoke by full `dt`
10. Apply smoke dissipation

The reflection step `vHat = 2*v - vTilde` is a Richardson extrapolation that
increases the order of accuracy from first to second order in time. The key detail
is that the tracer `vDiv0` is captured AFTER the first projection — it is the
divergence-free velocity field without boundary contamination.

### Single-pass mode (`stepNoReflect`, first-order)

1. Emit smoke, apply buoyancy, enforce wall boundaries
2. Advect all velocity components by full `dt`
3. Project with full iterations
4. Advect smoke by full `dt`
5. Apply smoke dissipation

The iteration count for reflection is `std::max(1, projectIterations / 2)` per
half-step (roughly equal total work to the single-pass case).

---

## Velocity `idZ` indexing detail

The `vz` array is `n × n × (n+1)` elements (z-faces with z ∈ [0, n]). The index
formula `x + n*(y + n*z)` is **correct**:

- x ∈ [0, n-1], stride = 1
- y ∈ [0, n-1], stride = n
- z ∈ [0, n], stride = n*n

Maximum index = `(n-1) + n*(n-1) + n*n*n` = `n³ + n² - 1`, which equals the
total allocation `n*n*(n+1) - 1`. The formula is identical to `idC` because the
extra element is in the **outermost** (z) dimension, whose stride is `n*n`
regardless of whether z stops at `n-1` or `n`.

---

## Box resizing

When `boxSize` changes via the UI slider, all simulation fields (smoke, velocity,
pressure, divergence) are **cleared to zero** — a fresh start at the new domain
size. The emitter stays at its relative [0,1] position and all other parameters
are preserved.

---

## Volume rendering

A 3D texture (`GL_R16F`, half-float) is uploaded to the GPU every frame with
`glTexSubImage3D`. The texture wrapping is `GL_CLAMP_TO_EDGE` on all three axes.

### Ray-marching shader

The fragment shader:
1. Casts a ray from the camera through each pixel
2. Intersects the domain box `[0, boxSize]³` (uses front-face culling with
   `glCullFace(GL_FRONT)` so interior entry faces are rasterised)
3. Steps through the volume with step size `u_StepScale / u_GridRes`
4. At each sample where `density > 0.001`, accumulates alpha-composited smoke
   with an optional self-shadowing term
5. Exits early when opacity reaches 0.99

**Shadow calculation**: for each lit sample, a secondary shadow ray is cast
toward the light direction. Opacity is accumulated along the shadow ray with an
exponential falloff: `shadow = exp(-accumulated_opacity)`.

**Non-smoke pixels** are discarded (`discard`), revealing the clear color
background `(0.03, 0.04, 0.07)`.

### Camera

Orbit-style camera around the box centre. Mouse drag rotates (azimuth /
elevation), scroll wheel zooms. Controls work when the mouse is not captured by
Dear ImGui.

---

## ASCII framebuffer dump (debug)

Every 5th frame (when "Debug print" is checked), the RGB framebuffer is read back
via `glReadPixels`, downscaled to 80×24, and printed as ASCII art to the log file
using the ramp ` .:-=+*#%@` (dark to bright).

---

## Build

Requires CMake ≥ 3.20, a C++20 compiler, and OpenGL 4.3. GLFW and Dear ImGui are
included as git submodules under `dependency/`.

```powershell
git submodule update --init --recursive
cmake -S . -B build
cmake --build build --config Release
.\build\Release\vapor.exe
```

Visual Studio 2022 builds with `MSBuild` are also supported (tested with
`/p:Configuration=Release`).

---

## UI controls

| Control          | Default    | Range      | Description                                      |
| ---------------- | ---------- | ---------- | ------------------------------------------------ |
| Pause / Resume   | —          | —          | Freeze or advance the simulation                 |
| Reset            | —          | —          | Zero all simulation fields, keep all parameters  |
| Reflection       | true       | bool       | Toggle second-order midpoint-reflect vs single-pass |
| MC smoke         | true       | bool       | Toggle MacCormack correction for smoke advection |
| MC vel           | true       | bool       | Toggle MacCormack correction for velocity        |
| Debug print      | true       | bool       | Log diagnostic stats and ASCII framebuffer       |
| GPU SOR          | false      | bool       | Run SOR pressure solve on GPU via compute shaders|
| GPU Advection    | false      | bool       | Run semi-Lagrangian advection on GPU via compute shaders |
| Buoyancy         | 2.5        | 0–10       | Upward force strength                            |
| Source           | 1.0        | 0–3        | Emitter injection multiplier                     |
| Smoke decay      | 0.06       | 0–2        | Per-frame density dissipation rate               |
| Box size         | 4.5        | 0.5–8      | Domain extent `[0, boxSize]³`                    |
| Emitt radius     | 0.035      | 0.01–0.2   | Emitter sphere radius (normalised [0,1])         |
| Stir             | 0.5        | 0–2        | Random turbulence jitter in emitter              |
| Alpha mul        | 30.0       | 5–100      | Volume rendering opacity multiplier              |
| SOR its          | n*2 (64)   | 10–500     | SOR pressure solve iterations per projection     |
| MC CFL           | 3.0        | 0.5–10     | MacCormack CFL guard threshold                   |
| Step             | 0.5        | 0.1–3      | Ray-marching step scale                          |
| Shadow str       | 0.3        | 0–5        | Shadow self-occlusion strength                   |
| Shadow step      | 0.1        | 0.05–1     | Shadow ray step size                             |
| Light dir        | (0.3,0.4,1)| [-1,1] each | Directional light (auto-normalised)            |

---

## Logging

When running, `vapor_console.log` is written with per-frame statistics:

```
f    1 | ke=0.012345 maxDiv=0.000123 maxV=1.234 smoke=1.000 | REFLECT
```

Every 5 frames, an 80×24 ASCII frame dump is appended.

---

## Files

| File            | Purpose                                             |
| --------------- | --------------------------------------------------- |
| `main.cpp`      | `SmokeSim3D` (solver), `GpuSORSolver` (GPU pressure + advection), `VolumeRenderer`, UI loop |
| `gl_loader.h`   | Loads OpenGL 4.3 entry points at runtime (compute shader support) |
| `dependency/`   | GLFW and Dear ImGui submodules                      |
| `build/`        | CMake build output                                  |
| `vapor_console.log` | Runtime diagnostic log (git-ignored)           |
