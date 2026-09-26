# Vapor

Interactive 2D and 3D smoke simulation built with C++20, OpenGL 4.3, and Dear ImGui.

## Build and run

Requires CMake 3.20+, a C++20 compiler, and an OpenGL 4.3-capable graphics driver.
GLFW and Dear ImGui are included as Git submodules.

```powershell
git submodule update --init --recursive
cmake -S . -B build
cmake --build build --config Release
./build/Release/vapor.exe
```

The executable path above is for Windows with Visual Studio. With a single-configuration
generator, set `-DCMAKE_BUILD_TYPE=Release` when configuring and run `build/vapor`.

Choose a solver in the startup picker and click **Start simulation**.
Close the simulation window to return to the picker. Settings are controlled in
the application; no configuration file is needed. CMake copies required meshes
beside the executable.

Windows builds request the high-performance GPU through GLFW's
Optimus/PowerXpress exports, preferring the NVIDIA GPU on hybrid NVIDIA laptops.
Restart the executable after rebuilding. The 3D box controls display the actual
compute device. An explicit Windows graphics preference can override this request.

## Solvers

| Solver | Domain | Execution |
| --- | --- | --- |
| MAC 2D Simple | Regular grid | CPU, semi-Lagrangian |
| MAC 2D Reflection | Regular grid | CPU, semi-Lagrangian or MacCormack smoke transport |
| MAC 3D Simple / Reflection | Regular grid | CPU or GPU, semi-Lagrangian or MacCormack |
| Simplicial2D-Teapot | Triangular teapot silhouette | CPU |
| Simplicial2D-Box | Square filled with triangles | CPU |
| Simplicial3D-Box | Cube filled with tetrahedra | GPU advection/recovery + CPU support |
| Simplicial3D-Bunny | Tetrahedral bunny volume | CPU |

Simplicial solvers transport circulation and recover incompressible flow.
They have no viscosity. The bunny uses 9,921 tetrahedra approximating the imported
bunny shape. The 3D box scene uses OpenGL compute for circulation
advection, smoke advection and the potential solve. Boundary reconstruction,
forces, mesh setup and volume resampling remain on the CPU; the controls show
the active GPU device. This is a hybrid solver, not a fully GPU-resident pipeline.
The experimental bunny GPU backend remains outside the picker because its
sustained CPU/GPU comparison currently fails.

## Controls

- **Pause / Reset:** stop the simulation or restart it.
- **Source, buoyancy, decay:** control smoke concentration, upward force, and fading.
- **Emitter stir:** perturb the teapot flow or apply a local rotating force in the bunny.
  Zero disables stirring. The bunny also has stir radius and frequency controls.
- **Emitter position and radius:** move and resize the source.
- **Simulation speed:** adjust time advanced per frame in the simplicial scenes.
- **Advanced:** solver iterations, tolerances, and other scene-specific settings.

The bunny opens with smoke and a thin exterior outline, leaving the interior visible.
Mesh inspection modes show primal edges in white and dual edges in green, with
independent toggles. Drag to orbit and scroll to zoom in 3D.
**Save render (.ppm)** exports the simplicial view without the controls.

Smoke is replenished continuously while the source is enabled. A steady-looking
plume does not mean emission has stopped. **Persistent smoke preset** disables decay,
but scalar transport is not mass-conservative and does not guarantee uniform filling.

## Tests

```powershell
ctest --test-dir build -C Release --output-on-failure
```

The simplicial 3D box uses a flux-derived slip-wall closure without a circulation
limiter or physical viscosity. Tangential wall velocities are projected from the
recovered interior field; boundary circulation is derived from that trace rather
than evolved as an independent constraint. This changes the boundary discretization:
it does not reproduce the paper's independently transported wall circulation.
Box timesteps use CFL subdivision including existing smoke and incoming buoyancy.
The plume, sustained unforced transport, and timestep-refinement regressions check
the resulting scheme; they do not establish unconditional stability.
The box scenes use a square mesh with 32 subdivisions per side (2,048 triangles)
and a side-length-3 tetrahedral cube generated at 10 subdivisions per axis
(1,331 vertices, 6,177 tetrahedra, `examples/simplicial3d/box-dense.tet`).
The finer 16- and 32-subdivision tetrahedral assets remain available for tests.
The 3D box emitter uses the same capped radial profile as 2D:
`min(1, 3 * sourceStrength * exp(-3.5 * (distance / radius)^2))` inside the
source sphere. The cap controls injected concentration; existing density is
preserved when the source shrinks or turns off.
The older 343-vertex cube is retained as
a numerical regression fixture. The 3D box opens with
smoke and the full domain outline. Mesh inspection and its optional X cutaway
remain available in the visualization controls. The 2D simplicial scene starts
with smoke decay disabled; linear scalar interpolation still introduces blur.
The tetrahedral asset is validated for full box volume, boundary faces, and
in-box circumcenters when loaded.

For numerical tests without the graphical application:

```powershell
cmake -S . -B build-headless -DVAPOR_BUILD_APP=OFF
cmake --build build-headless --config Release
ctest --test-dir build-headless -C Release --output-on-failure
```

## Code layout

- `src/app/`: startup picker and simulation controls.
- `src/solvers/`: MAC and simplicial solvers, including GPU experiments.
- `src/renderer/`: smoke rendering, outlines, and mesh inspection.
- `tests/`: numerical and graphics checks.
- `examples/`: mesh assets used by simulations and tests.
