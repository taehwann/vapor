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

## Solvers

| Solver | Domain | Execution |
| --- | --- | --- |
| MAC 2D Simple | Regular grid | CPU, semi-Lagrangian |
| MAC 2D Reflection | Regular grid | CPU, semi-Lagrangian or MacCormack smoke transport |
| MAC 3D Simple / Reflection | Regular grid | CPU or GPU, semi-Lagrangian or MacCormack |
| Simplicial2D-Teapot | Triangular teapot silhouette | CPU |
| Simplicial3D-Bunny | Tetrahedral bunny volume | CPU |

Simplicial solvers transport circulation and recover incompressible flow.
They have no viscosity. The bunny uses 9,921 tetrahedra approximating the imported
bunny shape. The experimental simplicial GPU backend is not enabled in the picker.

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

The older simplicial box stability and transport tests have known failures.
The square and box remain test fixtures, not picker scenes.

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
