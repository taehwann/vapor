# Stanford bunny tetrahedral test

443 vertices, 1,996 tetrahedra. Original Stanford bunny scan: Stanford
University Computer Graphics Laboratory (Greg Turk and Marc Levoy).

Volume mesh: Barbara Cutler et al., *Simplification and Improvement of
Tetrahedral Models for Simulation* (SGP 2004),
https://people.csail.mit.edu/bmcutler/PROJECTS/SGP04/meshes/index.html
(`bunny_2K.tet`). Both source material regions are retained.
Stanford source and usage terms: https://graphics.stanford.edu/data/3Dscanrep/

Reproduce with `python tools/import_bunny_mesh.py PATH/TO/bunny_2K.tet`.
The adjacent JSON records the source hash and coordinate transform.

## Fluid mesh

`bunny.tet` is the preserved imported source mesh. It is not used directly by
the fluid solver: its signed circumcentric metrics fail the positivity checks.

`bunny-fluid.tet` is the simulation mesh: 2,317 vertices and 9,921 tetrahedra.
Reproduce with `python tools/remesh_bunny.py --resolution 18` using NumPy/SciPy
from `tools/simplicial-mesher-requirements.txt`. It selects whole, well-centered
BCC lattice tetrahedra by containment in the source volume, removes nonmanifold
surface edge/vertex contacts, and retains the largest face-connected component.
It approximates the bunny surface with lattice triangles; it does not preserve
every original boundary triangle.

Select **Simplicial3D-Bunny** in Vapor for the smoke simulation. White wires show
unique primal edges. Green wires use the actual solver dual: shared tetrahedron
circumcenters, boundary triangle circumcenters, and boundary edge midpoints.
Surface dual edges are subdivided for tangential interpolation; this does not
change their geometric paths. Every dual face closes; the dual graph is connected.

Toggle both layers independently. Disable Boundary surface wire for a pure
dual-only render. Save render writes `simplicial3d-render.ppm` in the working
directory. See [validation and method differences](../../docs/simplicial-bunny-validation.md).
