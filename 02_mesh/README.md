**Standing: an illustration. No equation of the flow is solved on this grid, and no load depends on it.**

# 02_mesh

The load model (package `unistall`) is a set of ordinary differential equations marched in time. It has no
spatial mesh and does not read this grid.

`generate_mesh.py` writes a structured O-grid round the closed NACA 0012 outline of `01_geometry` (the same
outline as the panels of the reconstruction and the drawings), with its quality measures and four pictures.
`field_on_mesh.py` evaluates the potential-flow reconstruction of `unistall/flowfield.py` at the grid nodes,
for each case at the instant of peak lift, for one picture: the wall line takes the values on the surface
itself, and nodes nearer the wall than the distance within which the field shows the panel junctions carry
no value. That is the one use the grid has.

Run from the repository root: `PYTHONPATH=. python3 02_mesh/generate_mesh.py`, and after the solution stage
`PYTHONPATH=. python3 02_mesh/field_on_mesh.py`.
