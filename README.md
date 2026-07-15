# Volume remeshing

This code implements a variation of the volume meshing algorithm described in "**Convex Polyhedral Meshing for Robust Solid Modeling**" by Lorenzo Diazzi and Marco Attene (ACM Trans Graphics Vol 40, N. 6, Procs of SIGGRAPH Asia 2021). You may download a copy of the paper preprint here: http://arxiv.org/abs/2109.14434
The original code is available here: https://github.com/MarcoAttene/VolumeMesher

This fork adds, on top of the original algorithm:

* **Optional, combinable inputs.** A run takes any non-empty combination of an input **surface** (triangles), a set of **edges**, and a set of **points** — see [Inserting edges and points](#inserting-edges-and-points).
* **Exact surface tracking.** For each output element you can recover which input primitive it came from: output faces &rarr; input triangles, output edges &rarr; input edges, output vertices &rarr; input points — see [Tracking output back to the input](#tracking-output-back-to-the-input).
* **Embedding a surface (and edges/points) in an existing tetrahedral mesh**, programmatically via `embed.h`.

## Building

Clone this repository and configure with CMake:
```
git clone https://github.com/MarcoAttene/VolumeRemesher
cd VolumeRemesher
mkdir build
cd build
cmake ..
```
This produces a build configuration for your system (a `mesh_generator.sln` on Windows/MSVC, a `Makefile` on Linux/macOS). Build the `mesh_generator` target as usual. The integration tests are built by default (`VOLUMEREMESHER_BUILD_TESTS=ON`) — run them with `ctest`, or configure with `-DVOLUMEREMESHER_BUILD_TESTS=OFF` to skip them. (Some heavy tracking tests only run in a Release build.)

## Command-line usage

```
mesh_generator [options] [inputfile_A.off] [bool_opcode inputfile_B.off] [-e edges.off] [-p points.off]
```

The three inputs — the surface `inputfile_A.off`, the `-e` edges and the `-p` points — are **all optional**; provide any non-empty combination (only a surface, only edges, only points, or any mix). Run `mesh_generator` with no arguments to print the built-in help.

### Options

| Option | Meaning |
| --- | --- |
| `-v` | Verbose mode (prints progress and timings). |
| `-l` | Logging mode (appends a line of timings to `mesh_generator.log`). |
| `-s` | Save the mesh bounding surface ("skin") to `skin.off`. |
| `-b` | Save the subdivided input constraints ("black faces") to `black_faces.off`. |
| `-t` | Tetrahedralize the output: write a tet mesh to `volume.tet` (instead of the default polyhedral `volume.msh`). |
| `-r` | With `-t`, also write exact **rational** vertex coordinates to `volume.tet.rational`. |
| `-a` | Keep **all** cells: skip the inside/outside min-cut classification and tetrahedralize the whole domain. This makes the `-t` tracking **exact** (no cell is deleted, so every input primitive is realized in the output). Recommended whenever you use the tracking output. |
| `-f` | Fill holes: cap open boundaries with ear-clipped triangles before meshing (off by default; only affects the solid-output min-cut path — `-a` needs no capping). |
| `-e edges.off` | Also insert these edges into the output (implies `-a`). See [input formats](#input-formats). |
| `-p points.off` | Also insert these points into the output (implies `-a`). See [input formats](#input-formats). |
| `bool_opcode` | With two surfaces, combine them with a boolean op: `U` union (A∪B), `I` intersection (A∩B), `D` difference (A\B). |

### Input formats

* **Surface** (`inputfile_A.off`, and `inputfile_B.off` for booleans): a standard triangle **OFF** file.
* **Edges** (`-e`): an OFF-style file whose header is `nv ne 0`, followed by `nv` vertex lines `x y z`, then `ne` edge records `2 i j` (the `2` is the number of endpoints, mirroring how OFF faces start with `3`). See `models/cube_subdiv_edges.off`.
* **Points** (`-p`): an OFF-style file whose header is `np 0 0`, followed by `np` vertex lines `x y z`. See `models/cube_subdiv_points.off`.
* **Existing tetrahedral mesh** (for embedding a surface in it): a `.tet` or `.msh` (MEDIT) file passed as `inputfile_B` (ASCII). It is recognized by its extension, and a `bool_opcode` must still be given as a separator (it is ignored for embedding): `mesh_generator surface.off U mesh.tet`.

### Output files

| File | Written when | Contents |
| --- | --- | --- |
| `volume.msh` | default | Polyhedral volume mesh. |
| `volume.tet` | `-t` | Tetrahedral mesh: `<nv> vertices` / `<nt> tets`, then `nv` coordinate lines, then `nt` lines `4 a b c d` (0-based vertex indices). |
| `skin.off` | `-s` | Bounding surface. |
| `black_faces.off` | `-b` | Subdivided input constraints. |
| `volume.tet.rational` | `-t -r` | Exact rational coordinates, same order/indexing as the `volume.tet` vertex block. |
| `volume.tet.groups` | `-t` | Coplanar grouping of the input triangles. |
| `volume.tet.triangleprov` | `-t` | Triangle (face) provenance: per input coplanar group, the output faces tiling it. |
| `volume.tet.edgeprov` | `-t -e` | Edge provenance: per input edge, the output tet edges tiling it. |
| `volume.tet.pointprov` | `-t -p` | Point provenance: per input point, the equal output vertex. |

## Inserting edges and points

Edges and points are forced into the tetrahedralization so they appear as genuine tet **edges** and **vertices** of the output. Internally, each input edge `AB` is pinned by two small triangles sharing `AB`, and each input point `P` by three small triangles meeting at `P`; these helper triangles are excluded from the surface tracking. Insertion implies `-a` (keep all cells), because the pinned edges/points do not bound a solid and would otherwise be discarded by the inside/outside classification.

```
# a surface plus extra edges and points, with exact tracking
mesh_generator models/sphere.off -t -r -a -e models/sphere_edges.off -p models/sphere_points.off

# only edges (no surface); the domain is the bounding box of the edges
mesh_generator -t -r -a -e models/cube_subdiv_edges.off

# only points
mesh_generator -t -r -a -p models/cube_subdiv_points.off
```

## Tracking output back to the input

With `-t` the mesher writes sidecar files that map each input primitive to the output primitives it becomes. Use `-a` so the tracking is **exact** (every input primitive is realized). The three provenance sidecars share the **same shape** — keyed by the input primitive (coplanar triangle group / edge / point), each output primitive is given as an **output tet index plus that tet's vertex indices** (a face is `tet v0 v1 v2`, an edge is `tet v0 v1`, a point is `tet vertex`). All indices are 0-based into `volume.tet`'s tet and vertex blocks; `volume.tet.rational` gives the exact coordinates of each output vertex in the same order.

### `volume.tet.groups` — coplanar groups of the input triangles

The input triangles are partitioned into maximal **edge-connected, exactly-coplanar groups** (a flat region is one group; a triangle with no coplanar neighbour is a singleton group). Format:
```
# input_triangle_id coplanar_group_id
<num_input_triangles> <num_groups>
<input_triangle_id> <coplanar_group_id>
...
```
`input_triangle_id` is the triangle's index in the **input OFF file** (degenerate/collinear input triangles are dropped by the mesher, so this is not simply `0..n-1`). A group is a 1-to-1 map to a single input triangle iff it has exactly one triangle.

### `volume.tet.triangleprov` — triangle (face) provenance

One line per coplanar group, listing the output tet faces that tile it:
```
# group_id num_faces tet v0 v1 v2 ...
<num_groups>
<group_id> <num_faces> <tet v0 v1 v2> <tet v0 v1 v2> ...
```
Each face is given as an output tet index `tet` and the three output vertex indices `v0 v1 v2` of that tet's face (`v0 v1 v2` are three of tet `tet`'s four vertices). A face on two exactly-coplanar surfaces that overlap on the same plane is listed under both groups. Summing the (distinct) face areas of a group reconstructs the group's area.

### `volume.tet.edgeprov` — edge provenance

One line per input edge, listing the output tet edges that tile it (their union is exactly the input segment):
```
# edge_id num_out_edges tet v0 v1 ...
<num_input_edges>
<edge_id> <num_out_edges> <tet v0 v1> <tet v0 v1> ...
```
`edge_id` is the input edge's index in the `-e` file; each output edge is a tet index `tet` and the two output vertex indices `v0 v1` of that tet's edge.

### `volume.tet.pointprov` — point provenance

One line per input point:
```
# point_id tet out_vertex(-1 -1 if absent)
<num_input_points>
<point_id> <tet> <out_vertex>
```
`out_vertex` is the output vertex equal to the input point and `tet` an output tet containing it (both `-1` if the point did not survive into the output — should not happen with `-a`).

### Independent verifier

`tests/verify_tracking.cpp` re-derives all of the above in **exact rational arithmetic** and checks it (surface groups tile exactly, each edge is tiled by its output edges, each point equals its output vertex). It is a good reference for how to consume the sidecars:
```
mesh_generator models/cube_subdiv.off -t -r -a -e models/cube_subdiv_edges.off -p models/cube_subdiv_points.off
verify_tracking models/cube_subdiv.off volume.tet --strict \
    --edges models/cube_subdiv_edges.off --points models/cube_subdiv_points.off
# for an edges/points-only run (no surface), pass "none" as the surface:
verify_tracking none volume.tet --strict --edges models/cube_subdiv_edges.off
```

## Examples on the bundled meshes

```
# 1) Volume mesh (polyhedral) of a closed surface
mesh_generator models/bunny.off                       # -> volume.msh

# 2) Tetrahedral mesh of the enclosed solid
mesh_generator models/bunny.off -t                    # -> volume.tet

# 3) Boolean union of two overlapping surfaces
mesh_generator models/cube_subdiv.off U models/sphere.off   # -> volume.msh

# 4) Exact surface tracking of a tetrahedralization
mesh_generator models/cube_subdiv.off -t -r -a        # -> volume.tet(.prov/.groups/.rational)

# 5) Insert edges and points, with tracking
mesh_generator models/sphere.off -t -r -a \
    -e models/sphere_edges.off -p models/sphere_points.off

# 6) Embed a surface into an existing tetrahedral mesh (.tet/.msh as inputfile_B).
#    The bool_opcode (here U) is a required separator and is ignored for embedding;
#    the surface must lie inside the tet mesh's domain.
mesh_generator surface.off U input_tetmesh.tet        # -> volume.msh
```

## Programmatic API

To embed a set of triangles (and, optionally, edges and points) into an existing tetrahedral mesh from C++, see `embed.h` (`embed_tri_in_poly_mesh`). It returns the polyhedral/tetrahedral output plus the same provenance information described above (faces on input, and per-edge / per-point provenance). Note: in the embedding path the edges/points — and the small helper-triangle apexes they generate — must lie inside the input tet mesh's domain (a debug assertion flags violations).

## Notes

* The pipeline uses exact indirect predicates, so the output is deterministic and byte-identical across platforms/compilers (this is checked in the integration tests).

* We tested this code on macOS and Windows; it should also work on Linux.

|:warning: WARNING: Apparently, CLANG does not support a fully IEEE compliant floating point environment which is necessary to guarantee that indirect predicates work as expected. The only way we found to guarantee correctness on this compiler was to disable all optimizations. Please be aware of this fact should you notice a performance degradation in your experiments. |
| --- |
