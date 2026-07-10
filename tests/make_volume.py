#!/usr/bin/env python3
"""Create a background volume mesh (MEDIT .msh) from a surface OFF file.

The volume is built from the surface *points only* (the triangles are ignored):
we take the input vertices, add the 8 corners of a bounding box enlarged by 10%
so the whole surface sits strictly inside, then tetrahedralize the combined
point set with scipy's Delaunay (Qhull). The result is a convex background mesh
that fills the enlarged box -- exactly the kind of volume mesh_generator embeds a
surface into (mesh_generator surface.off volume.msh).

The .msh is written in the minimal MEDIT dialect that src/VolumeRemesher/main.cpp
read_MEDIT_file() understands: header tokens up to "Vertices", the vertex block,
then immediately the "Tetrahedra" block (1-indexed), nothing in between.

Usage:
    make_volume.py input.off [-o output.msh] [--scale 1.10]
    make_volume.py a.off b.off c.off        # batch; each -> <name>.msh next to it
"""

import argparse
import sys

import numpy as np
from scipy.spatial import Delaunay


def read_off_points(path):
    """Return the (N, 3) vertex array of an OFF file. Faces are ignored.

    Handles the 'OFF' magic on its own line, comment lines (starting with '#',
    as emitted by meshio) and the 'nv nf ne' counts line, mirroring the tolerant
    parser in main.cpp read_OFF_file()."""
    tokens = []
    with open(path, "r") as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            tokens.extend(line.split())

    if not tokens:
        raise ValueError(f"{path}: empty file")
    # First token is the OFF magic (OFF / COFF / NOFF / STOFF / ...).
    if not tokens[0].endswith("OFF"):
        raise ValueError(f"{path}: not an OFF file (first token '{tokens[0]}')")

    idx = 1
    nv, nf, ne = (int(tokens[idx]), int(tokens[idx + 1]), int(tokens[idx + 2]))
    idx += 3

    coords = np.array(tokens[idx : idx + 3 * nv], dtype=np.float64)
    if coords.size != 3 * nv:
        raise ValueError(f"{path}: expected {nv} vertices, found {coords.size // 3}")
    return coords.reshape(nv, 3)


def bounding_box_corners(points, scale):
    """8 corners of the point cloud's axis-aligned box enlarged by (scale-1)*100%."""
    mn = points.min(axis=0)
    mx = points.max(axis=0)
    center = 0.5 * (mn + mx)
    half = 0.5 * (mx - mn) * scale
    # Guard against a flat model (zero extent on some axis) producing a
    # degenerate, coplanar corner set that Qhull cannot tetrahedralize.
    half = np.where(half <= 0.0, 1e-4, half)
    signs = np.array(
        [[sx, sy, sz] for sx in (-1.0, 1.0) for sy in (-1.0, 1.0) for sz in (-1.0, 1.0)],
        dtype=np.float64,
    )
    return center + signs * half


def make_volume(points, scale):
    """Return (vertices (M,3), tets (T,4) 0-indexed) tetrahedralizing points+box."""
    corners = bounding_box_corners(points, scale)
    verts = np.vstack([points, corners])
    tets = Delaunay(verts).simplices  # (T, 4), indices into verts
    return verts, tets


def write_medit_msh(path, verts, tets):
    """Write the minimal MEDIT .msh read_MEDIT_file() expects (1-indexed tets)."""
    with open(path, "w") as f:
        f.write("MeshVersionFormatted 1\n")
        f.write("Dimension 3\n")
        f.write("Vertices\n")
        f.write(f"{len(verts)}\n")
        for x, y, z in verts:
            # %.17g round-trips a double exactly, so mesh_generator reads back
            # the same coordinates scipy computed.
            f.write(f"{x:.17g} {y:.17g} {z:.17g} 0\n")
        f.write("Tetrahedra\n")
        f.write(f"{len(tets)}\n")
        for a, b, c, d in tets + 1:  # MEDIT is 1-indexed
            f.write(f"{a} {b} {c} {d} 0\n")
        f.write("End\n")


def process(in_path, out_path, scale):
    points = read_off_points(in_path)
    verts, tets = make_volume(points, scale)
    write_medit_msh(out_path, verts, tets)
    print(
        f"{in_path} -> {out_path}: "
        f"{len(points)} surface pts + 8 box corners = {len(verts)} verts, "
        f"{len(tets)} tets"
    )


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("inputs", nargs="+", help="input surface .off file(s)")
    ap.add_argument("-o", "--output", help="output .msh (only valid with a single input)")
    ap.add_argument("--scale", type=float, default=1.10, help="bounding-box scale factor (default 1.10 = 10%% larger)")
    args = ap.parse_args(argv)

    if args.output and len(args.inputs) != 1:
        ap.error("-o/--output can only be used with a single input file")

    for in_path in args.inputs:
        out_path = args.output if args.output else (in_path.rsplit(".", 1)[0] + ".msh")
        process(in_path, out_path, args.scale)
    return 0


if __name__ == "__main__":
    sys.exit(main())
