#!/usr/bin/env python3
"""Convert a mesh_generator tetrahedral mesh (.tet) into a ParaView file (.vtu).

mesh_generator writes this .tet format when run with -t (BSPcomplex::saveMesh in
tetrahedralize mode):

    <N> vertices
    <M> tets
    x y z                 (N lines: vertex coordinates, doubles)
    4 a b c d              (M lines: one tetrahedron; leading 4 = vertices per
                            cell, then four 0-based vertex indices)

The result is a VTK UnstructuredGrid (.vtu) with tetra cells that opens directly
in ParaView. VTU is binary+compressed, so it is typically smaller than the .tet.

Usage:
    tet_to_vtu.py volume.tet [out.vtu]      # default out = <input>.vtu
    tet_to_vtu.py a.tet b.tet ...           # batch; each -> <name>.vtu

Requires meshio (pip install meshio); uses pandas for a fast parse if available,
otherwise falls back to numpy.
"""

import argparse
import sys

import numpy as np
import meshio


def read_tet(path):
    """Return (points (N,3) float64, tets (M,4) int) from a .tet file."""
    with open(path) as f:
        nverts = int(f.readline().split()[0])
        ntets = int(f.readline().split()[0])

    # The vertex block (3 columns) and tet block (5 columns) have different
    # widths, so they must be read separately. pandas' C parser is much faster
    # than numpy on the tens-of-millions of rows these meshes can reach.
    try:
        import pandas as pd

        pts = pd.read_csv(
            path, sep=r"\s+", header=None, skiprows=2, nrows=nverts, dtype=np.float64
        ).to_numpy()
        tets = pd.read_csv(
            path, sep=r"\s+", header=None, skiprows=2 + nverts, nrows=ntets, dtype=np.int64
        ).to_numpy()[:, 1:5]
    except ImportError:
        from itertools import islice

        with open(path) as f:
            f.readline()
            f.readline()
            pts = np.loadtxt(islice(f, nverts), dtype=np.float64)
            tets = np.loadtxt(islice(f, ntets), dtype=np.int64)[:, 1:5]

    if pts.shape != (nverts, 3):
        raise ValueError(f"{path}: expected {nverts} vertices, parsed {pts.shape}")
    if tets.shape != (ntets, 4):
        raise ValueError(f"{path}: expected {ntets} tets, parsed {tets.shape}")
    return pts, tets


def convert(in_path, out_path):
    pts, tets = read_tet(in_path)
    meshio.write_points_cells(out_path, pts, [("tetra", tets)])
    print(f"{in_path}: {len(pts)} vertices, {len(tets)} tets -> {out_path}")


def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("inputs", nargs="+", help="input .tet file(s)")
    ap.add_argument("-o", "--output", help="output .vtu (only with a single input)")
    args = ap.parse_args(argv)

    if args.output and len(args.inputs) != 1:
        ap.error("-o/--output can only be used with a single input file")

    for in_path in args.inputs:
        out_path = args.output if args.output else (in_path.rsplit(".", 1)[0] + ".vtu")
        convert(in_path, out_path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
