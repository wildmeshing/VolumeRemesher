#!/usr/bin/env python3
"""Export a 2D arrangement (input segments + output triangulation) for ParaView.

Given the input OBJ line soup and the output of `mesh_generator <in.obj> -2d`, writes three
ASCII .vtu files that are meant to be opened together and overlaid:

    <prefix>_input.vtu         the input segments, one line cell each
                               cell data: input_segment_id  (0 .. n_segments-1)

    <prefix>_output_tris.vtu   the output triangulation
                               cell data: triangle_id

    <prefix>_output_edges.vtu  the output triangle edges that carry provenance, one line cell
                               per (input segment, output edge) pair
                               cell data: input_segment_id   -- the segment it came from
                                          n_covering         -- how many input segments cover
                                                                this edge (>1 means the input
                                                                had overlapping collinear
                                                                segments)
                                          piece_index        -- position along its segment,
                                                                0..k-1 from the segment's first
                                                                endpoint

The point of the shared `input_segment_id` is that colouring `_input` and `_output_edges` by the
same array makes the provenance visible directly: every input segment and the chain of output
edges that tiles it get the same colour. Set the colour map to a categorical / random lookup
table so neighbouring ids are distinguishable.

By default the input segments are lifted slightly in z so they render on top of the triangulation
instead of z-fighting with it. Pass --lift 0 to keep everything coplanar.

Usage:
    arrangement_to_vtu.py <input.obj> <triangulation.off> [-o PREFIX] [--lift FRACTION]

Example:
    mesh_generator models/2d/england.obj -2d
    tests/arrangement_to_vtu.py models/2d/england.obj triangulation.off -o england
    # then open england_input.vtu, england_output_tris.vtu, england_output_edges.vtu

No third-party dependencies: the .vtu files are written directly (unlike tet_to_vtu.py, which
uses meshio). ParaView reads ASCII VTU natively.
"""

import argparse
import os
import sys

VTK_LINE = 3
VTK_TRIANGLE = 5


def read_obj_segments(path):
    """Read an OBJ line soup. Returns (points, segments); z is dropped."""
    pts, segs = [], []
    with open(path) as f:
        for line in f:
            s = line.split()
            if not s:
                continue
            if s[0] == "v":
                pts.append((float(s[1]), float(s[2])))
            elif s[0] == "l":
                # "l" may hold a polyline; emit consecutive pairs. Indices are 1-based and may
                # be negative (relative to the end).
                idx = []
                for tok in s[1:]:
                    v = int(tok.split("/")[0])
                    idx.append(v - 1 if v > 0 else len(pts) + v)
                segs.extend(zip(idx[:-1], idx[1:]))
    return pts, segs


def read_off(path):
    """Read the OFF written by -2d. Returns (points, triangles); z is dropped."""
    with open(path) as f:
        tok = []
        for line in f:
            line = line.split("#")[0]
            tok.extend(line.split())
    if not tok or tok[0] != "OFF":
        sys.exit(f"{path}: not an OFF file")
    nv, nf = int(tok[1]), int(tok[2])
    p = 4
    pts = []
    for _ in range(nv):
        pts.append((float(tok[p]), float(tok[p + 1])))
        p += 3
    tris = []
    for _ in range(nf):
        n = int(tok[p])
        if n != 3:
            sys.exit(f"{path}: face with {n} vertices, expected triangles")
        tris.append((int(tok[p + 1]), int(tok[p + 2]), int(tok[p + 3])))
        p += 4
    return pts, tris


def read_segmentprov(path):
    """Read the .segmentprov sidecar -> list (per input segment) of (tri, v0, v1) triples."""
    with open(path) as f:
        tok = f.read().split()
    n = int(tok[0])
    p = 1
    out = []
    for _ in range(n):
        p += 1  # the segment id, which is just the line index
        m = int(tok[p])
        p += 1
        edges = []
        for _ in range(m):
            edges.append((int(tok[p]), int(tok[p + 1]), int(tok[p + 2])))
            p += 3
        out.append(edges)
    return out


def write_vtu(path, points, cells, cell_type, cell_data, lift=0.0):
    """Write an ASCII VTK UnstructuredGrid. `points` are 2D; z is set to `lift`."""
    npts, ncells = len(points), len(cells)
    with open(path, "w") as f:
        f.write('<?xml version="1.0"?>\n')
        f.write('<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">\n')
        f.write("  <UnstructuredGrid>\n")
        f.write(f'    <Piece NumberOfPoints="{npts}" NumberOfCells="{ncells}">\n')

        f.write("      <Points>\n")
        f.write('        <DataArray type="Float64" NumberOfComponents="3" format="ascii">\n')
        for x, y in points:
            f.write(f"          {x!r} {y!r} {lift!r}\n")
        f.write("        </DataArray>\n      </Points>\n")

        f.write("      <Cells>\n")
        f.write('        <DataArray type="Int64" Name="connectivity" format="ascii">\n')
        for c in cells:
            f.write("          " + " ".join(str(i) for i in c) + "\n")
        f.write("        </DataArray>\n")
        f.write('        <DataArray type="Int64" Name="offsets" format="ascii">\n')
        off = 0
        for c in cells:
            off += len(c)
            f.write(f"          {off}\n")
        f.write("        </DataArray>\n")
        f.write('        <DataArray type="UInt8" Name="types" format="ascii">\n')
        for _ in cells:
            f.write(f"          {cell_type}\n")
        f.write("        </DataArray>\n      </Cells>\n")

        if cell_data:
            first = next(iter(cell_data))
            f.write(f'      <CellData Scalars="{first}">\n')
            for name, values in cell_data.items():
                f.write(f'        <DataArray type="Int64" Name="{name}" format="ascii">\n')
                for v in values:
                    f.write(f"          {v}\n")
                f.write("        </DataArray>\n")
            f.write("      </CellData>\n")

        f.write("    </Piece>\n  </UnstructuredGrid>\n</VTKFile>\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("obj", help="the input OBJ line soup")
    ap.add_argument("off", help="triangulation.off written by -2d")
    ap.add_argument("-o", "--out-prefix", default=None,
                    help="output prefix (default: the OBJ's basename)")
    ap.add_argument("--lift", type=float, default=0.005,
                    help="lift the input segments above the triangulation by this fraction of "
                         "the bounding-box diagonal, so they do not z-fight (default 0.005; "
                         "use 0 to keep everything coplanar)")
    args = ap.parse_args()

    prefix = args.out_prefix or os.path.splitext(os.path.basename(args.obj))[0]

    in_pts, in_segs = read_obj_segments(args.obj)
    out_pts, out_tris = read_off(args.off)
    prov = read_segmentprov(args.off + ".segmentprov")

    if len(prov) != len(in_segs):
        sys.exit(f"provenance has {len(prov)} entries but the OBJ has {len(in_segs)} segments")

    xs = [p[0] for p in out_pts]
    ys = [p[1] for p in out_pts]
    diag = ((max(xs) - min(xs)) ** 2 + (max(ys) - min(ys)) ** 2) ** 0.5 if out_pts else 1.0
    lift = args.lift * diag

    # --- input ---
    write_vtu(f"{prefix}_input.vtu", in_pts, in_segs, VTK_LINE,
              {"input_segment_id": list(range(len(in_segs)))}, lift=lift)

    # --- output triangulation ---
    write_vtu(f"{prefix}_output_tris.vtu", out_pts, out_tris, VTK_TRIANGLE,
              {"triangle_id": list(range(len(out_tris)))})

    # --- output edges, tagged with the segment they came from ---
    # An edge covered by several input segments (overlapping collinear input) is emitted once per
    # covering segment; n_covering says how many.
    cover = {}
    for edges in prov:
        for _, v0, v1 in edges:
            cover[(min(v0, v1), max(v0, v1))] = cover.get((min(v0, v1), max(v0, v1)), 0) + 1

    cells, seg_ids, n_cov, piece = [], [], [], []
    for sid, edges in enumerate(prov):
        for k, (_, v0, v1) in enumerate(edges):
            cells.append((v0, v1))
            seg_ids.append(sid)
            n_cov.append(cover[(min(v0, v1), max(v0, v1))])
            piece.append(k)

    write_vtu(f"{prefix}_output_edges.vtu", out_pts, cells, VTK_LINE,
              {"input_segment_id": seg_ids, "n_covering": n_cov, "piece_index": piece})

    n_split = sum(1 for e in prov if len(e) > 1)
    n_over = sum(1 for v in cover.values() if v > 1)
    print(f"{prefix}_input.vtu         {len(in_pts)} points, {len(in_segs)} segments")
    print(f"{prefix}_output_tris.vtu   {len(out_pts)} points, {len(out_tris)} triangles")
    print(f"{prefix}_output_edges.vtu  {len(cells)} edges "
          f"({n_split} input segments were split, {n_over} output edges covered by >1 segment)")
    print(f"\nOpen all three in ParaView and colour _input and _output_edges by "
          f"'input_segment_id' with a categorical lookup table.")


if __name__ == "__main__":
    main()
