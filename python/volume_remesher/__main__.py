"""``python -m volume_remesher``: the ``mesh_generator`` command line, argument for argument.

Same flags, same validation, same output file names and the same bytes in them. The pipeline
is the one the package wraps, so this is the direct check that the two cannot drift apart:
``python -m volume_remesher bunny.off -t`` and ``mesh_generator bunny.off -t`` write identical
``volume.tet`` files.
"""

from __future__ import annotations

import sys
from typing import Optional

from . import _core

USAGE = """
Usage: mesh_generator [-v | -l | -s | -b | -t] [inputfile_A.off] [bool_opcode inputfile_B.off] [-e edges.off] [-p points.off]

Defines the volume enclosed by the input OFF file(s) and saves a volume mesh to 'volume.msh'. The surface (inputfile_A.off), -e edges, and -p points are all optional; provide any non-empty combination.

Command line arguments:
-v = verbose mode
-l = logging mode (appends a line to mesh_generator.log)
-s = save the mesh bounding surface to 'skin.off'
-b = save the subdivided constraints to 'black_faces.off'
-t = triangulate/tetrahedrize output
-r = with -t, also write exact-rational vertex coords to 'volume.tet.rational'
-2d = 2D mode: read an OBJ line soup (v/l records) and write the arrangement
     of the input segments as a triangulation to 'triangulation.off'
-a = keep all cells: skip the in/out min-cut and tetrahedralize the whole
     domain, so the -t surface tracking is exact (no cells are deleted)
-f = fill holes: cap open boundaries with ear-clipped triangles (off by
     default; only affects the solid-output min-cut path)
-e edges.off  = also insert these edges into the output (implies -a); OFF
     with 'nv ne 0' then vertices then edge records '2 i j'
-p points.off = also insert these points into the output (implies -a); OFF
     with 'np 0 0' then vertex lines
bool_opcode: {U, I, D}
  U -> union (AuB),
  I -> intersection (A^B),
  D -> difference (A\\B)

The parallel passes run on VOLUMEREMESHER_NUM_THREADS threads if set, else on
the CPUs the process may run on.

Example:
mesh_generator ant.off U pig.off
"""


class _Exit(Exception):
    def __init__(self, message: str):
        super().__init__(message)
        self.message = message


def _error(message: str) -> None:
    # ip_error: the message on stderr and a failure exit status.
    raise _Exit(message)


def _run(argv: list[str]) -> int:
    triangulate = verbose = logging = surfmesh = blackfaces = False
    export_rational = fill_holes = two_d = keep_all_cells = False
    edge_file: Optional[str] = None
    point_file: Optional[str] = None
    file_a: Optional[str] = None
    file_b: Optional[str] = None
    bool_opcode = "0"

    if len(argv) < 2:
        sys.stdout.write(USAGE)
        return 0

    i = 1
    while i < len(argv):
        arg = argv[i]
        if arg.startswith("-"):
            if arg == "-2d":
                two_d = True
            elif arg[1:2] == "v":
                verbose = True
            elif arg[1:2] == "t":
                triangulate = True
            elif arg[1:2] == "l":
                logging = True
            elif arg[1:2] == "b":
                blackfaces = True
            elif arg[1:2] == "s":
                surfmesh = True
            elif arg[1:2] == "r":
                export_rational = True
            elif arg[1:2] == "f":
                fill_holes = True
            elif arg[1:2] == "a":
                keep_all_cells = True
            elif arg[1:2] == "e":
                if i + 1 >= len(argv):
                    _error("-e requires an edge file argument\n")
                i += 1
                edge_file = argv[i]
            elif arg[1:2] == "p":
                if i + 1 >= len(argv):
                    _error("-p requires a point file argument\n")
                i += 1
                point_file = argv[i]
            else:
                _error("Unknown option\n")
        elif file_a is None:
            file_a = arg
        elif bool_opcode == "0":
            bool_opcode = arg[0]
        elif file_b is None:
            file_b = arg
        else:
            _error("Too many args passed\n")
        i += 1

    two_input = bool_opcode != "0"

    if two_d:
        if file_a is None:
            _error("-2d needs an input .obj line soup\n")
        points, segments = _core.read_obj_segments(file_a)
        if verbose:
            print(f"Loaded {len(points)} points and {len(segments)} segments from {file_a}")
        arr = _core.arrange_2d(points, segments, export_rational, verbose)
        arr.save("triangulation.off")
        if verbose:
            print(f"Wrote triangulation.off: {len(arr.coords)} vertices, {len(arr.triangles)} triangles")
        return 0

    if verbose:
        if file_b is None:
            print("\nResolve auto-intersections and/or repair.\n")
            print(f"Loading {file_a}")
        else:
            print("\nBoolean operator: ", end="")
            if bool_opcode == "U":
                print(" union.\n")
            elif bool_opcode == "I":
                print(" intersection.\n")
            elif bool_opcode == "D":
                print(" difference.\n")
            else:
                print("INVALID\n")
                return 0
            print(f"Loading {file_a} and {file_b}.\n")

    if file_a is None and edge_file is None and point_file is None:
        _error("No input: provide a surface .off and/or -e edges.off and/or -p points.off\n")
    if file_a is None and two_input:
        _error("A boolean operation needs a surface input\n")

    file_b_is_tet = file_b is not None and file_b.endswith(".tet")
    file_b_is_msh = file_b is not None and file_b.endswith(".msh")
    embedsurf = file_b_is_tet or file_b_is_msh
    if embedsurf:
        two_input = False
    if two_input and file_b is None:
        _error("A boolean operation needs a second surface\n")
    if bool_opcode not in ("0", "U", "I", "D"):
        _error("Unknown boolean operation\n")

    output = (
        _core.OutputKind.BLACK_FACES
        if blackfaces
        else (_core.OutputKind.SKIN if surfmesh else _core.OutputKind.VOLUME)
    )

    A = _core.read_off(file_a, verbose) if file_a else None
    B = _core.read_off(file_b, verbose) if two_input else None
    edges = _core.read_edge_off(edge_file, verbose) if edge_file else None
    points = _core.read_points_off(point_file, verbose) if point_file else None

    if embedsurf:
        T = _core.read_tet(file_b, verbose) if file_b_is_tet else _core.read_medit(file_b, verbose)
        mesh = _core.embed_mesh(
            A, T, bool_opcode, edges, points, verbose, triangulate, export_rational, output
        )
    else:
        mesh = _core.make_mesh(
            A, file_a, bool_opcode, B, file_b, edges, points,
            verbose, logging, triangulate, export_rational, fill_holes, keep_all_cells, output,
        )

    print("Writing output file...")
    if output == _core.OutputKind.BLACK_FACES:
        mesh.save("black_faces.off")
    elif output == _core.OutputKind.SKIN:
        mesh.save("skin.off")
    else:
        mesh.save("volume.tet" if triangulate else "volume.msh")
    print("Done.")
    return 0


def main(argv: Optional[list[str]] = None) -> int:
    argv = ["mesh_generator", *(sys.argv[1:] if argv is None else argv)]
    try:
        return _run(argv)
    except (_Exit, _core.VolumeRemesherError) as e:
        sys.stdout.flush()
        sys.stderr.write(e.message if isinstance(e, _Exit) else str(e))
        return 1


if __name__ == "__main__":
    sys.exit(main())
