"""The user-facing API over ``volume_remesher._core``.

Shape conventions: vertices are ``(n, 3)`` float64, triangles ``(m, 3)``, tetrahedra ``(m, 4)``,
edges ``(m, 2)`` and 2D points ``(n, 2)``. Index arrays are uint32 on the C++ side; any integer
array is accepted and converted.
"""

from __future__ import annotations

import os
from collections.abc import Sequence
from dataclasses import dataclass, field
from fractions import Fraction
from functools import cached_property
from typing import Optional, Union

import numpy as np

from . import _core

PathLike = Union[str, "os.PathLike[str]"]
ArrayLike = Union[np.ndarray, Sequence]
#: A surface: a path to an OFF file, or ``(vertices, triangles)`` arrays.
SurfaceLike = Union[PathLike, tuple[ArrayLike, ArrayLike]]
#: Extra edges: a path to an edge OFF file, or ``(vertices, edges)`` arrays.
EdgesLike = Union[PathLike, tuple[ArrayLike, ArrayLike]]
#: Extra points: a path to a point OFF file, or an ``(n, 3)`` array.
PointsLike = Union[PathLike, ArrayLike]
#: A tetrahedral mesh: a path to a ``.tet`` or MEDIT ``.msh``/``.mesh`` file, or
#: ``(vertices, tets)`` arrays.
TetMeshLike = Union[PathLike, tuple[ArrayLike, ArrayLike]]
#: A 2D segment soup: a path to an OBJ line file, or ``(points, segments)`` arrays.
SegmentsLike = Union[PathLike, tuple[ArrayLike, ArrayLike]]

_OUTPUT_KINDS = {
    "volume": _core.OutputKind.VOLUME,
    "skin": _core.OutputKind.SKIN,
    "black_faces": _core.OutputKind.BLACK_FACES,
}
_BOOLEAN_OPS = {
    "union": "U",
    "intersection": "I",
    "difference": "D",
    "U": "U",
    "I": "I",
    "D": "D",
}
_ABSENT = np.uint32(0xFFFFFFFF)


# --------------------------------------------------------------------------------------------
# Input coercion
# --------------------------------------------------------------------------------------------
def _is_path(x) -> bool:
    return isinstance(x, (str, bytes, os.PathLike))


def _f64(a, cols: int, what: str) -> np.ndarray:
    arr = np.ascontiguousarray(a, dtype=np.float64)
    if arr.ndim != 2 or arr.shape[1] != cols:
        raise ValueError(f"{what} must have shape (n, {cols}), got {arr.shape}")
    return arr


def _u32(a, cols: int, what: str) -> np.ndarray:
    arr = np.asarray(a)
    if arr.ndim != 2 or arr.shape[1] != cols:
        raise ValueError(f"{what} must have shape (m, {cols}), got {arr.shape}")
    if not np.issubdtype(arr.dtype, np.integer):
        raise TypeError(f"{what} must be an integer array, got dtype {arr.dtype}")
    if arr.size and (arr.min() < 0 or arr.max() > 0xFFFFFFFF):
        raise ValueError(f"{what} indices must fit in uint32")
    return np.ascontiguousarray(arr, dtype=np.uint32)


def _surface(x: SurfaceLike, verbose: bool):
    """-> ((V, F), name)."""
    if _is_path(x):
        path = os.fspath(x)
        return _core.read_off(path, verbose), path
    V, F = x
    return (_f64(V, 3, "surface vertices"), _u32(F, 3, "surface triangles")), None


def _edges(x: Optional[EdgesLike], verbose: bool):
    if x is None:
        return None
    if _is_path(x):
        return _core.read_edge_off(os.fspath(x), verbose)
    V, E = x
    return (_f64(V, 3, "edge vertices"), _u32(E, 2, "edges"))


def _points(x: Optional[PointsLike], verbose: bool):
    if x is None:
        return None
    if _is_path(x):
        return _core.read_points_off(os.fspath(x), verbose)
    return _f64(x, 3, "points")


def _tetmesh(x: TetMeshLike, verbose: bool):
    if _is_path(x):
        return read_tet_mesh(x, verbose=verbose)
    V, T = x
    return (_f64(V, 3, "tet mesh vertices"), _u32(T, 4, "tetrahedra"))


def _segments(x: SegmentsLike):
    if _is_path(x):
        return _core.read_obj_segments(os.fspath(x))
    P, S = x
    return (_f64(P, 2, "2D points"), _u32(S, 2, "segments"))


def _output_kind(output: str):
    try:
        return _OUTPUT_KINDS[output]
    except KeyError:
        raise ValueError(
            f"output must be one of {sorted(_OUTPUT_KINDS)}, got {output!r}"
        ) from None


_CLASSIFIERS = {
    "winding": _core.CellClassifier.WINDING,
    "graphcut": _core.CellClassifier.GRAPH_CUT,
}


def _classifier(classifier: str):
    try:
        return _CLASSIFIERS[classifier]
    except KeyError:
        raise ValueError(
            f"classifier must be one of {sorted(_CLASSIFIERS)}, got {classifier!r}"
        ) from None


def has_graphcut() -> bool:
    """Whether ``classifier="graphcut"`` is available in this build.

    The legacy min-cut classifier depends on third-party code (GCO, maxflow) whose license is
    not LGPL compatible, so it is compiled only when the module is built from source with
    ``-DVOLUMEREMESHER_WITH_GRAPHCUT=ON``; the wheels on PyPI do not include it.
    """
    return bool(_core.has_graphcut())


# --------------------------------------------------------------------------------------------
# Readers
# --------------------------------------------------------------------------------------------
def read_off(path: PathLike, *, verbose: bool = False) -> tuple[np.ndarray, np.ndarray]:
    """Read an OFF triangle mesh exactly as ``mesh_generator`` does.

    Returns ``(vertices, triangles)`` as ``(n, 3)`` float64 and ``(m, 3)`` uint32 arrays.
    """
    return _core.read_off(os.fspath(path), verbose)


def read_edges_off(path: PathLike, *, verbose: bool = False) -> tuple[np.ndarray, np.ndarray]:
    """Read an edge OFF file (``nv ne 0`` header, then ``2 i j`` records) -> ``(vertices, edges)``."""
    return _core.read_edge_off(os.fspath(path), verbose)


def read_points_off(path: PathLike, *, verbose: bool = False) -> np.ndarray:
    """Read a point OFF file (``np 0 0`` header) -> ``(n, 3)`` float64 array."""
    return _core.read_points_off(os.fspath(path), verbose)


def read_tet_mesh(path: PathLike, *, verbose: bool = False) -> tuple[np.ndarray, np.ndarray]:
    """Read a tetrahedral mesh -> ``(vertices, tets)``.

    ``.tet`` is the format :meth:`VolumeMesh.save` writes; anything else is read as MEDIT
    (``.msh``/``.mesh``), like ``mesh_generator`` does for its second input.
    """
    path = os.fspath(path)
    if path.endswith(".tet"):
        return _core.read_tet(path, verbose)
    return _core.read_medit(path, verbose)


def read_obj_segments(path: PathLike) -> tuple[np.ndarray, np.ndarray]:
    """Read an OBJ line soup (``v x y z`` + ``l i j``) -> ``(points (n, 2), segments (m, 2))``."""
    return _core.read_obj_segments(os.fspath(path))


# --------------------------------------------------------------------------------------------
# Results
# --------------------------------------------------------------------------------------------
class Ragged(Sequence):
    """A sequence of integer arrays of varying length, stored flat.

    ``r[i]`` is the ``i``-th array (a view). ``r.flat`` is the concatenation and ``r.offsets``
    (length ``len(r) + 1``) delimits the pieces, so ``r[i] is r.flat[r.offsets[i]:r.offsets[i+1]]``.
    """

    __slots__ = ("offsets", "flat")

    def __init__(self, offsets: np.ndarray, flat: np.ndarray):
        self.offsets = np.asarray(offsets)
        self.flat = np.asarray(flat)

    def __len__(self) -> int:
        return max(len(self.offsets) - 1, 0)

    def __getitem__(self, i):
        if isinstance(i, slice):
            return [self[k] for k in range(*i.indices(len(self)))]
        n = len(self)
        if i < 0:
            i += n
        if not 0 <= i < n:
            raise IndexError(i)
        return self.flat[self.offsets[i] : self.offsets[i + 1]]

    def sizes(self) -> np.ndarray:
        """The length of every piece."""
        return np.diff(self.offsets)

    def to_array(self) -> np.ndarray:
        """All pieces as one 2-D array; raises ``ValueError`` unless they have equal length."""
        n = len(self)
        if n == 0:
            return self.flat.reshape((0,) + self.flat.shape[1:])
        sizes = self.sizes()
        if not np.all(sizes == sizes[0]):
            raise ValueError("pieces have different lengths")
        return self.flat.reshape((n, int(sizes[0])) + self.flat.shape[1:])

    def __repr__(self) -> str:
        return f"Ragged({len(self)} pieces, {len(self.flat)} items)"


def _fractions(strings: list[str], dim: int) -> list[tuple[Fraction, ...]]:
    it = iter(Fraction(s) for s in strings)
    return [tuple(next(it) for _ in range(dim)) for _ in range(len(strings) // dim)]


def _signed(a: np.ndarray) -> np.ndarray:
    """uint32 with UINT32_MAX as 'absent' -> int64 with -1."""
    out = a.astype(np.int64)
    out[a == _ABSENT] = -1
    return out


@dataclass(frozen=True, repr=False)
class VolumeMesh:
    """A volume mesh: polyhedral cells, or tetrahedra with provenance.

    Produced by :func:`remesh`, :func:`boolean` and :func:`embed` with ``output="volume"``.
    """

    #: ``(n, 3)`` float64. The double approximation of the exact vertex positions, the same
    #: numbers the ``.tet``/``.msh`` file holds; with ``round_vertices=True`` the nearest doubles.
    vertices: np.ndarray
    #: ``(m, 4)`` uint32, or ``None`` for a polyhedral mesh.
    tets: Optional[np.ndarray]
    #: The polyhedral cells as vertex lists, or ``None`` for a tetrahedral mesh.
    cells: Optional[Ragged]
    #: Per coplanar group of input triangles, the ``(k, 4)`` ``{tet, v0, v1, v2}`` output faces
    #: tiling it. Tetrahedral meshes with ``face_provenance="groups"`` only.
    triangle_provenance: Ragged
    #: ``(c, 2)`` ``{input_triangle_id, coplanar_group_id}`` per non-degenerate input triangle.
    constraint_groups: np.ndarray
    num_coplanar_groups: int
    #: ``edge_ids[i]`` is the input edge whose output edges are ``edge_provenance[i]``
    #: (``(k, 3)`` ``{tet, v0, v1}``).
    edge_ids: np.ndarray
    edge_provenance: Ragged
    #: ``(p, 3)`` int64 ``{point_id, tet, vertex}`` per input point; ``-1`` where it did not
    #: survive.
    point_provenance: np.ndarray
    #: ``(m,)`` uint8 with ``face_provenance="faces"``, else ``None``: bit ``i`` of entry ``t``
    #: is set iff the face of tet ``t`` opposite its vertex ``i`` tiles input triangles (the
    #: faces ``triangle_provenance`` would list, marked on both their tets).
    surface_faces: Optional[np.ndarray]
    _core_mesh: _core.Mesh = field(repr=False, compare=False)
    _exact: list[str] = field(default_factory=list, repr=False, compare=False)

    @cached_property
    def exact_vertices(self) -> Optional[list[tuple[Fraction, Fraction, Fraction]]]:
        """The exact rational vertex positions, or ``None`` unless ``exact=True`` was asked."""
        if not self._exact:
            return None
        return _fractions(self._exact, 3)

    @property
    def tetrahedral(self) -> bool:
        return self.tets is not None

    def save(self, path: PathLike) -> None:
        """Write the mesh exactly as ``mesh_generator`` would.

        A tetrahedral mesh gets the ``.tet`` format plus the ``.triangleprov``, ``.groups``,
        ``.edgeprov`` and ``.pointprov`` sidecars (and ``.rational`` when built with
        ``exact=True``); a polyhedral mesh gets the ``.msh`` format. Byte-identical to the tool.
        With ``face_provenance="faces"`` or ``"none"`` the ``.triangleprov`` and ``.groups``
        sidecars are not written (the tool always computes the grouped provenance).
        """
        self._core_mesh.save(os.fspath(path))

    def __repr__(self) -> str:
        if self.tets is not None:
            return f"VolumeMesh({len(self.vertices)} vertices, {len(self.tets)} tets)"
        return f"VolumeMesh({len(self.vertices)} vertices, {len(self.cells)} cells)"


@dataclass(frozen=True, repr=False)
class SurfaceMesh:
    """A polygon mesh: the bounding skin (``output="skin"``) or the subdivided input constraints
    (``output="black_faces"``)."""

    #: ``(n, 3)`` float64.
    vertices: np.ndarray
    #: The faces as vertex lists (triangles when built with ``tetrahedralize=True``).
    faces: Ragged
    _core_mesh: _core.Mesh = field(repr=False, compare=False)

    @property
    def triangles(self) -> np.ndarray:
        """``(m, 3)`` uint32; raises ``ValueError`` if a face is not a triangle."""
        return self.faces.to_array()

    def save(self, path: PathLike) -> None:
        """Write the OFF file exactly as ``mesh_generator -s`` / ``-b`` would."""
        self._core_mesh.save(os.fspath(path))

    def __repr__(self) -> str:
        return f"SurfaceMesh({len(self.vertices)} vertices, {len(self.faces)} faces)"


@dataclass(frozen=True, repr=False)
class Arrangement:
    """The arrangement of a 2D segment soup as a triangulation, from :func:`arrange_2d`."""

    #: ``(n, 2)`` float64.
    vertices: np.ndarray
    #: ``(m, 3)`` uint32, counterclockwise.
    triangles: np.ndarray
    #: Per input segment, the ``(k, 3)`` ``{triangle, v0, v1}`` output edges tiling it, ordered
    #: from its first endpoint; empty for a dropped (zero-length) segment.
    segment_provenance: Ragged
    #: ``(p, 2)`` int64 ``{triangle, vertex}`` per input point; ``-1`` where it did not survive.
    point_provenance: np.ndarray
    _core_arr: _core.Arrangement = field(repr=False, compare=False)
    _exact: list[str] = field(default_factory=list, repr=False, compare=False)

    @cached_property
    def exact_vertices(self) -> Optional[list[tuple[Fraction, Fraction]]]:
        """The exact rational vertex positions, or ``None`` unless ``exact=True`` was asked."""
        if not self._exact:
            return None
        return _fractions(self._exact, 2)

    def save(self, path: PathLike) -> None:
        """Write the OFF triangulation at ``path`` with its ``.segmentprov``/``.pointprov``
        (and ``.rational`` with ``exact=True``), byte-identical to ``mesh_generator -2d``."""
        self._core_arr.save(os.fspath(path))

    def __repr__(self) -> str:
        return f"Arrangement({len(self.vertices)} vertices, {len(self.triangles)} triangles)"


@dataclass(frozen=True, repr=False)
class PolyhedralEmbedding:
    """The exact output of embedding triangles (and edges/points) into a tetrahedral mesh, from
    :func:`embed_exact`: the C++ ``embed_tri_in_poly_mesh`` API.

    Vertices are exact rationals in lowest terms. ``facets``/``cells`` describe the polyhedral
    complex whose cells are (weakly) convex; ``tets`` is its tetrahedralization.
    """

    #: Polygonal facets, each a list of vertex indices.
    facets: Ragged
    #: Polyhedral cells, each a list of facet indices.
    cells: Ragged
    #: ``(m, 4)`` uint32 tetrahedra, positively oriented.
    tets: np.ndarray
    #: ``(m,)`` the cell each tet was carved from.
    tets_parent: np.ndarray
    #: Per tet, the indices (into ``facets``) of its parent cell's facets it lies on.
    tets_parent_faces: Ragged
    #: Facets that overlap the input triangles.
    facets_on_input: np.ndarray
    #: ``(c,)`` bool, cells having a facet on the input.
    cells_with_faces_on_input: np.ndarray
    #: Per coplanar group, the ``(k, 4)`` ``{tet, v0, v1, v2}`` faces tiling it.
    triangle_provenance: Ragged
    #: ``(t,)`` int64, the coplanar group of each input triangle, ``-1`` for a dropped one.
    triangle_group: np.ndarray
    #: Per input edge, the ``(k, 3)`` ``{tet, v0, v1}`` output edges tiling it.
    edge_provenance: Ragged
    #: ``(p, 2)`` int64 ``{tet, vertex}`` per input point, ``-1`` if it did not survive.
    point_provenance: np.ndarray
    #: The raw ``'[-]num[/den]'`` coordinate strings, 3 per vertex.
    exact_coords: list[str] = field(repr=False, compare=False)

    @property
    def num_vertices(self) -> int:
        return len(self.exact_coords) // 3

    @cached_property
    def vertices(self) -> list[tuple[Fraction, Fraction, Fraction]]:
        """``(n, 3)`` tuples of :class:`fractions.Fraction`, exact and in lowest terms. Built on
        first access: on a mesh of millions of vertices this takes a while and a lot of memory."""
        return _fractions(self.exact_coords, 3)

    def __repr__(self) -> str:
        return (f"PolyhedralEmbedding({self.num_vertices} vertices, {len(self.cells)} cells, "
                f"{len(self.tets)} tets)")


@dataclass(frozen=True, repr=False)
class SegmentEmbedding:
    """The exact 2D arrangement from :func:`arrange_2d_exact`: the C++ ``embed_seg_in_tri_mesh``
    API. Like :class:`Arrangement` but with exact vertices in lowest terms and no file output."""

    #: ``(m, 3)`` uint32, counterclockwise.
    triangles: np.ndarray
    #: Per input segment, the ``(k, 3)`` ``{triangle, v0, v1}`` output edges tiling it.
    segment_provenance: Ragged
    #: ``(p, 2)`` int64 ``{triangle, vertex}`` per input point, ``-1`` if it did not survive.
    point_provenance: np.ndarray
    #: The raw ``'[-]num[/den]'`` coordinate strings, 2 per vertex.
    exact_coords: list[str] = field(repr=False, compare=False)

    @property
    def num_vertices(self) -> int:
        return len(self.exact_coords) // 2

    @cached_property
    def vertices(self) -> list[tuple[Fraction, Fraction]]:
        """``(n, 2)`` tuples of :class:`fractions.Fraction`, exact and in lowest terms. Built on
        first access."""
        return _fractions(self.exact_coords, 2)

    def __repr__(self) -> str:
        return f"SegmentEmbedding({self.num_vertices} vertices, {len(self.triangles)} triangles)"


_FACE_PROVENANCE = ("groups", "faces", "none")


def _check_tet_options(tetrahedralize: bool, exact: bool, face_provenance: str, round_vertices: bool):
    if exact and not tetrahedralize:
        raise ValueError("exact=True needs tetrahedralize=True")
    if face_provenance not in _FACE_PROVENANCE:
        raise ValueError(f"face_provenance must be one of {_FACE_PROVENANCE}, got {face_provenance!r}")
    if (face_provenance != "groups" or round_vertices) and not tetrahedralize:
        raise ValueError("face_provenance and round_vertices need tetrahedralize=True")


def _wrap_mesh(m: _core.Mesh, face_provenance: str = "groups") -> Union[VolumeMesh, SurfaceMesh]:
    if m.kind != _core.OutputKind.VOLUME:
        return SurfaceMesh(vertices=m.coords, faces=Ragged(m.face_offsets, m.face_vertices), _core_mesh=m)
    tetrahedral = m.tetrahedral
    return VolumeMesh(
        vertices=m.coords,
        tets=m.tets if tetrahedral else None,
        cells=None if tetrahedral else Ragged(m.cell_offsets, m.cell_vertices),
        triangle_provenance=Ragged(m.triangle_provenance_offsets, m.triangle_provenance),
        constraint_groups=m.constraint_groups,
        num_coplanar_groups=int(m.num_coplanar_groups),
        edge_ids=m.edge_ids,
        edge_provenance=Ragged(m.edge_provenance_offsets, m.edge_provenance),
        point_provenance=_signed(m.point_provenance),
        surface_faces=m.surface_faces if tetrahedral and face_provenance == "faces" else None,
        _core_mesh=m,
        _exact=list(m.exact_coords),
    )


# --------------------------------------------------------------------------------------------
# Pipeline
# --------------------------------------------------------------------------------------------
def set_num_threads(n: int) -> None:
    """Run the parallel passes on ``n`` threads; ``0`` restores the default (:func:`get_num_threads`)."""
    if n < 0:
        raise ValueError(f"n must be >= 0, got {n}")
    _core.set_num_threads(n)


def get_num_threads() -> int:
    """The threads of the parallel passes: :func:`set_num_threads`, else the
    ``VOLUMEREMESHER_NUM_THREADS`` environment variable, else the CPUs this process may run on
    (its affinity mask on Linux, e.g. under ``taskset`` or a batch scheduler)."""
    return _core.get_num_threads()


def remesh(
    surface: Optional[SurfaceLike] = None,
    *,
    edges: Optional[EdgesLike] = None,
    points: Optional[PointsLike] = None,
    tetrahedralize: bool = False,
    keep_all_cells: bool = False,
    fill_holes: bool = False,
    exact: bool = False,
    face_provenance: str = "groups",
    round_vertices: bool = False,
    output: str = "volume",
    classifier: str = "winding",
    verbose: bool = False,
    log: bool = False,
    name: Optional[str] = None,
) -> Union[VolumeMesh, SurfaceMesh]:
    """Mesh the volume enclosed by a surface, resolving self-intersections and repairing it.

    Equivalent to ``mesh_generator surface.off [-e edges] [-p points] [-t] [-a] [-f] [-r]
    [-s | -b] [-c classifier]``.

    Args:
        surface: an OFF path or ``(vertices, triangles)``. Optional if ``edges`` or ``points``
            are given, in which case the domain is their bounding box.
        edges: extra edges to force into the output (implies ``keep_all_cells``).
        points: extra points to force into the output (implies ``keep_all_cells``).
        tetrahedralize: ``-t``. Tetrahedral instead of polyhedral output; for ``output="skin"``
            or ``"black_faces"``, triangulated faces.
        keep_all_cells: ``-a``. Skip the in/out classification and keep the whole domain, which
            makes the surface tracking exact.
        fill_holes: ``-f``. Cap open boundaries (solid-output path only).
        exact: ``-r``. Also compute the exact rational vertex coordinates
            (:attr:`VolumeMesh.exact_vertices`, and the ``.rational`` sidecar on save). Needs
            ``tetrahedralize``.
        face_provenance: with ``tetrahedralize``, how the output faces on the input surface are
            reported: ``"groups"`` (default) per coplanar group of input triangles
            (:attr:`VolumeMesh.triangle_provenance`, ``constraint_groups``); ``"faces"`` only
            which tet faces they are (:attr:`VolumeMesh.surface_faces`), which skips the
            grouping; ``"none"`` not at all. No command line counterpart.
        round_vertices: with ``tetrahedralize``, :attr:`VolumeMesh.vertices` are the nearest
            doubles of the exact vertex positions instead of approximations (which can be an
            ulp or so off). No command line counterpart.
        output: ``"volume"`` (default), ``"skin"`` (``-s``) or ``"black_faces"`` (``-b``).
        classifier: ``-c``. How cells are decided inside or outside: ``"winding"`` (default),
            the winding number of the surface, exact for a closed and consistently oriented
            input and the generalized winding number otherwise; or ``"graphcut"``, the legacy
            area-weighted min-cut, available only when the module was built with
            ``VOLUMEREMESHER_WITH_GRAPHCUT`` (see :func:`has_graphcut`; the PyPI wheels were
            not, its third-party code is not LGPL compatible).
        verbose: ``-v``. Progress on the C stdout.
        log: ``-l``. Append a timing line to ``mesh_generator.log`` in the working directory.
        name: how the input is called in the log; defaults to its path, or ``"<array>"``.

    Returns:
        :class:`VolumeMesh` for ``output="volume"``, else :class:`SurfaceMesh`.

    Raises:
        VolumeRemesherError: on a malformed input or a failure inside the pipeline.
    """
    if surface is None and edges is None and points is None:
        raise ValueError("provide a surface and/or edges and/or points")
    _check_tet_options(tetrahedralize, exact, face_provenance, round_vertices)
    A, path = (None, None) if surface is None else _surface(surface, verbose)
    core = _core.make_mesh(
        A, name or path or ("<array>" if A is not None else None), "0", None, None,
        _edges(edges, verbose), _points(points, verbose),
        verbose, log, tetrahedralize, exact, fill_holes, keep_all_cells, _output_kind(output),
        _classifier(classifier), face_provenance, round_vertices,
    )
    return _wrap_mesh(core, face_provenance)


def boolean(
    a: SurfaceLike,
    b: SurfaceLike,
    op: str,
    *,
    edges: Optional[EdgesLike] = None,
    points: Optional[PointsLike] = None,
    tetrahedralize: bool = False,
    keep_all_cells: bool = False,
    fill_holes: bool = False,
    exact: bool = False,
    face_provenance: str = "groups",
    round_vertices: bool = False,
    output: str = "volume",
    classifier: str = "winding",
    verbose: bool = False,
    log: bool = False,
    names: Optional[tuple[str, str]] = None,
) -> Union[VolumeMesh, SurfaceMesh]:
    """Mesh the boolean of two surfaces: ``mesh_generator a.off {U|I|D} b.off ...``.

    ``op`` is ``"union"``, ``"intersection"`` or ``"difference"`` (``a`` minus ``b``), or the
    tool's letters ``"U"``, ``"I"``, ``"D"``. The other arguments are those of :func:`remesh`.
    """
    try:
        opcode = _BOOLEAN_OPS[op]
    except KeyError:
        raise ValueError(f"op must be union/intersection/difference (or U/I/D), got {op!r}") from None
    _check_tet_options(tetrahedralize, exact, face_provenance, round_vertices)
    A, path_a = _surface(a, verbose)
    B, path_b = _surface(b, verbose)
    name_a, name_b = names if names else (path_a or "<array A>", path_b or "<array B>")
    core = _core.make_mesh(
        A, name_a, opcode, B, name_b,
        _edges(edges, verbose), _points(points, verbose),
        verbose, log, tetrahedralize, exact, fill_holes, keep_all_cells, _output_kind(output),
        _classifier(classifier), face_provenance, round_vertices,
    )
    return _wrap_mesh(core, face_provenance)


def embed(
    surface: SurfaceLike,
    tet_mesh: TetMeshLike,
    *,
    edges: Optional[EdgesLike] = None,
    points: Optional[PointsLike] = None,
    tetrahedralize: bool = False,
    exact: bool = False,
    face_provenance: str = "groups",
    round_vertices: bool = False,
    output: str = "volume",
    classifier: str = "winding",
    verbose: bool = False,
) -> Union[VolumeMesh, SurfaceMesh]:
    """Embed a surface (and optional edges/points) into an existing tetrahedral mesh, cutting
    its tetrahedra so every input triangle is a union of output faces.

    Equivalent to ``mesh_generator surface.off U mesh.tet [-e] [-p] [-t] [-r] [-s | -b]
    [-c classifier]`` (the boolean letter is a separator the tool requires and does not use).
    The tet mesh must be convex and every input vertex must lie strictly inside it, or coincide
    with one of its vertices. As with the tool, the cut cells are classified against the
    surface (``classifier``, see :func:`remesh`) and only the inside ones are returned. The other arguments are those of :func:`remesh`.
    """
    _check_tet_options(tetrahedralize, exact, face_provenance, round_vertices)
    A, _ = _surface(surface, verbose)
    T = _tetmesh(tet_mesh, verbose)
    core = _core.embed_mesh(
        A, T, "0", _edges(edges, verbose), _points(points, verbose),
        verbose, tetrahedralize, exact, _output_kind(output), _classifier(classifier),
        face_provenance, round_vertices,
    )
    return _wrap_mesh(core, face_provenance)


def arrange_2d(
    segments: SegmentsLike,
    *,
    exact: bool = False,
    verbose: bool = False,
) -> Arrangement:
    """The arrangement of a soup of 2D segments as a triangulation: ``mesh_generator -2d``.

    ``segments`` is an OBJ line-soup path or ``(points (n, 2), segments (m, 2))``. Segments may
    cross, overlap, repeat or be degenerate; every input segment is exactly a union of output
    edges (see :attr:`Arrangement.segment_provenance`). The triangulation covers the bounding
    box expanded by 10%.
    """
    P, S = _segments(segments)
    core = _core.arrange_2d(P, S, exact, verbose)
    return Arrangement(
        vertices=core.coords,
        triangles=core.triangles,
        segment_provenance=Ragged(core.segment_provenance_offsets, core.segment_provenance),
        point_provenance=_signed(core.point_provenance),
        _core_arr=core,
        _exact=list(core.exact_coords),
    )


def embed_exact(
    surface: SurfaceLike,
    tet_mesh: TetMeshLike,
    *,
    edges: Optional[EdgesLike] = None,
    points: Optional[PointsLike] = None,
    verbose: bool = False,
) -> PolyhedralEmbedding:
    """Embed triangles (and optional edges/points) into a tetrahedral mesh with exact output.

    The C++ ``embed_tri_in_poly_mesh`` API: cuts the tetrahedra so that every input triangle is
    a union of output facets and every cell is (weakly) convex, and returns the polyhedral
    complex, its tetrahedralization and the provenance with exact rational coordinates in
    lowest terms. The tet mesh must be convex; every vertex of the surface and of the extra
    features must lie strictly inside it or coincide with one of its vertices.

    Raises:
        VolumeRemesherError: for a vertex outside or on the boundary of the tet mesh, a
            zero-length edge, or any other failure.
    """
    A, _ = _surface(surface, verbose)
    T = _tetmesh(tet_mesh, verbose)
    r = _core.embed_polyhedral(A, T, _edges(edges, verbose), _points(points, verbose), verbose)
    exact = list(r.exact_coords)
    return PolyhedralEmbedding(
        facets=Ragged(r.facet_offsets, r.facet_vertices),
        cells=Ragged(r.cell_offsets, r.cell_facets),
        tets=r.tets,
        tets_parent=r.tets_parent,
        tets_parent_faces=Ragged(r.tets_parent_faces_offsets, r.tets_parent_faces),
        facets_on_input=r.facets_on_input,
        cells_with_faces_on_input=r.cells_with_faces_on_input.astype(bool),
        triangle_provenance=Ragged(r.triangle_provenance_offsets, r.triangle_provenance),
        triangle_group=_signed(r.triangle_group),
        edge_provenance=Ragged(r.edge_provenance_offsets, r.edge_provenance),
        point_provenance=_signed(r.point_provenance),
        exact_coords=exact,
    )


def arrange_2d_exact(segments: SegmentsLike, *, verbose: bool = False) -> SegmentEmbedding:
    """The arrangement of a 2D segment soup with exact vertices in lowest terms.

    The C++ ``embed_seg_in_tri_mesh`` API; see :func:`arrange_2d` for the file-writing form.
    """
    P, S = _segments(segments)
    r = _core.embed_segments(P, S, verbose)
    exact = list(r.exact_coords)
    return SegmentEmbedding(
        triangles=r.triangles,
        segment_provenance=Ragged(r.segment_provenance_offsets, r.segment_provenance),
        point_provenance=_signed(r.point_provenance),
        exact_coords=exact,
    )
