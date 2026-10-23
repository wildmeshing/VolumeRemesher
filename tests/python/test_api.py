"""The Python-facing behaviour: inputs as paths or arrays, result objects, errors, process state
(rounding mode, repeatability inside one interpreter) and the package metadata."""

from __future__ import annotations

import os
import subprocess
import sys
from fractions import Fraction

import numpy as np
import pytest

import volume_remesher as vr
from volume_remesher import _core

from conftest import model
from helpers import sha256_file


def test_version_and_exports():
    assert isinstance(vr.__version__, str)
    for name in vr.__all__:
        assert hasattr(vr, name)


def test_path_and_array_inputs_give_the_same_mesh():
    p = model("sphere.off")
    V, F = vr.read_off(p)
    assert V.dtype == np.float64 and V.shape[1] == 3
    assert F.dtype == np.uint32 and F.shape[1] == 3
    a = vr.remesh(p, tetrahedralize=True)
    b = vr.remesh((V, F), tetrahedralize=True)
    c = vr.remesh((V.tolist(), F.astype(np.int64)), tetrahedralize=True)
    assert np.array_equal(a.vertices, b.vertices) and np.array_equal(a.tets, b.tets)
    assert np.array_equal(a.vertices, c.vertices) and np.array_equal(a.tets, c.tets)
    # pathlib and str both work
    assert np.array_equal(vr.remesh(str(p)).vertices, vr.remesh(p).vertices)


def test_volume_mesh_fields():
    m = vr.remesh(model("sphere.off"), tetrahedralize=True, exact=True)
    assert isinstance(m, vr.VolumeMesh) and m.tetrahedral
    assert m.vertices.dtype == np.float64 and m.vertices.shape == (len(m.vertices), 3)
    assert m.tets.dtype == np.uint32 and m.tets.shape[1] == 4 and m.tets.max() < len(m.vertices)
    assert m.cells is None
    assert len(m.exact_vertices) == len(m.vertices)
    assert all(isinstance(x, Fraction) for x in m.exact_vertices[0])
    assert np.allclose(m.vertices, [[float(x) for x in v] for v in m.exact_vertices])
    assert len(m.triangle_provenance) == m.num_coplanar_groups
    assert m.constraint_groups.shape[1] == 2
    assert m.constraint_groups[:, 1].max() < m.num_coplanar_groups
    assert len(m.edge_ids) == len(m.edge_provenance) == 0
    assert m.point_provenance.shape == (0, 3)
    assert "VolumeMesh" in repr(m) and "tets" in repr(m)
    # arrays are read-only views that keep the result alive
    with pytest.raises(ValueError):
        m.vertices[0, 0] = 1.0


def test_polyhedral_mesh_fields():
    m = vr.remesh(model("cube_subdiv.off"))
    assert not m.tetrahedral and m.tets is None
    assert isinstance(m.cells, vr.Ragged) and len(m.cells) > 0
    assert m.cells.flat.max() < len(m.vertices)
    assert np.all(m.cells.sizes() >= 4)
    assert m.exact_vertices is None
    assert "cells" in repr(m)


def test_exact_needs_tetrahedralize():
    with pytest.raises(ValueError):
        vr.remesh(model("sphere.off"), exact=True)


def test_surface_outputs():
    p = model("sphere.off")
    s = vr.remesh(p, output="skin", tetrahedralize=True)
    assert isinstance(s, vr.SurfaceMesh)
    assert s.triangles.shape[1] == 3
    poly = vr.remesh(p, output="skin")
    assert isinstance(poly.faces, vr.Ragged)
    b = vr.remesh(p, output="black_faces")
    assert isinstance(b, vr.SurfaceMesh) and len(b.faces) > 0
    with pytest.raises(ValueError):
        vr.remesh(p, output="nope")


def test_edges_and_points_as_arrays_and_paths():
    p = model("cube_subdiv.off")
    e = model("cube_subdiv_edges.off")
    q = model("cube_subdiv_points.off")
    a = vr.remesh(p, edges=e, points=q, tetrahedralize=True)
    EV, EI = vr.read_edges_off(e)
    P = vr.read_points_off(q)
    b = vr.remesh(p, edges=(EV, EI), points=P, tetrahedralize=True)
    assert np.array_equal(a.tets, b.tets)
    assert len(a.edge_provenance) == len(EI) and np.array_equal(a.edge_ids, np.arange(len(EI)))
    assert a.point_provenance.shape == (len(P), 3)
    assert np.all(a.point_provenance[:, 2] >= 0)
    # every tracked point is an output vertex at its own position
    for pid, tet, v in a.point_provenance.tolist():
        assert np.array_equal(a.vertices[v], P[pid])
    # only edges, no surface
    c = vr.remesh(edges=e, tetrahedralize=True)
    assert len(c.edge_provenance) == len(EI)
    with pytest.raises(ValueError):
        vr.remesh()


def test_boolean_ops():
    a, b = model("cube_on_cube.off"), model("cube_subdiv.off")
    u = vr.boolean(a, b, "union", tetrahedralize=True)
    i = vr.boolean(a, b, "I", tetrahedralize=True)
    d = vr.boolean(a, b, "difference", tetrahedralize=True)
    assert len(u.tets) > len(d.tets) > 0 and len(i.tets) > 0
    with pytest.raises(ValueError):
        vr.boolean(a, b, "xor")


def test_embed_into_tet_mesh():
    s = model("Octocat.off")
    m = vr.embed(s, model("Octocat.bg.tet"), tetrahedralize=True)
    assert isinstance(m, vr.VolumeMesh) and len(m.tets) > 0
    V, T = vr.read_tet_mesh(model("Octocat.bg.tet"))
    m2 = vr.embed(s, (V, T), tetrahedralize=True)
    assert np.array_equal(m.tets, m2.tets)
    V3, T3 = vr.read_tet_mesh(model("Octocat_medit.msh"))
    assert T3.shape[1] == 4 and len(V3) > 0


def test_arrangement_fields():
    P, S = vr.read_obj_segments(model("2d/siggraph_logo.obj"))
    assert P.shape[1] == 2 and S.shape[1] == 2
    A = vr.arrange_2d((P, S), exact=True)
    assert isinstance(A, vr.Arrangement)
    assert A.vertices.shape[1] == 2 and A.triangles.shape[1] == 3
    assert len(A.segment_provenance) == len(S)
    assert A.point_provenance.shape == (len(P), 2)
    assert len(A.exact_vertices) == len(A.vertices)
    assert "Arrangement" in repr(A)


def test_errors_are_exceptions_with_the_tools_message(tmp_path):
    with pytest.raises(vr.VolumeRemesherError, match="cannot open input file"):
        vr.remesh(tmp_path / "missing.off")
    bad = tmp_path / "bad.off"
    bad.write_text("PLY\n")
    with pytest.raises(vr.VolumeRemesherError, match="different from OFF"):
        vr.remesh(bad)
    with pytest.raises(ValueError):
        vr.remesh((np.zeros((4, 2)), np.zeros((1, 3), dtype=np.uint32)))
    with pytest.raises(ValueError):
        vr.remesh((np.zeros((4, 3)), np.array([[0, 1, -1]])))
    # a degenerate input is reported, not fatal
    with pytest.raises(vr.VolumeRemesherError):
        vr.remesh((np.zeros((3, 3)), np.array([[0, 1, 2]], dtype=np.uint32)))
    assert issubclass(vr.VolumeRemesherError, Exception)


def test_rounding_mode_is_restored_and_results_repeat():
    """The exact predicates switch the FPU to round-up while they run. The pipeline runs on
    its own thread, so the caller's thread must never see that; and NFG's thread-local pools
    must leave no state that changes a later result."""
    assert _core.rounding_mode_is_nearest()
    p = model("sphere.off")
    first = vr.remesh(p, tetrahedralize=True, exact=True)
    assert _core.rounding_mode_is_nearest()
    assert 1.0 / 3.0 == 0.3333333333333333  # round-to-nearest in Python arithmetic
    for _ in range(3):
        again = vr.remesh(p, tetrahedralize=True, exact=True)
        assert np.array_equal(first.vertices, again.vertices)
        assert np.array_equal(first.tets, again.tets)
        assert first.exact_vertices == again.exact_vertices


def test_num_threads():
    """set_num_threads overrides the default thread count and 0 restores it; the result does
    not depend on the thread count."""
    default = vr.get_num_threads()
    assert default >= 1
    if hasattr(os, "sched_getaffinity") and "VOLUMEREMESHER_NUM_THREADS" not in os.environ:
        assert default == len(os.sched_getaffinity(0))
    p = model("bunny.off")
    try:
        vr.set_num_threads(1)
        assert vr.get_num_threads() == 1
        one = vr.remesh(p, tetrahedralize=True, exact=True)
        vr.set_num_threads(4)
        assert vr.get_num_threads() == 4
        four = vr.remesh(p, tetrahedralize=True, exact=True)
    finally:
        vr.set_num_threads(0)
    assert vr.get_num_threads() == default
    assert np.array_equal(one.vertices, four.vertices) and np.array_equal(one.tets, four.tets)
    assert one.exact_vertices == four.exact_vertices
    with pytest.raises(ValueError):
        vr.set_num_threads(-1)


def test_num_threads_environment_variable():
    env = dict(os.environ, VOLUMEREMESHER_NUM_THREADS="3")
    code = "import volume_remesher as vr; print(vr.get_num_threads())"
    out = subprocess.run([sys.executable, "-c", code], env=env, capture_output=True, text=True)
    assert out.returncode == 0, out.stderr
    assert out.stdout.strip() == "3"


def test_save_writes_the_sidecars(tmp_path):
    p = model("cube_subdiv.off")
    m = vr.remesh(p, edges=model("cube_subdiv_edges.off"), points=model("cube_subdiv_points.off"),
                  tetrahedralize=True, exact=True, keep_all_cells=True)
    out = tmp_path / "volume.tet"
    m.save(out)
    for suffix in ("", ".rational", ".groups", ".triangleprov", ".edgeprov", ".pointprov"):
        assert (tmp_path / ("volume.tet" + suffix)).is_file()
    # the files round-trip through the readers
    V, T = vr.read_tet_mesh(out)
    assert np.array_equal(T, m.tets)
    assert np.array_equal(V, m.vertices)
    # saving twice gives the same bytes
    m.save(tmp_path / "again.tet")
    assert sha256_file(out) == sha256_file(tmp_path / "again.tet")
    # a surface and an arrangement also save
    vr.remesh(p, output="skin").save(tmp_path / "skin.off")
    assert (tmp_path / "skin.off").read_text().startswith("OFF\n")
    vr.arrange_2d(model("2d/siggraph_logo.obj")).save(tmp_path / "tri.off")
    assert (tmp_path / "tri.off.segmentprov").is_file()


def test_ragged():
    r = vr.Ragged(np.array([0, 2, 2, 5]), np.arange(5))
    assert len(r) == 3
    assert r[0].tolist() == [0, 1] and r[1].tolist() == [] and r[-1].tolist() == [2, 3, 4]
    assert [x.tolist() for x in r[1:]] == [[], [2, 3, 4]]
    assert r.sizes().tolist() == [2, 0, 3]
    with pytest.raises(IndexError):
        r[3]
    with pytest.raises(ValueError):
        r.to_array()
    assert vr.Ragged(np.array([0, 2, 4]), np.arange(4)).to_array().tolist() == [[0, 1], [2, 3]]
    assert len(vr.Ragged(np.array([], dtype=np.uint32), np.array([], dtype=np.uint32))) == 0


def test_verbose_and_log(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    vr.remesh(model("sphere.off"), verbose=True, log=True)
    assert (tmp_path / "mesh_generator.log").is_file()
    assert "sphere.off" in (tmp_path / "mesh_generator.log").read_text()
