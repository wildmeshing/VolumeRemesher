"""The tetrahedralization's output options without command line counterparts: face_provenance
("groups", "faces", "none") and round_vertices. Each must leave the tets unchanged and report
exactly what the default output implies."""

from __future__ import annotations

from fractions import Fraction

import numpy as np
import pytest

import volume_remesher as vr

from conftest import model

# flat regions (coplanar groups), curved surfaces, touching and nested parts, an open surface
MODELS = ["cube_subdiv.off", "cube_on_cube.off", "two_cubes.off", "sphere.off",
          "double_sphere.off", "open_box.off"]


def bits_from_groups(m: vr.VolumeMesh) -> np.ndarray:
    """surface_faces as implied by triangle_provenance: each listed face marked on both tets."""
    tets = m.tets.astype(np.int64)
    # tet face 4 t + i (opposite vertex i) as a sorted triple
    keys = np.sort(np.stack([np.delete(tets, i, axis=1) for i in range(4)], axis=1), axis=2)
    keys = keys.reshape(-1, 3)
    rows = m.triangle_provenance.flat.reshape(-1, 4).astype(np.int64)
    listed = np.zeros(len(keys), bool)
    for t, *abc in rows:
        i = [k for k in range(4) if tets[t, k] not in abc]
        assert len(i) == 1
        listed[4 * t + i[0]] = True
    # mark the mate of every listed face too
    _, inverse = np.unique(keys, axis=0, return_inverse=True)
    on = np.zeros(inverse.max() + 1, bool)
    on[inverse[listed]] = True
    slot_on = on[inverse].reshape(-1, 4)
    return (slot_on * (1 << np.arange(4))).sum(axis=1).astype(np.uint8)


@pytest.mark.parametrize("name", MODELS)
@pytest.mark.parametrize("keep_all_cells", [False, True])
def test_surface_faces_are_the_grouped_faces(name, keep_all_cells):
    groups = vr.remesh(model(name), tetrahedralize=True, keep_all_cells=keep_all_cells)
    faces = vr.remesh(model(name), tetrahedralize=True, keep_all_cells=keep_all_cells,
                      face_provenance="faces")
    assert np.array_equal(groups.vertices, faces.vertices)
    assert np.array_equal(groups.tets, faces.tets)
    assert groups.surface_faces is None
    assert faces.surface_faces.dtype == np.uint8 and faces.surface_faces.shape == (len(faces.tets),)
    assert faces.surface_faces.any()
    assert np.array_equal(faces.surface_faces, bits_from_groups(groups))
    # the grouping is skipped
    assert len(faces.triangle_provenance) == 0 and len(faces.constraint_groups) == 0
    assert faces.num_coplanar_groups == 0


def test_surface_faces_of_a_boolean_and_an_embedding():
    a, b = model("cube_subdiv.off"), model("sphere.off")
    groups = vr.boolean(a, b, "union", tetrahedralize=True)
    faces = vr.boolean(a, b, "union", tetrahedralize=True, face_provenance="faces")
    assert np.array_equal(groups.tets, faces.tets)
    assert np.array_equal(faces.surface_faces, bits_from_groups(groups))
    s, t = model("Octocat.off"), model("Octocat.bg.tet")
    groups = vr.embed(s, t, tetrahedralize=True)
    faces = vr.embed(s, t, tetrahedralize=True, face_provenance="faces")
    assert np.array_equal(groups.tets, faces.tets)
    assert np.array_equal(faces.surface_faces, bits_from_groups(groups))


def test_no_face_provenance():
    groups = vr.remesh(model("cube_subdiv.off"), tetrahedralize=True, keep_all_cells=True)
    none = vr.remesh(model("cube_subdiv.off"), tetrahedralize=True, keep_all_cells=True,
                     face_provenance="none")
    assert np.array_equal(groups.tets, none.tets)
    assert len(none.triangle_provenance) == 0 and none.surface_faces is None


@pytest.mark.parametrize("name", ["double_sphere.off", "cube_on_cube.off", "112856simplified.off"])
def test_round_vertices_are_the_nearest_doubles(name):
    approx = vr.remesh(model(name), tetrahedralize=True, keep_all_cells=True, exact=True)
    rounded = vr.remesh(model(name), tetrahedralize=True, keep_all_cells=True, exact=True,
                        round_vertices=True)
    assert np.array_equal(approx.tets, rounded.tets)
    assert approx.exact_vertices == rounded.exact_vertices
    # float(Fraction) is correctly rounded
    nearest = np.array([[float(x) for x in v] for v in rounded.exact_vertices])
    assert np.array_equal(rounded.vertices, nearest)
    assert np.allclose(approx.vertices, nearest, rtol=1e-14, atol=0)


def test_option_errors():
    p = model("sphere.off")
    with pytest.raises(ValueError):
        vr.remesh(p, tetrahedralize=True, face_provenance="all")
    with pytest.raises(ValueError):
        vr.remesh(p, face_provenance="faces")
    with pytest.raises(ValueError):
        vr.remesh(p, round_vertices=True)


def test_save_skips_the_provenance_sidecars_it_did_not_compute(tmp_path):
    """With face_provenance="faces" or "none" there is no grouped provenance, so save() writes no
    .triangleprov / .groups (the tool always computes them, so its files are unaffected)."""
    import os

    import volume_remesher as vr
    from conftest import model

    p = model("sphere.off")
    vr.remesh(p, tetrahedralize=True).save(tmp_path / "groups.tet")
    vr.remesh(p, tetrahedralize=True, face_provenance="faces").save(tmp_path / "faces.tet")
    vr.remesh(p, tetrahedralize=True, face_provenance="none").save(tmp_path / "none.tet")
    names = sorted(os.listdir(tmp_path))
    assert "groups.tet.triangleprov" in names and "groups.tet.groups" in names
    assert "faces.tet" in names and "faces.tet.triangleprov" not in names and "faces.tet.groups" not in names
    assert "none.tet" in names and "none.tet.triangleprov" not in names and "none.tet.groups" not in names
