"""Byte-for-byte parity between mesh_generator and the Python module, on every output file
(including the provenance sidecars) and on the exit status. Needs the C++ binary from a CMake
build; skipped otherwise."""

from __future__ import annotations

import filecmp
from pathlib import Path

import pytest

from conftest import MODELS, model, run_python_cli, run_tool

pytestmark = pytest.mark.needs_cli

M = MODELS
# (argv, the files the run produces) -- every flag and input combination the tool has.
CASES = [
    (["sphere.off"], ["volume.msh"]),
    (["sphere.off", "-t"], ["volume.tet", "volume.tet.triangleprov", "volume.tet.groups"]),
    (["sphere.off", "-t", "-r"], ["volume.tet", "volume.tet.rational"]),
    (["sphere.off", "-t", "-a"], ["volume.tet"]),
    (["sphere.off", "-c", "winding"], ["volume.msh"]),
    (["bunny.off", "-t", "-c", "winding"], ["volume.tet", "volume.tet.triangleprov", "volume.tet.groups"]),
    (["open_box.off", "-f"], ["volume.msh"]),
    (["open_box.off", "-f", "-t"], ["volume.tet"]),
    (["sphere.off", "-s"], ["skin.off"]),
    (["sphere.off", "-s", "-t"], ["skin.off"]),
    (["sphere.off", "-b"], ["black_faces.off"]),
    (["sphere.off", "-b", "-t"], ["black_faces.off"]),
    (["cube_subdiv.off", "-t", "-r", "-a", "-e", "cube_subdiv_edges.off", "-p", "cube_subdiv_points.off"],
     ["volume.tet", "volume.tet.rational", "volume.tet.edgeprov", "volume.tet.pointprov",
      "volume.tet.groups", "volume.tet.triangleprov"]),
    (["-t", "-r", "-a", "-e", "sphere_edges.off"], ["volume.tet", "volume.tet.edgeprov"]),
    (["-t", "-a", "-p", "sphere_points.off"], ["volume.tet", "volume.tet.pointprov"]),
    (["two_cubes.off", "U", "cube_on_cube.off"], ["volume.msh"]),
    (["two_cubes.off", "I", "cube_on_cube.off", "-t"], ["volume.tet"]),
    (["two_cubes.off", "D", "cube_on_cube.off", "-s"], ["skin.off"]),
    (["Octocat.off", "U", "Octocat.bg.tet", "-t"], ["volume.tet", "volume.tet.triangleprov"]),
    (["Octocat.off", "U", "Octocat_medit.msh"], ["volume.msh"]),
    (["2d/siggraph_logo.obj", "-2d"], ["triangulation.off", "triangulation.off.segmentprov",
                                       "triangulation.off.pointprov"]),
    (["2d/england.obj", "-2d", "-r"], ["triangulation.off", "triangulation.off.rational"]),
]


def _resolve(argv):
    out = []
    for a in argv:
        if a.endswith((".off", ".tet", ".msh", ".obj")):
            out.append(str(model(a)))
        else:
            out.append(a)
    return out


@pytest.mark.parametrize("argv,files", CASES, ids=[" ".join(a) for a, _ in CASES])
def test_same_files_as_mesh_generator(mesh_generator, argv, files, tmp_path):
    argv = _resolve(argv)
    cli = tmp_path / "cli"
    py = tmp_path / "py"
    cli.mkdir()
    py.mkdir()
    rc = run_tool(mesh_generator, argv, cwd=cli)
    rp = run_python_cli(argv, cwd=py)
    assert rc.returncode == 0, rc.stderr
    assert rp.returncode == 0, rp.stderr
    produced_cli = sorted(p.name for p in cli.iterdir())
    produced_py = sorted(p.name for p in py.iterdir())
    assert produced_py == produced_cli
    for f in files:
        assert f in produced_cli
    for f in produced_cli:
        assert filecmp.cmp(cli / f, py / f, shallow=False), f"{f} differs"


@pytest.mark.parametrize("argv", [
    ["does_not_exist.off"],
    ["-x"],
    ["sphere.off", "-e"],
    ["sphere.off", "U"],
    ["-2d"],
    ["sphere.off", "-c"],
    ["sphere.off", "-c", "mincut"],
])
def test_same_failure_as_mesh_generator(mesh_generator, argv, tmp_path):
    argv = [str(MODELS / a) if a == "sphere.off" else a for a in argv]
    rc = run_tool(mesh_generator, argv, cwd=tmp_path)
    rp = run_python_cli(argv, cwd=tmp_path)
    assert rc.returncode != 0
    assert rp.returncode != 0
    assert rp.stderr.strip() == rc.stderr.strip()


def test_graphcut_classifier_parity(mesh_generator, tmp_path):
    """-c graphcut: the same volume.msh from the command line tool and the Python module."""
    argv = [str(MODELS / "sphere.off"), "-c", "graphcut"]
    cli = tmp_path / "cli"
    py = tmp_path / "py"
    cli.mkdir()
    py.mkdir()
    rc = run_tool(mesh_generator, argv, cwd=cli)
    rp = run_python_cli(argv, cwd=py)
    assert rc.returncode == 0 and rp.returncode == 0, rc.stderr + rp.stderr
    assert filecmp.cmp(cli / "volume.msh", py / "volume.msh", shallow=False)


def test_usage_text(mesh_generator, tmp_path):
    rc = run_tool(mesh_generator, [], cwd=tmp_path)
    rp = run_python_cli([], cwd=tmp_path)
    assert rc.returncode == 0 and rp.returncode == 0
    assert rp.stdout == rc.stdout
