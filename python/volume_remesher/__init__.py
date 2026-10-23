"""volume_remesher: exact volume meshing of triangle soups.

The Python interface of `VolumeRemesher <https://github.com/wildmeshing/VolumeRemesher>`_.
Everything the ``mesh_generator`` command line tool does is available as a function, with the
same options and bit-identical output files:

>>> import volume_remesher as vr
>>> mesh = vr.remesh("bunny.off", tetrahedralize=True)
>>> mesh.vertices.shape, mesh.tets.shape
((..., 3), (..., 4))
>>> mesh.save("volume.tet")          # the same bytes as `mesh_generator bunny.off -t`

Inputs are paths or numpy arrays; results are frozen dataclasses of numpy arrays with a
``save()`` method. Errors raise :class:`VolumeRemesherError`.
"""

from ._api import (
    Arrangement,
    PolyhedralEmbedding,
    Ragged,
    SegmentEmbedding,
    SurfaceMesh,
    VolumeMesh,
    arrange_2d,
    arrange_2d_exact,
    boolean,
    embed,
    embed_exact,
    get_num_threads,
    read_edges_off,
    read_obj_segments,
    read_off,
    read_points_off,
    read_tet_mesh,
    remesh,
    set_num_threads,
)
from ._core import VolumeRemesherError

try:
    from ._version import version as __version__
except ImportError:  # the module was built straight from the source tree by CMake
    __version__ = "0.0.0"

__all__ = [
    "Arrangement",
    "PolyhedralEmbedding",
    "Ragged",
    "SegmentEmbedding",
    "SurfaceMesh",
    "VolumeMesh",
    "VolumeRemesherError",
    "__version__",
    "arrange_2d",
    "arrange_2d_exact",
    "boolean",
    "embed",
    "embed_exact",
    "get_num_threads",
    "read_edges_off",
    "read_obj_segments",
    "read_off",
    "read_points_off",
    "read_tet_mesh",
    "remesh",
    "set_num_threads",
]
