#include "driver.h"

#include "2d/arrangement2d.h"
#include "BSP.h"

#include <memory>

namespace vol_rem {
using NFG::ip_error;

const char* default_output_filename(const MeshOptions& opts)
{
    if (opts.output == MeshOutput::BLACK_FACES) return "black_faces.off";
    if (opts.output == MeshOutput::SKIN) return "skin.off";
    return opts.triangulate ? "volume.tet" : "volume.msh";
}

namespace {

// -e/-p inputs as the mesher wants them. Keeps no ownership: edges/points must outlive it.
bool fill_extra(const EdgeData* edges, const PointData* points, extra_features_t& extra)
{
    bool any = false;
    if (edges) {
        extra.edge_verts = edges->coords.data();
        extra.n_edge_verts = edges->num_vertices();
        extra.edge_idx = edges->edges.data();
        extra.n_edges = edges->num_edges();
        any = true;
    }
    if (points) {
        extra.point_verts = points->coords.data();
        extra.n_points = points->num_points();
        any = true;
    }
    return any;
}

// The output mesh_generator picks for these options (black faces win over skin, as in the
// command line parser), taken from the complex, which is then released.
MeshOutput extract_and_free(BSPcomplex* complex, char bool_opcode, const MeshOptions& opts)
{
    std::unique_ptr<BSPcomplex> owner(complex);
    MeshOutput out;
    out.kind = opts.output;
    if (opts.output == MeshOutput::BLACK_FACES)
        complex->extractSurface(bool_opcode, opts.triangulate, true, out.surface);
    else if (opts.output == MeshOutput::SKIN)
        complex->extractSurface(bool_opcode, opts.triangulate, false, out.surface);
    else
        complex->extractVolume(
            bool_opcode,
            opts.triangulate,
            opts.export_rational,
            out.volume,
            opts.face_provenance,
            opts.round_vertices);
    return out;
}

} // namespace

MeshOutput make_mesh(
    TriMeshData* A,
    const char* name_A,
    char bool_opcode,
    TriMeshData* B,
    const char* name_B,
    const EdgeData* edges,
    const PointData* points,
    const MeshOptions& opts)
{
    extra_features_t extra;
    const bool has_extra = fill_extra(edges, points, extra);

    BSPcomplex* complex = makePolyhedralMesh(
        name_A ? name_A : "",
        A ? A->coords.data() : nullptr,
        A ? A->num_vertices() : 0,
        A ? A->tris.data() : nullptr,
        A ? A->num_triangles() : 0,
        name_B,
        B ? B->coords.data() : nullptr,
        B ? B->num_vertices() : 0,
        B ? B->tris.data() : nullptr,
        B ? B->num_triangles() : 0,
        bool_opcode,
        false, // the TriMeshData own their memory
        opts.verbose,
        opts.logging,
        opts.fill_holes,
        opts.keep_all_cells,
        has_extra ? &extra : nullptr);
    return extract_and_free(complex, bool_opcode, opts);
}

MeshOutput embed_mesh(
    const TriMeshData& A,
    const TetMeshData& B,
    char bool_opcode,
    const EdgeData* edges,
    const PointData* points,
    const MeshOptions& opts)
{
    extra_features_t extra;
    const bool has_extra = fill_extra(edges, points, extra);

    BSPcomplex* complex = remakePolyhedralMesh(
        A.coords.data(),
        A.num_vertices(),
        A.tris.data(),
        A.num_triangles(),
        B.coords.data(),
        B.num_vertices(),
        B.tets.data(),
        B.num_tets(),
        opts.verbose,
        false,
        has_extra ? &extra : nullptr);
    return extract_and_free(complex, bool_opcode, opts);
}

vr2d::ArrangementOutput arrange_2d(
    const std::vector<double>& seg_coords,
    const std::vector<uint32_t>& seg_indexes,
    bool export_rational,
    bool verbose)
{
    vr2d::Arrangement2D arr;
    if (!vr2d::build_arrangement(seg_coords, seg_indexes, arr, verbose))
        ip_error("Failed to build the 2D arrangement\n");
    vr2d::ArrangementOutput out;
    vr2d::extract_arrangement(arr, export_rational, out);
    return out;
}

} // namespace vol_rem
