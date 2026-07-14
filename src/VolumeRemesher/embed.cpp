#include "embed.h"
#include <vector>
#include "BSP.h"

namespace vol_rem {
void embed_tri_in_poly_mesh(
    const std::vector<double>& tri_vrt_coords,
    const std::vector<uint32_t>& triangle_indexes,
    const std::vector<double>& tet_vrt_coords,
    const std::vector<uint32_t>& tet_indexes,
    std::vector<bigrational>& out_vrt_coords,
    std::vector<uint32_t>& out_poly_vindexes,
    std::vector<uint32_t>& out_cell_findexes,
    std::vector<std::array<uint32_t, 4>>& out_tets,
    std::vector<uint32_t>& final_tets_parent,
    std::vector<uint32_t>& facets_on_input,
    std::vector<bool>& cells_with_faces_on_input,
    std::vector<std::vector<uint32_t>>& final_tets_parent_faces,
    const std::vector<double>& edge_vrt_coords,
    const std::vector<uint32_t>& edge_indexes,
    const std::vector<double>& point_coords,
    std::vector<std::vector<std::array<uint32_t, 2>>>& out_edge_provenance,
    std::vector<uint32_t>& out_point_provenance,
    bool verbose)
{
    // Optional extra edges/points to force into the embedding.
    extra_features_t extra;
    extra.edge_verts = edge_vrt_coords.data();
    extra.n_edge_verts = (uint32_t)(edge_vrt_coords.size() / 3);
    extra.edge_idx = edge_indexes.data();
    extra.n_edges = (uint32_t)(edge_indexes.size() / 2);
    extra.point_verts = point_coords.data();
    extra.n_points = (uint32_t)(point_coords.size() / 3);
    const bool has_extra = extra.n_edges > 0 || extra.n_points > 0;

    // Make a conformal polyhedralization
    BSPcomplex* complex = remakePolyhedralMesh(
        tri_vrt_coords.data(),
        (uint32_t)tri_vrt_coords.size() / 3,
        triangle_indexes.data(),
        (uint32_t)triangle_indexes.size() / 3,
        tet_vrt_coords.data(),
        (uint32_t)tet_vrt_coords.size() / 3,
        tet_indexes.data(),
        (uint32_t)tet_indexes.size() / 4,
        verbose,
        true,
        has_extra ? &extra : nullptr);

    // Triangulate every (convex) BSP face so the output facets are all triangles,
    // then tetrahedralize the whole complex. keep_all_cells=true is essential
    // here: this embedding keeps the entire background domain (both sides of the
    // inserted surface), unlike the boolean path which only keeps region A.
    // TODOfix: triangulateFace can cause degenerate triangles!!!
    for (size_t f_id = 0; f_id < complex->faces.size(); f_id++) {
        complex->triangulateFace(f_id);
    }
    complex->makeTetrahedra(verbose, /*keep_all_cells=*/true);

    if (verbose) printf("Producing vertices...\n");
    // Get exact vertex coordinates
    out_vrt_coords.resize(complex->vertices.size() * 3);
    for (uint64_t v_id = 0; v_id < complex->vertices.size(); v_id++) {
        if (!complex->vertices[v_id]->getExactXYZCoordinates(
                out_vrt_coords[v_id * 3],
                out_vrt_coords[v_id * 3 + 1],
                out_vrt_coords[v_id * 3 + 2]))
            ip_error(
                "embed_tri_in_poly_mesh: could not compute exact coordinates. Should not "
                "happen!\n");
    }

    if (verbose) printf("Producing facets...\n");
    // Get facets
    for (size_t f_id = 0; f_id < complex->faces.size(); f_id++) {
        BSPface& face = complex->faces[f_id];
        std::vector<uint32_t> face_vrts(face.edges.size(), 0);
        complex->list_faceVertices(face, face_vrts);
        out_poly_vindexes.push_back((uint32_t)face_vrts.size());
        for (uint32_t cvi : face_vrts) {
            out_poly_vindexes.push_back(cvi);
        }
        if (face.colour == BLACK_A) {
            facets_on_input.push_back((uint32_t)f_id);
        }
    }

    if (verbose) printf("Producing cells...\n");
    // Get polyhedra
    for (uint64_t c_id = 0; c_id < complex->cells.size(); c_id++) {
        BSPcell& cell = complex->cells[c_id];
        out_cell_findexes.push_back((uint32_t)cell.faces.size());
        for (uint64_t cfi : cell.faces) {
            out_cell_findexes.push_back((uint32_t)cfi);
        }
    }

    if (verbose) printf("Producing tets...\n");
    // Get tets
    assert(complex->final_tets.size() % 4 == 0);
    out_tets.resize(complex->final_tets.size() / 4);
    for (uint64_t t_id = 0; t_id < complex->final_tets.size(); t_id += 4) {
        for (uint64_t i = 0; i < 4; ++i) {
            out_tets[t_id / 4][i] = complex->final_tets[t_id + i];
        }
    }

    // Surface-tracking metadata, so the caller can carry the input-surface tags
    // from the polygonal faces onto the output tets.

    // Parent BSP cell of each output tet (recorded by makeTetrahedra).
    final_tets_parent.assign(
        complex->final_tets_parent_cell.begin(), complex->final_tets_parent_cell.end());

    // Per-cell flag: does the cell touch the input surface? A face is on the
    // input surface iff its colour is BLACK_A.
    cells_with_faces_on_input.assign(complex->cells.size(), false);
    for (uint64_t c_id = 0; c_id < complex->cells.size(); c_id++) {
        for (uint64_t f_id : complex->cells[c_id].faces) {
            if (complex->faces[f_id].colour == BLACK_A) {
                cells_with_faces_on_input[c_id] = true;
                break;
            }
        }
    }

    // For each output tet, which of its parent cell's (triangular) faces actually
    // bound it: a triangle face is a face of the tet iff its three vertices are a
    // subset of the tet's four vertices (a tetrahedron's four faces are exactly
    // the four vertex-triples of its vertices). Faces introduced internally by
    // the decomposition (barycenter / cone apex) are not BSP faces and are thus
    // correctly excluded. Only computed for cells that touch the input surface;
    // for the rest the caller leaves every tet face untagged anyway.
    final_tets_parent_faces.assign(out_tets.size(), {});
    for (uint64_t t_id = 0; t_id < out_tets.size(); t_id++) {
        const uint32_t c_id = final_tets_parent[t_id];
        if (!cells_with_faces_on_input[c_id]) continue;
        const std::array<uint32_t, 4>& tet = out_tets[t_id];
        for (uint64_t f_id : complex->cells[c_id].faces) {
            std::vector<uint32_t> fv(complex->faces[f_id].edges.size(), 0);
            complex->list_faceVertices(complex->faces[f_id], fv);
            if (fv.size() != 3) continue; // only triangles bound a tet face
            bool subset = true;
            for (uint32_t w : fv) {
                if (w != tet[0] && w != tet[1] && w != tet[2] && w != tet[3]) {
                    subset = false;
                    break;
                }
            }
            if (subset) final_tets_parent_faces[t_id].push_back((uint32_t)f_id);
        }
    }

    // Provenance of the inserted edges/points. The indices are into 'vertices'/'out_tets'
    // (this path emits every complex vertex, so no remapping is needed).
    out_edge_provenance.assign(complex->edge_provenance.size(), {});
    for (const auto& ep : complex->edge_provenance)
        if (ep.edge_id < out_edge_provenance.size()) out_edge_provenance[ep.edge_id] = ep.out_edges;
    out_point_provenance.assign(complex->point_provenance.size(), UINT32_MAX);
    for (const auto& pp : complex->point_provenance)
        if (pp.point_id < out_point_provenance.size())
            out_point_provenance[pp.point_id] = pp.out_vertex;

    if (verbose) printf("Done\n");
}
} // namespace vol_rem