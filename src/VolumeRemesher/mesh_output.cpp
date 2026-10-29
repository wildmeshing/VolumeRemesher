#include "mesh_output.h"

#include <numerics.h> // ip_error

#include <fstream>

namespace vol_rem {

void write_surface_off(const char* filename, const SurfaceOutput& out)
{
    // Binary mode so line endings stay LF on every OS (byte-identical output).
    std::ofstream f(filename, std::ios::binary);
    if (!f) NFG::ip_error("write_surface_off: cannot open the file.\n");

    f << "OFF\n";
    f << out.num_vertices() << " ";
    f << out.num_faces() << " ";
    f << "0\n";

    for (uint64_t v = 0; v < out.num_vertices(); v++)
        f << out.coords[3 * v] << " " << out.coords[3 * v + 1] << " " << out.coords[3 * v + 2]
          << "\n";

    for (uint64_t i = 0; i < out.num_faces(); i++) {
        const uint32_t b = out.face_offsets[i], e = out.face_offsets[i + 1];
        f << (e - b);
        for (uint32_t k = b; k < e; k++) f << " " << out.face_vertices[k];
        f << "\n";
    }
    f.close();
}

void write_volume(const char* filename, const VolumeOutput& out)
{
    // Binary mode so line endings stay LF on every OS (byte-identical output).
    std::ofstream f(filename, std::ios::binary);
    if (!f) NFG::ip_error("\nBSPcomplex::[BSP.cpp]saveTetMesh: FATAL ERROR cannot open the file.\n");

    // Full double precision: the default ostream precision (6 significant figures)
    // collapses/inverts elements on fine or large-coordinate meshes.
    f.precision(17);

    const uint64_t nv = out.num_vertices();
    if (out.tetrahedral) {
        f << nv << " vertices\n";
        f << out.num_tets() << " tets\n";
        for (uint64_t v = 0; v < nv; v++)
            f << out.coords[3 * v] << " " << out.coords[3 * v + 1] << " " << out.coords[3 * v + 2]
              << "\n";
        for (uint64_t t = 0; t < out.tets.size(); t += 4)
            f << "4 " << out.tets[t] << " " << out.tets[t + 1] << " " << out.tets[t + 2] << " "
              << out.tets[t + 3] << "\n";

        // Sidecar 1: input-surface provenance, keyed by coplanar group (symmetric with the
        // edge/point sidecars below). One line per group listing the output tet faces tiling
        // it: "<group_id> <num_faces> <tet v0 v1 v2> <tet v0 v1 v2> ...", where tet is the
        // output tet index (0-based, matching the tet block above) and v0 v1 v2 are its
        // face's output vertex indices. (A face on two overlapping coplanar surfaces is
        // listed under both groups; group sizes are in the .groups sidecar.) Only when the
        // grouped provenance was computed: with FaceProvenance::FACES or NONE there is nothing
        // to list, and the tool itself always computes it.
        if (out.face_provenance == FaceProvenance::GROUPS) {
            std::ofstream pf((std::string(filename) + ".triangleprov").c_str(), std::ios::binary);
            pf << "# group_id num_faces tet v0 v1 v2 ...\n";
            pf << out.triangle_provenance.size() << "\n";
            for (uint32_t g = 0; g < out.triangle_provenance.size(); g++) {
                pf << g << " " << out.triangle_provenance[g].size();
                for (const auto& fc : out.triangle_provenance[g])
                    pf << " " << fc[0] << " " << fc[1] << " " << fc[2] << " " << fc[3];
                pf << "\n";
            }
        }

        // Sidecar 2: the coplanar group of each input triangle, so the group areas can
        // be reconstructed. "<num_constraints> <num_groups>" then, per real constraint,
        // "<input_triangle_id> <coplanar_group_id>". input_triangle_id is the triangle's
        // index in the INPUT file (degenerate/collinear input triangles are dropped by
        // the mesher, so this is not simply 0..n-1 when the input has any).
        if (out.face_provenance == FaceProvenance::GROUPS) {
            std::ofstream gf((std::string(filename) + ".groups").c_str(), std::ios::binary);
            gf << "# input_triangle_id coplanar_group_id\n";
            gf << out.constraint_groups.size() << " " << out.num_coplanar_groups << "\n";
            for (const auto& cg : out.constraint_groups) gf << cg[0] << " " << cg[1] << "\n";
        }

        // Sidecar 3: provenance of inserted edges. One line per input edge:
        // "<edge_id> <num_out_edges> <tet v0 v1> <tet v0 v1> ...", where tet is the output
        // tet index and v0 v1 are its edge's output vertex indices (matching the .tet block).
        // The listed output tet edges tile the input segment.
        if (!out.edge_provenance.empty()) {
            std::ofstream ef((std::string(filename) + ".edgeprov").c_str(), std::ios::binary);
            ef << "# edge_id num_out_edges tet v0 v1 ...\n";
            ef << out.edge_provenance.size() << "\n";
            for (const auto& ep : out.edge_provenance) {
                ef << ep.edge_id << " " << ep.out_edges.size();
                for (const auto& oe : ep.out_edges) ef << " " << oe[0] << " " << oe[1] << " " << oe[2];
                ef << "\n";
            }
        }

        // Sidecar 4: provenance of inserted points. One line per input point:
        // "<point_id> <tet> <out_vertex>", where out_vertex is the output vertex equal to the
        // point and tet is an output tet containing it (both -1 if the point did not survive).
        if (!out.point_provenance.empty()) {
            std::ofstream pf((std::string(filename) + ".pointprov").c_str(), std::ios::binary);
            pf << "# point_id tet out_vertex(-1 -1 if absent)\n";
            pf << out.point_provenance.size() << "\n";
            for (const auto& pp : out.point_provenance) {
                if (pp.out_vertex == UINT32_MAX)
                    pf << pp.point_id << " -1 -1\n";
                else
                    pf << pp.point_id << " " << pp.tet << " " << pp.out_vertex << "\n";
            }
        }

        // Optional, for exact verification: the output vertex coordinates as exact rationals,
        // same order/indexing as the .tet vertex block. Format per line: "<x> <y> <z>", each a
        // rational "[-]num[/den]" (portable decimal).
        if (!out.exact_coords.empty()) {
            std::ofstream rf((std::string(filename) + ".rational").c_str(), std::ios::binary);
            rf << nv << "\n";
            for (uint64_t v = 0; v < nv; v++)
                rf << out.exact_coords[3 * v] << " " << out.exact_coords[3 * v + 1] << " "
                   << out.exact_coords[3 * v + 2] << "\n";
        }
    } else {
        f << nv << " vertices\n";
        f << out.num_cells() << " cells\n";
        for (uint64_t v = 0; v < nv; v++)
            f << out.coords[3 * v] << " " << out.coords[3 * v + 1] << " " << out.coords[3 * v + 2]
              << "\n";
        for (uint64_t c = 0; c < out.num_cells(); c++) {
            const uint32_t b = out.cell_offsets[c], e = out.cell_offsets[c + 1];
            f << (e - b) << " ";
            for (uint32_t k = b; k < e; k++) f << out.cell_vertices[k] << " ";
            f << "\n";
        }
    }
    f.close();
}

void write_mesh_output(const char* filename, const MeshOutput& out)
{
    if (out.kind == MeshOutput::VOLUME)
        write_volume(filename, out.volume);
    else
        write_surface_off(filename, out.surface);
}

} // namespace vol_rem
