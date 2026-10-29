// In-memory form of everything mesh_generator writes, and the writers that turn it into the
// files. BSPcomplex::extract* fills these (BSP.cpp); write_* serializes them byte-for-byte as
// the old save* members did, so a file written from a snapshot is identical to one written
// straight from the complex. The Python module hands the snapshot to the caller as arrays and
// writes files from it on request, through these same writers, which is what makes its output
// bit-identical to the command line tool's.
//
// Everything here is plain data (doubles, integers, decimal strings): a snapshot may be built on
// one thread and read on another, which NFG's pool-allocated exact numbers do not allow (see
// parallel.h), so exact coordinates are carried as their portable decimal "[-]num[/den]" form.

#ifndef VOLUMEREMESHER_MESH_OUTPUT_H
#define VOLUMEREMESHER_MESH_OUTPUT_H

#include <array>
#include <cstdint>
#include <string>
#include <vector>

namespace vol_rem {

// A polygon soup written as OFF: the bounding skin (-s) or the subdivided constraints (-b).
struct SurfaceOutput
{
    std::vector<double> coords; // 3 per vertex
    std::vector<uint32_t> face_offsets; // num_faces + 1 entries into face_vertices
    std::vector<uint32_t> face_vertices;
    uint64_t num_vertices() const { return coords.size() / 3; }
    uint64_t num_faces() const { return face_offsets.empty() ? 0 : face_offsets.size() - 1; }
};

// How a tetrahedralization reports its faces on the input surface.
enum class FaceProvenance {
    GROUPS, // triangle_provenance per coplanar group of input triangles, and constraint_groups
    FACES, // only surface_faces: which tet faces tile input triangles, without the grouping
    NONE // neither
};

// The volume: polyhedral cells (volume.msh) or tetrahedra plus provenance (volume.tet and its
// sidecar files).
struct VolumeOutput
{
    bool tetrahedral = false;
    std::vector<double> coords; // 3 per vertex

    // tetrahedral == true
    std::vector<uint32_t> tets; // 4 per tet
    // 3 per vertex, each "[-]num[/den]"; empty unless exact coordinates were requested (-r).
    std::vector<std::string> exact_coords;
    // Per coplanar group of input triangles: the output tet faces tiling it, {tet, v0, v1, v2}.
    std::vector<std::vector<std::array<uint32_t, 4>>> triangle_provenance;
    // Per real input constraint: {input_triangle_id, coplanar_group_id}.
    std::vector<std::array<uint32_t, 2>> constraint_groups;
    uint32_t num_coplanar_groups = 0;
    // FaceProvenance::FACES: per tet, bit i set iff its face opposite vertex i tiles input
    // triangles (the faces triangle_provenance lists, marked on both their tets).
    std::vector<uint8_t> surface_faces;
    // Which of the above was computed. write_volume writes the .triangleprov and .groups
    // sidecars only for GROUPS (the tool's behaviour); they would otherwise be empty.
    FaceProvenance face_provenance = FaceProvenance::GROUPS;
    struct EdgeProvenance
    {
        uint32_t edge_id;
        std::vector<std::array<uint32_t, 3>> out_edges; // {tet, v0, v1}
    };
    std::vector<EdgeProvenance> edge_provenance;
    struct PointProvenance
    {
        uint32_t point_id;
        uint32_t tet; // UINT32_MAX if the point did not survive
        uint32_t out_vertex; // UINT32_MAX if the point did not survive
    };
    std::vector<PointProvenance> point_provenance;

    // tetrahedral == false: each cell as the list of its vertices.
    std::vector<uint32_t> cell_offsets; // num_cells + 1 entries into cell_vertices
    std::vector<uint32_t> cell_vertices;

    uint64_t num_vertices() const { return coords.size() / 3; }
    uint64_t num_tets() const { return tets.size() / 4; }
    uint64_t num_cells() const { return cell_offsets.empty() ? 0 : cell_offsets.size() - 1; }
};

// What one run of mesh_generator writes: exactly one of the three.
struct MeshOutput
{
    enum Kind { VOLUME, SKIN, BLACK_FACES };
    Kind kind = VOLUME;
    VolumeOutput volume;
    SurfaceOutput surface;
};

// Writers. Each reproduces the bytes the corresponding BSPcomplex::save* used to emit.
// The surface OFF writers keep the default stream precision (6 significant digits), as
// saveSkin/saveBlackFaces always did; the volume writers use 17.
void write_surface_off(const char* filename, const SurfaceOutput& out);
// volume.msh, or volume.tet with its .triangleprov/.groups/.edgeprov/.pointprov sidecars and,
// when exact_coords is filled, .rational.
void write_volume(const char* filename, const VolumeOutput& out);
void write_mesh_output(const char* filename, const MeshOutput& out);

} // namespace vol_rem

#endif // VOLUMEREMESHER_MESH_OUTPUT_H
