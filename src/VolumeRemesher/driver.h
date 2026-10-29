// The pipeline behind mesh_generator, as a library call: the same options the command line
// takes, the same inputs its readers produce (io.h), and the output as an in-memory snapshot
// (mesh_output.h) that write_mesh_output turns into the same files. main.cpp is a thin argv
// parser over this, and the Python module binds it, so the two cannot drift apart.

#ifndef VOLUMEREMESHER_DRIVER_H
#define VOLUMEREMESHER_DRIVER_H

#include "2d/io2d.h"
#include "BSP.h" // CellClassifier
#include "io.h"
#include "mesh_output.h"

#include <cstdint>
#include <vector>

namespace vol_rem {

// One flag per command line option (the letter is in the comment).
struct MeshOptions
{
    bool verbose = false; // -v
    bool logging = false; // -l  append a line to mesh_generator.log
    bool triangulate = false; // -t  tetrahedralize the volume / triangulate the surface output
    bool export_rational = false; // -r  with -t, also carry the exact rational vertex coordinates
    bool fill_holes = false; // -f  cap open boundaries (solid-output path only)
    bool keep_all_cells = false; // -a  skip the in/out classification, keep the whole domain
    MeshOutput::Kind output = MeshOutput::VOLUME; // default / -s (SKIN) / -b (BLACK_FACES)
    // -c winding|graphcut  how cells are classified inside/outside. GRAPH_CUT needs a build
    // with VOLUMEREMESHER_WITH_GRAPHCUT, otherwise it is an error (ip_error).
    CellClassifier classifier = CellClassifier::WINDING_NUMBER;
    // With triangulate (no command line flags):
    FaceProvenance face_provenance = FaceProvenance::GROUPS; // how faces on the input are reported
    bool round_vertices = false; // coords: nearest doubles of the exact coordinates
};

// The file mesh_generator writes for these options: volume.msh, volume.tet, skin.off or
// black_faces.off.
const char* default_output_filename(const MeshOptions& opts);

// `mesh_generator [A.off] [op B.off] [-e edges] [-p points] <flags>`: the volume enclosed by
// surface A, or by the boolean of A and B (bool_opcode 'U', 'I' or 'D'; '0' for none), with
// optional extra edges and points forced into the output (they imply keep_all_cells). A may be
// null only with edges and/or points; B only with an opcode. name_A/name_B label the inputs in
// the log file and in verbose output. Errors go through ip_error.
MeshOutput make_mesh(
    TriMeshData* A,
    const char* name_A,
    char bool_opcode,
    TriMeshData* B,
    const char* name_B,
    const EdgeData* edges,
    const PointData* points,
    const MeshOptions& opts);

// `mesh_generator A.off X B.tet`: embed surface A (and optional extra edges/points) into the
// tetrahedral mesh B; the cells of B cut by A are classified against A and only the inside ones
// are written, as the tool always did. bool_opcode is what the command line passed as the
// separator; the writers apply it to the cell labels exactly as the tool always has.
MeshOutput embed_mesh(
    const TriMeshData& A,
    const TetMeshData& B,
    char bool_opcode,
    const EdgeData* edges,
    const PointData* points,
    const MeshOptions& opts);

// `mesh_generator -2d soup.obj [-r]`: the arrangement of a 2D segment soup as a triangulation.
// seg_coords holds 2 doubles per point, seg_indexes 2 indices per segment.
vr2d::ArrangementOutput arrange_2d(
    const std::vector<double>& seg_coords,
    const std::vector<uint32_t>& seg_indexes,
    bool export_rational,
    bool verbose);

} // namespace vol_rem

#endif // VOLUMEREMESHER_DRIVER_H
