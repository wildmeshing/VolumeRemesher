#ifndef _CONFORMING_MESH_
#define _CONFORMING_MESH_

#include "delaunay3d_wrapper.h"

#include <vector>

#define CONSTR_GROUP_T uint32_t
#define CONSTR_A 0
#define CONSTR_B 1

namespace vol_rem {
// half-edges struct
struct half_edge_t
{
    uint32_t endpts[2]; // endpoints of a constraint(triangle) edge,
    // be sure that endpts[0] < endpts[1].
    uint32_t tri_ind; // id_number of the constraint: tri_vertices[i]/3, where
    // tri_vertices[i]=endpts[0] or tri_vertices[i]=endpts[1].
};

// constraints are triangles
class constraints_t
{
public:
    uint32_t* tri_vertices;
    uint32_t num_triangles;
    uint32_t* constr_group;
    uint32_t num_virtual_triangles; // virtual constraints are indexed
    // from num_triangles - num_virtual_triangles
    // to num_triangles-1.
    uint32_t* tri_original_index; // per (real) constraint: its index in the input
    // file's triangle list, before degenerate triangles were dropped. NULL if not
    // tracked. Length = number of real (non-virtual) constraints (input + hole caps).
    uint32_t num_input_triangles; // number of genuine input surface constraints. Extra
    // constraints are appended after these in a fixed order: hole caps, then edge-forcing
    // triangles (0 to 2 per input edge), then point-forcing triangles (3 per input point, only
    // when the points are the whole input), then
    // the virtual (bounding-box) constraints. UINT32_MAX means "unset" (no caps/edges/points
    // added); callers then treat every real constraint as input surface. It can legitimately
    // be 0 (edges/points but no surface triangles), which is why the sentinel is not 0.
    uint32_t num_edge_triangles; // edge-forcing triangles (see insert_edges_and_points)
    uint32_t num_point_triangles; // point-forcing triangles
    // Per forced edge its two endpoints, per forced point its vertex, as mesh vertex ids
    // (filled by insert_edges_and_points). The provenance reads the features from here: a
    // feature can have no forcing triangle to read them from.
    std::vector<uint32_t> feature_edge_vertices;
    std::vector<uint32_t> feature_point_vertices;

    constraints_t()
        : tri_vertices(NULL)
        , num_triangles(0)
        , constr_group(NULL)
        , num_virtual_triangles(0)
        , tri_original_index(NULL)
        , num_input_triangles(UINT32_MAX)
        , num_edge_triangles(0)
        , num_point_triangles(0) {};
    ~constraints_t()
    {
        if (tri_vertices) free(tri_vertices);
        if (constr_group) free(constr_group);
        if (tri_original_index) free(tri_original_index);
    }
};


void fill_half_edges(const constraints_t* constraints, half_edge_t* half_edges);
void sort_half_edges(half_edge_t* half_edges, uint32_t num_half_edges);

// Detect the open boundaries (holes) of the current input constraints and cap
// each one with an ear-clipped triangulation, so that every input triangle bounds
// a closed volume. The cap ("fake") triangles are appended right after the real
// input triangles (num_triangles grows); constraints->num_input_triangles is set
// to the original count so the surface-tracking can exclude the caps. Operates on
// the deduplicated mesh vertices; must be called before the Delaunay permutation
// and before place_virtual_constraints. Returns the number of cap triangles added.
uint32_t fill_holes_in_constraints(constraints_t* constraints, const TetMesh* mesh, bool verbose);

// Extra 1D/0D features to force into the tetrahedralization output.
//   edges:  edge_verts holds n_edge_verts vertices (x,y,z...); edge_idx holds n_edges
//           endpoint index pairs into edge_verts.
//   points: point_verts holds n_points vertices (x,y,z...).
struct extra_features_t
{
    const double* edge_verts = nullptr;
    uint32_t n_edge_verts = 0;
    const uint32_t* edge_idx = nullptr;
    uint32_t n_edges = 0;
    const double* point_verts = nullptr;
    uint32_t n_points = 0;
};

// Force the given edges and points to appear in the output tetrahedralization. A point is
// appended to the mesh as a vertex, which makes it a vertex of the arrangement. An edge AB
// gets the triangles that pin it, sharing AB:
//  - none if it is an edge of the constraints already (the surface, hole caps) that the
//    arrangement has: an edge of exactly one constraint (it gets a virtual constraint) or of
//    two that are not coplanar (their planes cross along it);
//  - one if it is an edge of constraints that are all coplanar (inside a flat region), with
//    its third corner along the axis of their plane's largest normal component;
//  - otherwise two, the two largest-area triangles among (A, B, A + L e_{x,y,z}).
// The third corners are at a displacement L = 0.1 x |AB| from A. When the points are the
// whole input (no surface, no edge triangles) each gets the "corner" triangles
// (P, P + L e_i, P + L e_j), L = 0.1 x the point spacing, as the arrangement needs a
// constraint. New vertices are appended to the mesh and the forcing triangles to the
// constraints (num_edge_triangles / num_point_triangles record the counts, and
// feature_edge_vertices / feature_point_vertices the features). Must run before the
// Delaunay permutation and place_virtual_constraints.
//
// num_tet_mesh_vertices: when embedding into a given tet mesh (remakePolyhedralMesh), the
// number of its vertices -- the mesh's first vertices; those after them were appended for
// the surface and are inserted later. 0 when the mesh has no tets yet (makePolyhedralMesh,
// where the Delaunay tetrahedralization makes room for everything). In a given tet mesh
// every appended vertex has to lie strictly inside the mesh: a forced point or edge endpoint
// that is not an existing vertex and lies on or outside the boundary is an error, and a
// corner that would fall outside is moved -- to -L along its axis, then +-L/2, +-L/4, ...;
// when no axis direction fits (the base point is a boundary vertex of the mesh whose inward
// cone holds none, e.g. the apex of a rotated tet) the corners are taken inside a tet that
// holds the base point instead. Either way the triangles still share only AB (only P), so
// they pin it the same way; a flat edge's one corner is also kept off its plane. The mesh must
// be convex: see locate_point.
void insert_edges_and_points(
    constraints_t* constraints,
    TetMesh* mesh,
    const extra_features_t& extra,
    bool verbose,
    uint32_t num_tet_mesh_vertices = 0);

// Where point q lies with respect to a tet mesh (one with ghost tets around its boundary, as
// TetMesh::init_tets builds): PT_INSIDE in the open interior of a tet or on a face/edge
// interior to the mesh, PT_ON_HULL on a face, edge or vertex of the mesh's boundary, else
// PT_OUTSIDE. Exact. *walk_tet is the (corner-base) tet the walk starts from, updated to the
// tet found so that a nearby query starts close; *tet, if given, receives that tet for a
// point inside or on the hull. Walks by adjacency like TetMesh::searchTetrahedron, so the
// mesh must be convex: in a non-convex one the walk can leave through a reentrant face and
// report a point in another part of the mesh as PT_OUTSIDE.
enum point_location_t { PT_OUTSIDE, PT_ON_HULL, PT_INSIDE };
point_location_t
locate_point(const TetMesh* mesh, uint64_t* walk_tet, const pointType& q, uint64_t* tet);
uint32_t
place_virtual_constraints(TetMesh* mesh, constraints_t* constraints, half_edge_t* half_edges);
void insert_constraints(
    TetMesh*,
    constraints_t*,
    uint32_t*,
    uint32_t**,
    uint32_t*,
    uint32_t**,
    uint32_t*,
    uint32_t**,
    uint32_t*,
    uint32_t**,
    uint32_t*,
    uint32_t**);
} // namespace vol_rem

#endif
