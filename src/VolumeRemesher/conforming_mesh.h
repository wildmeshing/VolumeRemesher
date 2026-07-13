#ifndef _CONFORMING_MESH_
#define _CONFORMING_MESH_

#include "delaunay.h"

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
    uint32_t num_input_triangles; // number of genuine input constraints. When holes
    // are filled, cap ("fake") constraints are appended after these and occupy
    // [num_input_triangles, num_triangles - num_virtual_triangles). 0 means "unset"
    // (no hole filling ran); callers then treat every real constraint as input.

    constraints_t()
        : tri_vertices(NULL)
        , num_triangles(0)
        , constr_group(NULL)
        , num_virtual_triangles(0)
        , tri_original_index(NULL)
        , num_input_triangles(0) {};
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
