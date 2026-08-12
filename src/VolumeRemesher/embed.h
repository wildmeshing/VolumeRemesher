#include <array>
#include <vector>

#include "BSP.h"

namespace vol_rem {
/// <summary>
/// embed_tri_in_poly_mesh
/// Takes a set of triangles T and a tetrahedral mesh M.
/// Cuts tetrahedra in M so as to form a polyhedral mesh M' such that:
/// 1) each triangle in T is represented by the union of facets in M';
/// 2) each polyhedral cell in M' is (possibly weakly) convex.
/// Since cut points may be not representable using floating point
/// coordinates, M' is returned with rational coordinates.
/// </summary>
///
/// INPUT (T and M)
/// <param name="tri_vrt_coords">Vertex coordinates of T (x1,y1,z1,x2,y2,z2,...,xn,yn,zn)</param>
/// <param name="triangle_indexes">Triangle indexes of T (t1_v1, t1_v2, t1_v3, t2_v1, t2_v2, t2_v3, ..., tn_v1, tn_v2, tn_v3)</param>
/// <param name="tet_vrt_coords">Vertex coordinates of M (x1,y1,z1,x2,y2,z2,...,xn,yn,zn)</param>
/// <param name="tet_indexes">Tet indexes of M (t1_v1, t1_v2, t1_v3, t1_v4, t2_v1, t2_v2, t2_v3, t2_v4, ...)</param>
///
/// OUTPUT (M')
/// <param name="vertices">Vertex coordinates of M' (x1,y1,z1,x2,y2,z2,...,xn,yn,zn)</param>
/// <param name="facets">Polygonal facets in M'. A facet with n vertices is a sequence (n, p_v1, p_v2, ..., p_vn).
/// The first number in the sequence is the number of vertices in the polygon, whereas the subsequent n numbers are its
/// vertex indexes.</param>
/// <param name="cells">Polyhedral cells in M'. A polyhedron with n facets is a sequence (n, p_f1, p_f2, ..., p_fn).
/// The first number in the sequence is the number of facets bounding the polyhedron, whereas the subsequent n numbers are its
/// facet indexes.</param>
/// <param name="facets_on_input">Indexes of facets that overlap with T.</param>
///
/// <param name="verbose">Set to TRUE to enable verbosity.</param>
///
/// EXTRA 1D/0D FEATURES (optional; pass empty vectors to skip)
/// NOTE: as with the input surface, the extra edges/points -- and the forcing-triangle
/// apexes they generate, which are offset by the edge length (edges) or the average input
/// edge length (points) -- must lie inside the tet mesh M's domain.
/// <param name="edge_vrt_coords">Vertex coordinates of the extra edges (x,y,z,...)</param>
/// <param name="edge_indexes">Endpoint index pairs into edge_vrt_coords (e1_v1,e1_v2,...)</param>
/// <param name="point_coords">Coordinates of the extra points (x,y,z,...)</param>
///
/// PROVENANCE OUTPUTS. Symmetric: each entry is an output tet index plus the output
/// vertex indices of that tet's face / edge / vertex (indices into 'out_tets'/'vertices').
/// <param name="out_triangle_provenance">Per input coplanar group (see the surface tracking):
/// the output faces tiling it, each {tet, v0, v1, v2}.</param>
/// <param name="out_triangle_group">Per input triangle (same indexing as 'triangle_indexes'):
/// the coplanar group it belongs to, i.e. its index into 'out_triangle_provenance', or
/// UINT32_MAX if the triangle was degenerate and therefore dropped before the arrangement.
/// Triangles are grouped by transitive edge-adjacency AND exact coplanarity, so a flat region
/// tiled by many input triangles is one group and a triangle with no coplanar neighbour is its
/// own singleton; note this means a group can span triangles from several input surfaces where
/// they meet coplanarly along a shared edge. Without this the group ids in
/// 'out_triangle_provenance' cannot be mapped back to the input, which is what a caller needs
/// to carry per-triangle attributes (material tags, which input file a surface came from)
/// through the arrangement.</param>
/// <param name="out_edge_provenance">Per input edge: the output edges tiling it, each
/// {tet, v0, v1}.</param>
/// <param name="out_point_provenance">Per input point: {tet, vertex} of the output vertex equal
/// to it, or {UINT32_MAX, UINT32_MAX} if it did not survive into the output.</param>

void embed_tri_in_poly_mesh(
    const std::vector<double>& tri_vrt_coords,
    const std::vector<uint32_t>& triangle_indexes,
    const std::vector<double>& tet_vrt_coords,
    const std::vector<uint32_t>& tet_indexes,
    std::vector<bigrational>& vertices,
    std::vector<uint32_t>& facets,
    std::vector<uint32_t>& cells,
    std::vector<std::array<uint32_t, 4>>& out_tets,
    std::vector<uint32_t>& final_tets_parent,
    std::vector<uint32_t>& facets_on_input,
    std::vector<bool>& cells_with_faces_on_input,
    std::vector<std::vector<uint32_t>>& final_tets_parent_faces,
    const std::vector<double>& edge_vrt_coords,
    const std::vector<uint32_t>& edge_indexes,
    const std::vector<double>& point_coords,
    std::vector<std::vector<std::array<uint32_t, 4>>>& out_triangle_provenance,
    std::vector<uint32_t>& out_triangle_group,
    std::vector<std::vector<std::array<uint32_t, 3>>>& out_edge_provenance,
    std::vector<std::array<uint32_t, 2>>& out_point_provenance,
    bool verbose);

//
// TO DO: CONSIDER INPUT POINTS AND SEGMENTS
} // namespace vol_rem