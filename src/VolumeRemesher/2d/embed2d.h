// Public wrapper for the 2D pipeline, mirroring embed_tri_in_poly_mesh in embed.h.
//
// Same contract as the 3D wrapper: flat input arrays of doubles, exact rational output
// coordinates, and provenance arrays indexed by INPUT id.

#ifndef VOLUMEREMESHER_2D_EMBED2D_H
#define VOLUMEREMESHER_2D_EMBED2D_H

#include <VolumeRemesher/numerics.h>

#include <array>
#include <cstdint>
#include <vector>

namespace vol_rem {

// Computes the arrangement of a soup of 2D segments and returns it as a triangulation.
//
// INPUT
//   seg_vrt_coords    2 doubles per point:  x0,y0, x1,y1, ...
//   segment_indexes   2 indices per segment, into seg_vrt_coords
//
//   Segments may intersect, overlap, be duplicated, or be degenerate. Duplicate points are
//   merged, zero-length segments are dropped, and both are still reported in the provenance
//   arrays (a dropped segment simply gets an empty edge list).
//
// OUTPUT
//   vertices          2 bigrationals per output vertex (x,y), EXACT
//   out_tris          output triangles, counterclockwise, indexing `vertices`
//   out_segment_provenance
//                     per INPUT segment, the output triangle edges tiling it, ordered from the
//                     segment's own first endpoint to its second. Each entry is
//                     {triangle, v0, v1} with v0 nearer the first endpoint.
//   out_point_provenance
//                     per INPUT point, {triangle, vertex}, or {UINT32_MAX, UINT32_MAX} if the
//                     point did not survive.
//
// The triangulation covers the input bounding box expanded by 10%; its four corners appear as
// extra output vertices. Returns false only if the input is empty or unusable.
bool embed_seg_in_tri_mesh(const std::vector<double>& seg_vrt_coords,
                           const std::vector<uint32_t>& segment_indexes,
                           std::vector<bigrational>& vertices,
                           std::vector<std::array<uint32_t, 3>>& out_tris,
                           std::vector<std::vector<std::array<uint32_t, 3>>>& out_segment_provenance,
                           std::vector<std::array<uint32_t, 2>>& out_point_provenance,
                           bool verbose);

} // namespace vol_rem

#endif
