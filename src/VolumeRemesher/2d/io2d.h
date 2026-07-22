// I/O for the 2D arrangement pipeline.
//
// Input is an OBJ line soup ("v x y z" + "l i j", 1-based, z ignored). Output mirrors the naming
// and the sidecar conventions of the 3D pipeline's saveMesh (BSP.cpp): an approximate OFF plus an
// exact ".rational" companion and provenance sidecars.

#ifndef VOLUMEREMESHER_2D_IO2D_H
#define VOLUMEREMESHER_2D_IO2D_H

#include "arrangement2d.h"

#include <cstdint>
#include <string>
#include <vector>

namespace vol_rem {
namespace vr2d {

// Reads "v"/"l" records. Returns false if the file cannot be opened or is not a line soup.
// Coordinates are returned as (x,y) pairs; any third coordinate is dropped.
bool read_OBJ_segments(const std::string& path, std::vector<double>& coords,
                       std::vector<uint32_t>& indexes);

// Writes the triangulation.
//
//   <base>            OFF, z = 0, approximate double coordinates
//   <base>.rational   exact rational coordinates, same vertex order (only if export_rational)
//   <base>.segmentprov   per input segment: <id> <n> then n triples "<tri> <v0> <v1>"
//   <base>.pointprov     per input point:   <id> <tri> <vertex>
//
// Only the finite triangles are emitted, and only the vertices they use; the four bounding-box
// corners and the virtual triangles are dropped.
bool write_arrangement(const std::string& base, const Arrangement2D& A, bool export_rational);

} // namespace vr2d
} // namespace vol_rem

#endif
