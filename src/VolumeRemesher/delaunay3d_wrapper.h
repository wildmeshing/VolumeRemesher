#pragma once

#include <VolumeRemesher/implicit_point.h>
#include <vector>
#include <iomanip>
#include <cstring>
#include <assert.h>
#include <iostream>
#include <fstream>

namespace vol_rem {

#include "delaunay.h"

typedef basicVec3d pointType;
typedef TetMesh_t<pointType> TetMesh;
} // namespace vol_rem
