#pragma once

#include <implicit_point.h>
#include <vector>
#include <iomanip>
#include <cstring>
#include <assert.h>
#include <iostream>
#include <fstream>

namespace vol_rem {

#include "delaunay.h"

using namespace Del3D;
using namespace IPs;
using namespace NFG;

typedef basicVec3d pointType;
typedef TetMesh_t<pointType> TetMesh;
} // namespace vol_rem
