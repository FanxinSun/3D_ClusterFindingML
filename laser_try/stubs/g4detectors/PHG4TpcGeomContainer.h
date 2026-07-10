#pragma once
#include "PHG4TpcGeom.h"
class PHG4TpcGeomContainer
{ public: PHG4TpcGeom* GetLayerCellGeom(int) { static PHG4TpcGeom g; return &g; } };
