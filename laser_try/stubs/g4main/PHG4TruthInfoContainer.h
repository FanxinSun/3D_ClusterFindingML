#pragma once
#include "PHG4Particlev3.h"
#include "PHG4VtxPointv1.h"
class PHG4TruthInfoContainer
{
 public:
  int maxtrkindex() const { return 0; }
  int mintrkindex() const { return 0; }
  void AddParticle(int, PHG4Particlev3*) {}
  int maxvtxindex() const { return 0; }
  void AddVertex(int, PHG4VtxPointv1*) {}
};
