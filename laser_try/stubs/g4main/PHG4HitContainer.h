#pragma once
#include "PHG4Hit.h"
#include <vector>
class PHG4HitContainer
{
 public:
  std::vector<PHG4Hit*> hits;
  PHG4Hit* AddHit(int, PHG4Hit* h) { hits.push_back(h); return h; }
};
