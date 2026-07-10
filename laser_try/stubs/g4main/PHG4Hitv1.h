#pragma once
#include "PHG4Hit.h"
class PHG4Hitv1 : public PHG4Hit
{
 public:
  PHG4Hitv1() = default;
  explicit PHG4Hitv1(const PHG4Hit* s) { *(PHG4Hit*)this = *s; }
};
