#pragma once
#include "PHCompositeNode.h"
class PHNodeIterator
{
 public:
  explicit PHNodeIterator(PHCompositeNode*) {}
  PHCompositeNode* findFirst(const char*, const char*) { static PHCompositeNode n; return &n; }
};
