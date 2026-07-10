#pragma once
#include "PHCompositeNode.h"
#include <string>
template <class T> class PHIODataNode : public PHNode
{
 public:
  template <class U>
  PHIODataNode(U*, const std::string&, const std::string&) {}
};
