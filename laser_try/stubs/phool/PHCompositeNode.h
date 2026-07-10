#pragma once
class PHNode { public: virtual ~PHNode() = default; };
#include <string>
class PHCompositeNode : public PHNode
{
 public:
  PHCompositeNode() = default;
  explicit PHCompositeNode(const std::string&) {}
  void addNode(PHNode*) {}
};
