#pragma once
#include <map>
#include <string>
class PHCompositeNode;
class PHParameterInterface
{
 public:
  explicit PHParameterInterface(const std::string&) {}
  virtual ~PHParameterInterface() = default;
  virtual void SetDefaultParameters() = 0;
  void InitializeParameters() { SetDefaultParameters(); }
  void set_default_double_param(const std::string& k, double v) { if (!dpar.count(k)) dpar[k] = v; }
  void set_default_int_param(const std::string& k, int v) { if (!ipar.count(k)) ipar[k] = v; }
  void set_default_string_param(const std::string& k, const std::string& v) { if (!spar.count(k)) spar[k] = v; }
  void set_double_param(const std::string& k, double v) { dpar[k] = v; }
  void set_int_param(const std::string& k, int v) { ipar[k] = v; }
  double get_double_param(const std::string& k) const { return dpar.at(k); }
  int get_int_param(const std::string& k) const { return ipar.at(k); }
  void UpdateParametersWithMacro() {}
  void SaveToNodeTree(PHCompositeNode*, const std::string&) {}
  void PutOnParNode(PHCompositeNode*, const std::string&) {}
 protected:
  std::map<std::string, double> dpar;
  std::map<std::string, int> ipar;
  std::map<std::string, std::string> spar;
};
