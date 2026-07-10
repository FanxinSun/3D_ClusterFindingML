#pragma once
#include <string>
class PHCompositeNode;
class SubsysReco
{
 public:
  explicit SubsysReco(const std::string& n = "") : m_name(n) {}
  virtual ~SubsysReco() = default;
  virtual int InitRun(PHCompositeNode*) { return 0; }
  virtual int process_event(PHCompositeNode*) { return 0; }
  int Verbosity() const { return m_verb; }
  void Verbosity(int v) { m_verb = v; }
  std::string Name() const { return m_name; }
 private:
  std::string m_name;
  int m_verb = 0;
};
