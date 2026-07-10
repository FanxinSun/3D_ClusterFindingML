#pragma once
#include <string>
#include "../g4main/PHG4HitContainer.h"
class PHCompositeNode;
extern PHG4HitContainer* g_stub_container;
class PHG4TpcGeomContainer;
class PHG4TruthInfoContainer;
class SvtxTrackMap;
namespace findNode
{
template <class T>
T* getClass(PHCompositeNode*, const std::string&) { return nullptr; }
template <>
inline PHG4HitContainer* getClass<PHG4HitContainer>(PHCompositeNode*, const std::string&) { return g_stub_container; }
}
#include "../g4detectors/PHG4TpcGeomContainer.h"
#include "../g4main/PHG4TruthInfoContainer.h"
#include "../trackbase_historic/SvtxTrackMap_v2.h"
namespace findNode
{
template <>
inline PHG4TpcGeomContainer* getClass<PHG4TpcGeomContainer>(PHCompositeNode*, const std::string&)
{ static PHG4TpcGeomContainer g; return &g; }
template <>
inline PHG4TruthInfoContainer* getClass<PHG4TruthInfoContainer>(PHCompositeNode*, const std::string&)
{ static PHG4TruthInfoContainer t; return &t; }
template <>
inline SvtxTrackMap* getClass<SvtxTrackMap>(PHCompositeNode*, const std::string&)
{ static SvtxTrackMap_v2 m; return &m; }
}
