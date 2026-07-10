#pragma once
#include "SvtxTrack_v2.h"
class SvtxTrackMap
{ public: virtual ~SvtxTrackMap() = default; void insert(SvtxTrack_v2*) {} unsigned size() const { return 0; } };
