#ifndef TRACKBASEHISTORIC_SVTXTRACKSEED_V3_H
#define TRACKBASEHISTORIC_SVTXTRACKSEED_V3_H

#include "TrackSeed.h"

#include <phool/PHObject.h>

#include <cmath>
#include <iostream>
#include <limits>

class SvtxTrackSeed_v3 : public TrackSeed
{
 public:
  SvtxTrackSeed_v3();
  ~SvtxTrackSeed_v3() override;

  SvtxTrackSeed_v3(const SvtxTrackSeed_v3&);
  SvtxTrackSeed_v3& operator=(const SvtxTrackSeed_v3& seed);

  void identify(std::ostream& os = std::cout) const override;
  void Reset() override { *this = SvtxTrackSeed_v3(); }
  int isValid() const override { return 1; }
  void CopyFrom(const TrackSeed&) override;
  void CopyFrom(TrackSeed* seed) override { CopyFrom(*seed); }
  PHObject* CloneMe() const override { return new SvtxTrackSeed_v3(*this); }

  unsigned int get_silicon_seed_index() const override { return m_silicon_seed; }
  unsigned int get_tpc_seed_index() const override { return m_tpc_seed; }
  short int get_crossing_estimate() const override { return m_crossing_estimate; }
  short int get_crossing() const override { return m_crossing; }
  void set_silicon_seed_index(const unsigned int index) override { m_silicon_seed = index; }
  void set_tpc_seed_index(const unsigned int index) override { m_tpc_seed = index; }
  void set_crossing_estimate(const short int cross) override { m_crossing_estimate = cross; }
  void set_crossing(const short int cross) override { m_crossing = cross; }

 private:
  unsigned int m_silicon_seed {std::numeric_limits<unsigned int>::max()};
  unsigned int m_tpc_seed {std::numeric_limits<unsigned int>::max()};
  short int m_crossing_estimate {std::numeric_limits<short>::max()};
  short int m_crossing {std::numeric_limits<short>::max()};

  ClassDefOverride(SvtxTrackSeed_v3, 1);
};

#endif
