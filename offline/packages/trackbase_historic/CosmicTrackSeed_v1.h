#ifndef TRACKBASEHISTORIC_COSMICTRACKSEED_V1_H
#define TRACKBASEHISTORIC_COSMICTRACKSEED_V1_H

#include <phool/PHObject.h>

#include "TrackSeed.h"

#include <cmath>
#include <iostream>

class CosmicTrackSeed_v1 : public TrackSeed
{
 public:
  CosmicTrackSeed_v1();
  ~CosmicTrackSeed_v1() override;

  CosmicTrackSeed_v1(const CosmicTrackSeed_v1&);
  CosmicTrackSeed_v1& operator=(const CosmicTrackSeed_v1& seed);

  void identify(std::ostream& os = std::cout) const override;
  void Reset() override { *this = CosmicTrackSeed_v1(); }
  int isValid() const override { return 1; }
  void CopyFrom(const TrackSeed&) override;
  void CopyFrom(TrackSeed* seed) override { CopyFrom(*seed); }
  PHObject* CloneMe() const override { return new CosmicTrackSeed_v1(*this); }

  unsigned int get_silicon_seed_index() const override { return m_silicon_seed1; }
  unsigned int get_tpc_seed_index() const override { return m_tpc_seed1; }
  void set_silicon_seed_index(const unsigned int index) override { m_silicon_seed1 = index; }
  void set_tpc_seed_index(const unsigned int index) override { m_tpc_seed1 = index; }
  unsigned int get_silicon_seed_index2() const override { return m_silicon_seed2; }
  unsigned int get_tpc_seed_index2() const override { return m_tpc_seed2; }
  void set_silicon_seed_index2(const unsigned int index) override { m_silicon_seed2 = index; }
  void set_tpc_seed_index2(const unsigned int index) override { m_tpc_seed2 = index; }

 private:
  unsigned int m_silicon_seed1 = std::numeric_limits<unsigned int>::max();
  unsigned int m_tpc_seed1 = std::numeric_limits<unsigned int>::max();
  unsigned int m_silicon_seed2 = std::numeric_limits<unsigned int>::max();
  unsigned int m_tpc_seed2 = std::numeric_limits<unsigned int>::max();

  ClassDefOverride(CosmicTrackSeed_v1, 1);
};

#endif
