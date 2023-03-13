#ifndef TRACKBASEHISTORIC_SVTXTRACKSEED_V1_H
#define TRACKBASEHISTORIC_SVTXTRACKSEED_V1_H

#include <phool/PHObject.h>

#include "TrackSeed.h"

#include <cmath>
#include <iostream>

class SvtxTrackSeed_v1 : public TrackSeed
{
 public: 
  SvtxTrackSeed_v1();
  ~SvtxTrackSeed_v1() override;
  
  SvtxTrackSeed_v1( const SvtxTrackSeed_v1& );
  SvtxTrackSeed_v1& operator=(const SvtxTrackSeed_v1& seed);
 
  void identify(std::ostream& os = std::cout) const override;
  void Reset() override { *this = SvtxTrackSeed_v1(); }
  int isValid() const override { return 1; }
  void CopyFrom( const TrackSeed&) override;
  void CopyFrom( TrackSeed* seed) override { CopyFrom( *seed ); }
  PHObject* CloneMe() const override { return new SvtxTrackSeed_v1(*this); }

  unsigned int get_silicon_seed_index() const override { return m_silicon_seed; }
  unsigned int get_tpc_seed_index() const override { return m_tpc_seed; }
  void set_silicon_seed_index(const unsigned int index) override { m_silicon_seed = index; }
  void set_tpc_seed_index(const unsigned int index) override { m_tpc_seed = index; }

 private:

  unsigned int m_silicon_seed = std::numeric_limits<unsigned int>::max();
  unsigned int m_tpc_seed = std::numeric_limits<unsigned int>::max();
  
  ClassDefOverride(SvtxTrackSeed_v1, 1);

};

#endif 
