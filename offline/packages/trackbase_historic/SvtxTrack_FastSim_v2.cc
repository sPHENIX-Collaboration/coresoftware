/*
 * SvtxTrack_FastSim_v2.C
 *
 *  Created on: Jul 28, 2016
 *      Author: yuhw
 */

#include "SvtxTrack_FastSim_v2.h"

#include "SvtxTrack.h"  // for SvtxTrack::ConstClusterIter, SvtxTrack

#include <map>      // for _Rb_tree_const_iterator
#include <ostream>  // for operator<<, basic_ostream, basic_ostream<>::_...

SvtxTrack_FastSim_v2::SvtxTrack_FastSim_v2(const SvtxTrack& source)
{ SvtxTrack_FastSim_v2::CopyFrom( source ); }

void SvtxTrack_FastSim_v2::CopyFrom( const SvtxTrack& source )
{
  // do nothing if copying onto oneself
  if( this == &source ) return;
 
  // parent class method
  SvtxTrack_v2::CopyFrom( source );

  // additional members
  _truth_track_id = source.get_truth_track_id();
  _nmeas = source.get_num_measurements();
  _g4hit_ids = source.g4hit_ids();
}

void SvtxTrack_FastSim_v2::identify(std::ostream& os) const
{
  SvtxTrack_v2::identify(os);

  os << "SvtxTrack_FastSim_v2 Object ";
  os << "_truth_Track_id: " << _truth_track_id << std::endl;
  os << "_nmeas: " << _nmeas << std::endl;

  os << "G4Hit IDs:" << std::endl;
  for( const auto& pair:_g4hit_ids )
  {
    os << "\thit container ID" << pair.first << " with hits: ";
    for( const auto& hitid:pair.second )
    { os << hitid << " "; }
    os << std::endl;
  }
  return;
}

int SvtxTrack_FastSim_v2::isValid() const
{
  return 1;
}
