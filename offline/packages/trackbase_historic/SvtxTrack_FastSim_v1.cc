/*
 * SvtxTrack_FastSim_v1.C
 *
 *  Created on: Jul 28, 2016
 *      Author: yuhw
 */

#include "SvtxTrack_FastSim_v1.h"

#include "SvtxTrack.h"  // for SvtxTrack::ConstClusterIter, SvtxTrack

#include <climits>
#include <map>      // for _Rb_tree_const_iterator
#include <ostream>  // for operator<<, basic_ostream, basic_ostream<>::_...

SvtxTrack_FastSim_v1::SvtxTrack_FastSim_v1(const SvtxTrack& source)
{ SvtxTrack_FastSim_v1::CopyFrom( source ); }

void SvtxTrack_FastSim_v1::CopyFrom( const SvtxTrack& source )
{
  // do nothing if copying onto oneself
  if( this == &source ) return;
 
  // parent class method
  SvtxTrack_FastSim::CopyFrom( source );
  
  // additional members
  _g4hit_ids = source.g4hit_ids();
}

void SvtxTrack_FastSim_v1::identify(std::ostream& os) const
{
  SvtxTrack_FastSim::identify(os);

  os << "SvtxTrack_FastSim_v1 Object ";
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

int SvtxTrack_FastSim_v1::isValid() const
{
  return 1;
}
