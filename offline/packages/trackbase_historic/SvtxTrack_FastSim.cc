/*
 * SvtxTrack_FastSim.C
 *
 *  Created on: Jul 28, 2016
 *      Author: yuhw
 */

#include "SvtxTrack_FastSim.h"

#include "SvtxTrack.h"  // for SvtxTrack::ConstClusterIter, SvtxTrack

#include <climits>
#include <map>          // for _Rb_tree_const_iterator
#include <ostream>      // for operator<<, basic_ostream, basic_ostream<>::_...

SvtxTrack_FastSim::SvtxTrack_FastSim(const SvtxTrack& source)
{ SvtxTrack_FastSim::CopyFrom( source ); }

void SvtxTrack_FastSim::CopyFrom( const SvtxTrack& source )
{
  
  // parent class method
  SvtxTrack_v1::CopyFrom( source );
  
  // additional members
  _truth_track_id = source.get_truth_track_id();
  _nmeas = source.get_num_measurements();

}

void SvtxTrack_FastSim::identify(std::ostream& os) const
{
  SvtxTrack_v1::identify( os );

  os << "SvtxTrack_FastSim Object ";
  os << "truth_track_id:" << get_truth_track_id();
  os << "id: " << get_id() << " ";
  os << "charge: " << get_charge() << " ";
  os << "chisq: " << get_chisq() << " ndf:" << get_ndf() << " ";
  os << std::endl;

  os << "(px,py,pz) = ("
     << get_px() << ","
     << get_py() << ","
     << get_pz() << ")" << std::endl;

  os << "(x,y,z) = (" << get_x() << "," << get_y() << "," << get_z() << ")" << std::endl;

  if (!empty_clusters())
  {
    os << "clusters: ";
    for (SvtxTrack::ConstClusterIter iter = begin_clusters();
         iter != end_clusters();
         ++iter)
    {
      unsigned int cluster_id = *iter;
      os << cluster_id << " ";
    }
    os << std::endl;
  }

  return;
}

int SvtxTrack_FastSim::isValid() const
{
  return 1;
}
