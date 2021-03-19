/*
 * SvtxTrack_FastSim_v1.C
 *
 *  Created on: Jul 28, 2016
 *      Author: yuhw
 */

#include "SvtxTrack_FastSim_v1.h"

#include "SvtxTrack.h"  // for SvtxTrack::ConstClusterIter, SvtxTrack

#include <climits>
#include <map>          // for _Rb_tree_const_iterator
#include <ostream>      // for operator<<, basic_ostream, basic_ostream<>::_...

using namespace std;

SvtxTrack_FastSim_v1::SvtxTrack_FastSim_v1()
{
}

SvtxTrack_FastSim_v1::~SvtxTrack_FastSim_v1()
{
}

void SvtxTrack_FastSim_v1::identify(std::ostream& os) const
{
  SvtxTrack_FastSim::identify(os);

  os << "G4Hit IDs" << endl;
  for (std::map<int, std::set<PHG4HitDefs::keytype> >::const_iterator iter =
           _g4hit_ids.begin();
       iter != _g4hit_ids.end(); ++iter)
  {
    for (std::set<PHG4HitDefs::keytype>::const_iterator jter =
             iter->second.begin();
         jter != iter->second.end(); ++jter)
    {
      os << *jter << " ";
    }
  }
  os << endl;
  return;
}

int SvtxTrack_FastSim_v1::isValid() const
{
  return 1;
}
