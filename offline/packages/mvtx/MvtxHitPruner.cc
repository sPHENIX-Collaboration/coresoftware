/**
 * @file mvtx/MvtxHitPruner.cc
 * @author D. McGlinchey
 * @date June 2018
 * @brief Implementation of MvtxHitPruner
 */
#include "MvtxHitPruner.h"
#include "CylinderGeom_Mvtx.h"

#include <g4detectors/PHG4CylinderGeom.h>
#include <g4detectors/PHG4CylinderGeomContainer.h>

#include <trackbase/TrkrClusterContainerv4.h>
#include <trackbase/TrkrClusterv3.h>
#include <trackbase/TrkrClusterv4.h>
#include <trackbase/TrkrDefs.h>                     // for hitkey, getLayer
#include <trackbase/MvtxDefs.h>                   
#include <trackbase/TrkrHitv2.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>
#include <trackbase/TrkrClusterHitAssocv3.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>                     // for SubsysReco

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>                           // for PHNode
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>                         // for PHObject
#include <phool/getClass.h>
#include <phool/phool.h>                            // for PHWHERE

#include <TMatrixFfwd.h>                            // for TMatrixF
#include <TMatrixT.h>                               // for TMatrixT, operator*
#include <TMatrixTUtils.h>                          // for TMatrixTRow
#include <TVector3.h>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <boost/graph/adjacency_list.hpp>
#pragma GCC diagnostic pop

#include <boost/graph/connected_components.hpp>

#include <array>
#include <cmath>
#include <cstdlib>                                 // for exit
#include <iostream>
#include <map>                                      // for multimap<>::iterator
#include <set>                                      // for set, set<>::iterator
#include <string>
#include <vector>                                   // for vector

using namespace boost;
using namespace std;

MvtxHitPruner::MvtxHitPruner(const string &name)
  : SubsysReco(name)
  , m_hits(nullptr)
{
}

int MvtxHitPruner::InitRun(PHCompositeNode * /*topNode*/)
{
  //-----------------
  // Add Cluster Node
  //-----------------

  //----------------
  // Report Settings
  //----------------

  if (Verbosity() > 0)
  {
    cout << "====================== MvtxHitPruner::InitRun() =====================" << endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int MvtxHitPruner::process_event(PHCompositeNode *topNode)
{
  // get node containing the digitized hits
  m_hits = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  if (!m_hits)
  {
    cout << PHWHERE << "ERROR: Can't find node TRKR_HITSET" << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // We want to combine all strobe values for a given hitset
  // Start by looping over all MVTX hitsets and making a map of physical sensor to hitsetkey-with-strobe
  //=============================================================================
  std::multimap<TrkrDefs::hitsetkey, TrkrDefs::hitsetkey> hitset_multimap;  // will map (bare hitset, hitset with strobe)
  std::set<TrkrDefs::hitsetkey> bare_hitset_set;  // list of all physical sensor hitsetkeys (i.e. with strobe set to zero)

  TrkrHitSetContainer::ConstRange hitsetrange =
    m_hits->getHitSets(TrkrDefs::TrkrId::mvtxId);
  for (TrkrHitSetContainer::ConstIterator hitsetitr = hitsetrange.first;
       hitsetitr != hitsetrange.second;
       ++hitsetitr)
    {
      auto hitsetkey = hitsetitr->first;

      // get the hitsetkey value for strobe 0
      unsigned int layer = TrkrDefs::getLayer(hitsetitr->first);
      unsigned int stave =  MvtxDefs::getStaveId(hitsetitr->first);
      unsigned int chip =  MvtxDefs::getChipId(hitsetitr->first);
      auto bare_hitsetkey =  MvtxDefs::genHitSetKey(layer, stave, chip, 0);

      hitset_multimap.insert(std::make_pair(bare_hitsetkey, hitsetkey));
      bare_hitset_set.insert(bare_hitsetkey);

      if(Verbosity() > 0) cout << " found hitsetkey " << hitsetkey << " for bare_hitsetkey " << bare_hitsetkey << endl;
    }

  // Now consolidate all hits into the hitset with strobe 0, and delete the other hitsets
  //==============================================================
  for(auto bare_it = bare_hitset_set.begin(); bare_it != bare_hitset_set.end(); ++bare_it)
    {
      auto bare_hitsetkey = *bare_it;
      TrkrHitSet* bare_hitset = (m_hits->findOrAddHitSet(bare_hitsetkey))->second;
      if(Verbosity() > 0) std::cout << "         bare_hitset " << bare_hitsetkey << " initially has " << bare_hitset->size() << " hits " << std::endl; 

      auto bare_hitsetrange= hitset_multimap.equal_range(bare_hitsetkey);
      for(auto it = bare_hitsetrange.first; it != bare_hitsetrange.second; ++ it)
	{ 
	  auto hitsetkey = it->second;

	  int strobe = MvtxDefs::getStrobeId(hitsetkey);
	  if(strobe != 0)
	    {
	      if(Verbosity() > 0)  cout << "            process hitsetkey " << hitsetkey << " for bare_hitsetkey " << bare_hitsetkey << endl;

	      // copy all hits to the hitset with strobe 0
	      TrkrHitSet* hitset = m_hits->findHitSet(hitsetkey);		

	       if(Verbosity() > 0) 
		 std::cout << "                hitsetkey " << hitsetkey << " has strobe " << strobe << " and has " << hitset->size() << " hits,  so copy it" << std::endl;

	      TrkrHitSet::ConstRange hitrangei = hitset->getHits();
	      for (TrkrHitSet::ConstIterator hitr = hitrangei.first;
		   hitr != hitrangei.second;
		   ++hitr)
		{
		  auto hitkey = hitr->first;
		  if(Verbosity() > 0) std::cout << "                 found hitkey " << hitkey << std::endl;		  
		  // if it is already there, leave it alone, this is a duplicate hit
		  auto tmp_hit = bare_hitset->getHit(hitkey);
		  if(tmp_hit) 
		    {
		      if(Verbosity() > 0) std::cout << "                          hitkey " << hitkey << " is already in bare hitsest, do not copy" << std::endl;
		      continue;
		    }

		  // otherwise copy the hit over 
		   if(Verbosity() > 0)  std::cout << "                          copying over hitkey " << hitkey << std::endl;
		  auto old_hit = hitr->second;
		  TrkrHit *new_hit = new TrkrHitv2();
		  new_hit->setAdc(old_hit->getAdc());
		  bare_hitset->addHitSpecificKey(hitkey, new_hit);
		}

	      // all hits are copied over to the strobe zero hitset, remove this hitset
	      m_hits->removeHitSet(hitsetkey);
	    }
	}
    }

  return Fun4AllReturnCodes::EVENT_OK;
}

