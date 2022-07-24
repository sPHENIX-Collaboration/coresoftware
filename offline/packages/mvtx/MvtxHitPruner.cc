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

#include <boost/graph/adjacency_list.hpp>
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

      cout << " found hitsetkey " << hitsetkey << " for bare_hitsetkey " << bare_hitsetkey << endl;
    }

  // Now consolidate all hits into the hitset with strobe 0, and delete the other hitsets
  //==============================================================
  for(auto bare_it = bare_hitset_set.begin(); bare_it != bare_hitset_set.end(); ++bare_it)
    {
      auto bare_hitsetkey = *bare_it;
      TrkrHitSet* bare_hitset = (m_hits->findOrAddHitSet(bare_hitsetkey))->second;
      std::cout << "         bare_hitset " << bare_hitsetkey << " initially has " << bare_hitset->size() << " hits " << std::endl; 

      auto bare_hitsetrange= hitset_multimap.equal_range(bare_hitsetkey);
      for(auto it = bare_hitsetrange.first; it != bare_hitsetrange.second; ++ it)
	{ 
	  auto hitsetkey = it->second;
	  cout << "            process hitsetkey " << hitsetkey << " for bare_hitsetkey " << bare_hitsetkey << endl;

	  int strobe = MvtxDefs::getStrobeId(hitsetkey);
	  if(strobe != 0)
	    {
	      // copy all hits to the hitset with strobe 0
	      TrkrHitSet* hitset = m_hits->findHitSet(hitsetkey);		

	      std::cout << "                hitsetkey " << hitsetkey << " has strobe " << strobe << " and has " << hitset->size() << " hits,  so copy it" << std::endl;

	      TrkrHitSet::ConstRange hitrangei = hitset->getHits();
	      for (TrkrHitSet::ConstIterator hitr = hitrangei.first;
		   hitr != hitrangei.second;
		   ++hitr)
		{
		  auto hitkey = hitr->first;
		  
		  // if it is already there, leave it alone, this is a duplicate hit
		  auto tmp_hit = bare_hitset->getHit(hitkey);
		  if(tmp_hit) continue;
		  
		  // otherwise copy the hit over 
		  std::cout << "                          copying over hitkey " << hitkey << std::endl;
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

  /*
  // Now loop over physical sensors one by one, make a vector of all hits from all hitsets in that sensor
  //===========================================================================
  
  for(auto bare_it = bare_hitset_set.begin(); bare_it != bare_hitset_set.end(); ++bare_it)
    {
      // bare hitset correponding to physical sensor
      // Make a vector of hits for the whole sensor to make things easier
      std::vector <std::pair< TrkrDefs::hitsetkey,  std::pair<TrkrDefs::hitkey, TrkrHit*> >> hitvec_initial;
      std::vector <std::pair<TrkrDefs::hitsetkey, std::pair< TrkrDefs::hitkey, TrkrHit*> >> hitvec;

      // get all hitsets (including all strobe bits) corresponding to this bare hitset
      auto bare_hitsetrange= hitset_multimap.equal_range(*bare_it);
      for(auto it = bare_hitsetrange.first; it != bare_hitsetrange.second; ++ it)
	{ 
	  auto hitsetkey = it->second;
	  auto bare_hitsetkey = it->first;

	  if(Verbosity() > 1) cout << " found hitsetkey " << hitsetkey << " for bare_hitsetkey " << bare_hitsetkey << endl;

	  TrkrHitSet *hitset = m_hits->findHitSet(hitsetkey);
	  if (Verbosity() > 2) hitset->identify();
	  
	  TrkrHitSet::ConstRange hitrangei = hitset->getHits();
	  for (TrkrHitSet::ConstIterator hitr = hitrangei.first;
	       hitr != hitrangei.second;
	       ++hitr)
	    {
	      // capture the hitkey and hit pointer for this hit
	      // All of this is for one sensor only
	      auto hitkey = hitr->first;
	      auto hit = hitr->second;
	      hitvec_initial.push_back(std::make_pair(hitsetkey, std::make_pair(hitkey, hit)));
	    }
	}
 
      //if (Verbosity() > 2)
      {
	TrkrDefs::hitsetkey hitsetkey = *bare_it;
	unsigned int layer = TrkrDefs::getLayer(hitsetkey);
	unsigned int stave = MvtxDefs::getStaveId(hitsetkey);
	unsigned int chip = MvtxDefs::getChipId(hitsetkey);
	unsigned int strobe = MvtxDefs::getStrobeId(hitsetkey);
	cout << "Bare hitsetkey " << hitsetkey << " layer " << layer << " stave " << stave << " chip " << chip << " strobe " << strobe << " hitvec_initial.size(): " << hitvec_initial.size() << endl;
      }  


      // All we need to know is the row and column of the hit in thhis sensor, we will assign it to a new hitset later
      // Capture row, column, hitkey, hit pointer
      //=================================================================
  
    // this multimap relates row and column of the physical sensor to the hit
      std::multimap<std::pair<unsigned int, unsigned int>, std::pair<TrkrDefs::hitsetkey, std::pair<TrkrDefs::hitkey, TrkrHit*>>> hitvec_map;
      // this set records the list of unique row and vector combinations found for the entire physical sensor
      std::set<std::pair<unsigned int, unsigned int>> hitvec_set;

      for(unsigned int i=0; i< hitvec_initial.size(); ++i)
	{
	  auto hitkey = hitvec_initial[i].second.first;
	  // these hits are all in the same chip, so we just detect repeat hits with the same col and row
	  auto col = MvtxDefs::getCol(hitkey);
	  auto row = MvtxDefs::getRow(hitkey);
	  auto row_col = std::make_pair(row, col);
	  hitvec_set.insert(row_col);
	  hitvec_map.insert(std::make_pair(row_col, hitvec_initial[i]));
	  std::cout << "    initial hits adding row " << row << " col " << col << std::endl;
	}

      // Choose only one version of each hit here and put it in hitvec
      //===============================================

      // go through hitvec_map and take just one hit per row and col value
      for(auto it = hitvec_set.begin(); it != hitvec_set.end(); ++it)
	{
	  auto ret = hitvec_map.equal_range(*it);
	  for(auto hitit = ret.first; hitit  != ret.second; ++hitit)
	    {
	      hitvec.push_back(hitit->second);

	      std::cout << "Add hit with row " << hitit->first.first << " col " << hitit->first.second << " hitsetkey " << hitit->second.first << " hitkey " << hitit->second.second.first << std::endl;
	      break;
	    }
	}
 
      //if (Verbosity() > 2) 
      cout << "hitvec.size() after pruning: " << hitvec.size() << endl;

      // copy the retained hits for this sensor to the hitset with strobe 0
      //================================================

      auto bare_hitsetkey =  *bare_it; 
      auto bare_hitset = (m_hits->findOrAddHitSet(bare_hitsetkey))->second;

      for(unsigned int i=0; i< hitvec.size(); ++i)
	{
	  auto hitkey = hitvec[i].second.first;

	  // if it is already there, leave it alone
	  auto tmp_hit = bare_hitset->getHit(hitkey);
	  if(tmp_hit) continue;

	  // otherwise copy the hit over 
	  TrkrHit* old_hit = hitvec[i].second.second;
	  TrkrHit *new_hit = new TrkrHitv2();
	  new_hit->setAdc(old_hit->getAdc());
	  bare_hitset->addHitSpecificKey(hitkey, new_hit);
	}
      
    }  // bare hitsets (i.e. physical sensor) loop
  

  // All chosen hits have been copied to the hitset with strobe 0 
  // now reset all hitsets with strobe != 0
  //============================
  std::vector<TrkrDefs::hitsetkey> remove_hitsets;  
  TrkrHitSetContainer::ConstRange mvtx_hitsets =
    m_hits->getHitSets(TrkrDefs::TrkrId::mvtxId);
  for (TrkrHitSetContainer::ConstIterator hitsetitr = mvtx_hitsets.first;
       hitsetitr != mvtx_hitsets.second;
       ++hitsetitr)
    {
      auto hitsetkey = hitsetitr->first;

      // get the strobe ID
      unsigned int strobe = MvtxDefs::getStrobeId(hitsetkey);
      if(strobe != 0)
	{
	  remove_hitsets.push_back(hitsetkey);
	}
    }

  for(unsigned int i = 0; i < remove_hitsets.size(); ++ i)
    {
      m_hits->removeHitSet(remove_hitsets[i]);
    }
  */

  return Fun4AllReturnCodes::EVENT_OK;
}

