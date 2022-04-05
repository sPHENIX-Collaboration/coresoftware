
#include "PHSiliconSeedMerger.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <phool/phool.h>
#include <phool/PHObject.h>
#include <phool/PHTimer.h>

#include <trackbase/TrkrDefs.h>

#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrack.h>


#include <phool/PHCompositeNode.h>

//____________________________________________________________________________..
PHSiliconSeedMerger::PHSiliconSeedMerger(const std::string &name):
 SubsysReco(name)
{}

//____________________________________________________________________________..
PHSiliconSeedMerger::~PHSiliconSeedMerger()
{
}

//____________________________________________________________________________..
int PHSiliconSeedMerger::Init(PHCompositeNode*)
{

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int PHSiliconSeedMerger::InitRun(PHCompositeNode *topNode)
{
  int ret = getNodes(topNode);
  return ret;
}

//____________________________________________________________________________..
int PHSiliconSeedMerger::process_event(PHCompositeNode *)
{

  std::multimap<unsigned int, std::set<TrkrDefs::cluskey>> matches;
  std::set<unsigned int> seedsToDelete;

  for(auto iter1 = m_siliconTrackMap->begin();
      iter1 != m_siliconTrackMap->end(); 
      ++iter1)
    {
      unsigned int track1ID = iter1->first;
      SvtxTrack* track1 = iter1->second;

      if(seedsToDelete.find(track1ID) != seedsToDelete.end())
	{ continue; }

      std::set<TrkrDefs::cluskey> mvtx1Keys;
      for (SvtxTrack::ConstClusterKeyIter iter = track1->begin_cluster_keys();
           iter != track1->end_cluster_keys();
           ++iter)
	{
	  TrkrDefs::cluskey ckey = *iter;
	  if(TrkrDefs::getTrkrId(ckey) == TrkrDefs::TrkrId::mvtxId)
	    { mvtx1Keys.insert(ckey); }
	}

      /// We can speed up the code by only iterating over the track seeds
      /// that are further in the map container from the current track,
      /// since the comparison of e.g. track 1 with track 2 doesn't need
      /// to be repeated with track 2 to track 1.
      for(auto iter2 = m_siliconTrackMap->find(track1ID);
	  iter2 != m_siliconTrackMap->end();
	  ++iter2) 
	{
	  unsigned int track2ID = iter2->first;
	  if(track1ID == track2ID) { continue; }
	  SvtxTrack* track2 = iter2->second;
	  std::set<TrkrDefs::cluskey> mvtx2Keys;
	  for (SvtxTrack::ConstClusterKeyIter iter = track2->begin_cluster_keys();
	       iter != track2->end_cluster_keys();
	       ++iter)
	    {
	      TrkrDefs::cluskey ckey = *iter;
	      if(TrkrDefs::getTrkrId(ckey) == TrkrDefs::TrkrId::mvtxId)
		{ mvtx2Keys.insert(ckey); }
	    }

	  std::vector<TrkrDefs::cluskey> intersection;
	  std::set_intersection(mvtx1Keys.begin(),
				mvtx1Keys.end(),
				mvtx2Keys.begin(),
				mvtx2Keys.end(),
				std::back_inserter(intersection));

	  /// If we have two clusters in common in the triplet, it is likely
	  /// from the same track
	  if(intersection.size() > 1) 
	    {
	      if(Verbosity() > 2) 
		{
		  std::cout << "Track " << track1ID << " keys " << std::endl;
		  for(auto& key : mvtx1Keys) 
		    { std::cout << "   ckey: " << key << std::endl; }
		  std::cout << "Track " << track2ID << " keys " << std::endl;
		  for(auto& key : mvtx2Keys) 
		    { std::cout << "   ckey: " << key << std::endl; }
		  std::cout << "Intersection keys " << std::endl;
		  for(auto& key : intersection)
		    { std::cout << "   ckey: " << key << std::endl; }
		}

	      for(auto& key : mvtx2Keys)
		{
		  if(mvtx1Keys.find(key) == mvtx1Keys.end())
		    { mvtx1Keys.insert(key); }
		}
	      
	      if(Verbosity() > 2)
		{ 
		  std::cout << "Match IDed"<<std::endl; 
		  for(auto& key : mvtx1Keys)
		    { std::cout << "  total track keys " << key << std::endl; }
		}

	      matches.insert(std::make_pair(track1ID, mvtx1Keys)); 
	      seedsToDelete.insert(track2ID);
	      break;
	    }
	}
    }

  for(const auto& [trackKey, mvtxKeys] : matches)
    {
      auto track = m_siliconTrackMap->get(trackKey);
      if(Verbosity() > 2)
	{ std::cout << "original track: " << std::endl; track->identify(); }

      for(auto& key : mvtxKeys) 
	{
	  if(track->find_cluster_key(key) == track->end_cluster_keys())
	    { 
	      track->insert_cluster_key(key); 
	      if(Verbosity() > 2) 
		std::cout << "adding " << key << std::endl;
	    }
	}
      
    }

  for(const auto& key : seedsToDelete) 
    {
      if(Verbosity() > 2 )
	{ std::cout << "Erasing track " << key << std::endl; }
      m_siliconTrackMap->erase(key);
    }

  if(Verbosity() > 2)
    {
      for(const auto& [key, track] : *m_siliconTrackMap)
	{ track->identify(); }
    }
	  
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int PHSiliconSeedMerger::ResetEvent(PHCompositeNode *)
{

  return Fun4AllReturnCodes::EVENT_OK;
}



//____________________________________________________________________________..
int PHSiliconSeedMerger::End(PHCompositeNode *)
{

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHSiliconSeedMerger::getNodes(PHCompositeNode *topNode)
{
  m_siliconTrackMap = findNode::getClass<SvtxTrackMap>(topNode, m_trackMapName.c_str());
  if(!m_siliconTrackMap)
    {
      std::cout << PHWHERE << "No silicon track map, can't merge seeds"
		<< std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  
return Fun4AllReturnCodes::EVENT_OK;
}

