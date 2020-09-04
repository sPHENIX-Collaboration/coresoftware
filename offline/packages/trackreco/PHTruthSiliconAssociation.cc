#include "PHTruthSiliconAssociation.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <phool/phool.h>

/// Tracking includes
#include <trackbase/TrkrClusterv1.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxVertexMap.h>

#include <g4eval/SvtxClusterEval.h>
#include <g4eval/SvtxEvalStack.h>
#include <g4eval/SvtxEvaluator.h>
#include <g4eval/SvtxTrackEval.h>

#include <g4main/PHG4Hit.h>  // for PHG4Hit
#include <g4main/PHG4Particle.h>  // for PHG4Particle
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4HitDefs.h>  // for keytype
#include <g4main/PHG4TruthInfoContainer.h>

#include "AssocInfoContainer.h"

using namespace std;

//____________________________________________________________________________..
PHTruthSiliconAssociation::PHTruthSiliconAssociation(const std::string &name):
 SubsysReco(name)
{
  //cout << "PHTruthSiliconAssociation::PHTruthSiliconAssociation(const std::string &name) Calling ctor" << endl;
}

//____________________________________________________________________________..
PHTruthSiliconAssociation::~PHTruthSiliconAssociation()
{

}

//____________________________________________________________________________..
int PHTruthSiliconAssociation::Init(PHCompositeNode *topNode)
{

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int PHTruthSiliconAssociation::InitRun(PHCompositeNode *topNode)
{

  GetNodes(topNode);

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int PHTruthSiliconAssociation::process_event(PHCompositeNode *topNode)
{
  if (Verbosity() >= 1) 
    cout << "PHTruthSiliconAssociation::process_event(PHCompositeNode *topNode) Processing Event" << endl;

  _svtxEvalStack = new SvtxEvalStack(topNode);
  _svtxEvalStack->next_event(topNode);

  SvtxTrackEval *trackeval = _svtxEvalStack->get_track_eval();
  SvtxClusterEval *clustereval = _svtxEvalStack->get_cluster_eval();

  // Loop over all SvtxTracks from the CA seeder
  // These should contain all TPC clusters already

  for (auto phtrk_iter = _track_map->begin();
       phtrk_iter != _track_map->end(); 
       ++phtrk_iter)
    {
      _tracklet = phtrk_iter->second;
      if (Verbosity() >= 1)
	{
	  std::cout
	    << __LINE__
	    << ": Processing seed itrack: " << phtrk_iter->first
	    << ": nhits: " << _tracklet-> size_cluster_keys()
	    << ": Total tracks: " << _track_map->size()
	    << ": phi: " << _tracklet->get_phi()
	    << endl;
	}


      if (Verbosity() >= 1) 
	_tracklet->identify(); 

      // identify the best truth track match for this seed track 

	PHG4Particle *g4particle = trackeval->max_truth_particle_by_nclusters(_tracklet);

	if (Verbosity() >= 1)
	  std::cout << "   g4particleID " << g4particle->get_track_id() 
		    << " px " <<  g4particle->get_px() 
		    << " py " <<  g4particle->get_py() 
		    << " phi " << atan2(g4particle->get_py(), g4particle->get_px() )
		    << std::endl;

	// identify the clusters that are associated with this g4particle
	std::set<TrkrDefs::cluskey> clusters = clustereval->all_clusters_from(g4particle);

	for (std::set<TrkrDefs::cluskey>::iterator jter = clusters.begin();
	     jter != clusters.end();
	     ++jter)
	  {
	    TrkrDefs::cluskey cluster_key = *jter;
	    unsigned int layer = TrkrDefs::getLayer(cluster_key);
	    unsigned int trkrid = TrkrDefs::getTrkrId(cluster_key);
	     if (Verbosity() >= 1)  std::cout << "     found cluster with key: " << cluster_key << " in layer " << layer << std::endl;
	    
	    // Identify the MVTX and INTT clusters and add them to the SvtxTrack cluster key list
	    if(trkrid == TrkrDefs::mvtxId)
	      {
		 if (Verbosity() >= 1)  std::cout << "            cluster belongs to MVTX, add to track " << std::endl;
		_tracklet->insert_cluster_key(cluster_key);
		_assoc_container->SetClusterTrackAssoc(cluster_key, _tracklet->get_id());
	      }
	    
	    if(trkrid == TrkrDefs::inttId)
	      {
		 if (Verbosity() >= 1) std::cout << "            cluster belongs to INTT, add to track " << std::endl;
		  _tracklet->insert_cluster_key(cluster_key);
		_assoc_container->SetClusterTrackAssoc(cluster_key, _tracklet->get_id());
	      }
	  }

	// the vertex id from the seeder is nonsense, use vertex 0
	unsigned int ivert = 0;
	_tracklet->set_vertex_id(ivert);

	// set the track position to the vertex position
	const SvtxVertex *svtxVertex = _vertex_map->get(ivert);
	
	_tracklet->set_x(svtxVertex->get_x());
	_tracklet->set_y(svtxVertex->get_y());
	_tracklet->set_z(svtxVertex->get_z());
 
      if (Verbosity() >= 1)
	{
	  _tracklet->identify(); 
	  std::cout << " new cluster keys size " << _tracklet->size_cluster_keys() << endl;
	}

      if (Verbosity() >= 1)
	std::cout << "Done with track " << phtrk_iter->first << std::endl;
    }

  if (Verbosity() >= 1)
    cout << "PHTruthSiliconAssociation::process_event(PHCompositeNode *topNode) Leaving process_event" << endl;  

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int PHTruthSiliconAssociation::ResetEvent(PHCompositeNode *topNode)
{
  //cout << "PHTruthSiliconAssociation::ResetEvent(PHCompositeNode *topNode) Resetting internal structures, prepare for next event" << endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int PHTruthSiliconAssociation::EndRun(const int runnumber)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int PHTruthSiliconAssociation::End(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int PHTruthSiliconAssociation::Reset(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
void PHTruthSiliconAssociation::Print(const std::string &what) const
{
  //cout << "PHTruthSiliconAssociation::Print(const std::string &what) const Printing info for " << what << endl;
}

int  PHTruthSiliconAssociation::GetNodes(PHCompositeNode* topNode)
{
  //---------------------------------
  // Get Objects off of the Node Tree
  //---------------------------------

  _cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  if (!_cluster_map)
  {
    cerr << PHWHERE << " ERROR: Can't find node TRKR_CLUSTER" << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  _vertex_map = findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMap");
  if (!_vertex_map)
    {
      cerr << PHWHERE << " ERROR: Can't find SvtxVertexMap." << endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  
  _track_map = findNode::getClass<SvtxTrackMap>(topNode,  "SvtxTrackMap");
  if (!_track_map)
  {
    cerr << PHWHERE << " ERROR: Can't find SvtxTrackMap: " << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  _assoc_container = findNode::getClass<AssocInfoContainer>(topNode, "AssocInfoContainer");
  if (!_assoc_container)
  {
    cerr << PHWHERE << " ERROR: Can't find AssocInfoContainer." << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}
