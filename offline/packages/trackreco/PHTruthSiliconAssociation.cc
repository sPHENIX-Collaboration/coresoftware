//____________________________________________________________________________..
//
// This is a template for a Fun4All SubsysReco module with all methods from the
// $OFFLINE_MAIN/include/fun4all/SubsysReco.h baseclass
// You do not have to implement all of them, you can just remove unused methods
// here and in PHTruthSiliconAssociation.h.
//
// PHTruthSiliconAssociation(const std::string &name = "PHTruthSiliconAssociation")
// everything is keyed to PHTruthSiliconAssociation, duplicate names do work but it makes
// e.g. finding culprits in logs difficult or getting a pointer to the module
// from the command line
//
// PHTruthSiliconAssociation::~PHTruthSiliconAssociation()
// this is called when the Fun4AllServer is deleted at the end of running. Be
// mindful what you delete - you do loose ownership of object you put on the node tree
//
// int PHTruthSiliconAssociation::Init(PHCompositeNode *topNode)
// This method is called when the module is registered with the Fun4AllServer. You
// can create historgrams here or put objects on the node tree but be aware that
// modules which haven't been registered yet did not put antyhing on the node tree
//
// int PHTruthSiliconAssociation::InitRun(PHCompositeNode *topNode)
// This method is called when the first event is read (or generated). At
// this point the run number is known (which is mainly interesting for raw data
// processing). Also all objects are on the node tree in case your module's action
// depends on what else is around. Last chance to put nodes under the DST Node
// We mix events during readback if branches are added after the first event
//
// int PHTruthSiliconAssociation::process_event(PHCompositeNode *topNode)
// called for every event. Return codes trigger actions, you find them in
// $OFFLINE_MAIN/include/fun4all/Fun4AllReturnCodes.h
//   everything is good:
//     return Fun4AllReturnCodes::EVENT_OK
//   abort event reconstruction, clear everything and process next event:
//     return Fun4AllReturnCodes::ABORT_EVENT; 
//   proceed but do not save this event in output (needs output manager setting):
//     return Fun4AllReturnCodes::DISCARD_EVENT; 
//   abort processing:
//     return Fun4AllReturnCodes::ABORT_RUN
// all other integers will lead to an error and abort of processing
//
// int PHTruthSiliconAssociation::ResetEvent(PHCompositeNode *topNode)
// If you have internal data structures (arrays, stl containers) which needs clearing
// after each event, this is the place to do that. The nodes under the DST node are cleared
// by the framework
//
// int PHTruthSiliconAssociation::EndRun(const int runnumber)
// This method is called at the end of a run when an event from a new run is
// encountered. Useful when analyzing multiple runs (raw data). Also called at
// the end of processing (before the End() method)
//
// int PHTruthSiliconAssociation::End(PHCompositeNode *topNode)
// This is called at the end of processing. It needs to be called by the macro
// by Fun4AllServer::End(), so do not forget this in your macro
//
// int PHTruthSiliconAssociation::Reset(PHCompositeNode *topNode)
// not really used - it is called before the dtor is called
//
// void PHTruthSiliconAssociation::Print(const std::string &what) const
// Called from the command line - useful to print information when you need it
//
//____________________________________________________________________________..

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
  cout << "PHTruthSiliconAssociation::PHTruthSiliconAssociation(const std::string &name) Calling ctor" << endl;
}

//____________________________________________________________________________..
PHTruthSiliconAssociation::~PHTruthSiliconAssociation()
{
  cout << "PHTruthSiliconAssociation::~PHTruthSiliconAssociation() Calling dtor" << endl;
}

//____________________________________________________________________________..
int PHTruthSiliconAssociation::Init(PHCompositeNode *topNode)
{
  cout << "PHTruthSiliconAssociation::Init(PHCompositeNode *topNode) Initializing" << endl;

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int PHTruthSiliconAssociation::InitRun(PHCompositeNode *topNode)
{
  cout << "PHTruthSiliconAssociation::InitRun(PHCompositeNode *topNode) Initializing for Run XXX" << endl;

  GetNodes(topNode);

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int PHTruthSiliconAssociation::process_event(PHCompositeNode *topNode)
{
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
      //if (Verbosity() >= 1)
	{
	  std::cout
	    << __LINE__
	    << ": Processing seed itrack: " << phtrk_iter->first
	    << ": nhits: " << _tracklet-> size_cluster_keys()
	    << ": Total tracks: " << _track_map->size()
	    << endl;
	}

      _tracklet->identify(); 
     
      // Reject short seed tracks
      if(_tracklet->size_cluster_keys() < _min_clusters_per_track)
	{
	  std::cout << " --- reject track  number " << phtrk_iter->first << std::endl;
	  //continue;
	}

      // identify the best truth track match for this seed track 

	PHG4Particle *g4particle = trackeval->max_truth_particle_by_nclusters(_tracklet);

	// identify the clusters that are associated with this g4particle

	std::set<TrkrDefs::cluskey> clusters = clustereval->all_clusters_from(g4particle);

	for (std::set<TrkrDefs::cluskey>::iterator jter = clusters.begin();
	     jter != clusters.end();
	     ++jter)
	  {
	    TrkrDefs::cluskey cluster_key = *jter;
	    unsigned int layer = TrkrDefs::getLayer(cluster_key);
	    unsigned int trkrid = TrkrDefs::getTrkrId(cluster_key);
	    std::cout << "     found cluster with key: " << cluster_key << " in layer " << layer << std::endl;
	    //TrkrCluster *cluster = clustermap->findCluster(cluster_key);
	    //cluster->identify();
	    
	    // Identify the MVTX and INTT clusters and add them to the SvtxTrack cluster key list
	    if(trkrid == TrkrDefs::mvtxId)
	      {
		std::cout << "            cluster belongs to MVTX, add to track " << std::endl;
		_tracklet->insert_cluster_key(cluster_key);
		_assoc_container->SetClusterTrackAssoc(cluster_key, _tracklet->get_id());
	      }
	    
	    if(trkrid == TrkrDefs::inttId)
	      {
		std::cout << "            cluster belongs to INTT, add to track " << std::endl;
		  _tracklet->insert_cluster_key(cluster_key);
		_assoc_container->SetClusterTrackAssoc(cluster_key, _tracklet->get_id());
	      }
	  }

	// the vertex id from the seeder is nonsense, use vertex 0
	unsigned int ivert = 0;
	_tracklet->set_vertex_id(ivert);

	_tracklet->identify(); 
	std::cout << " new cluster keys size " << _tracklet->size_cluster_keys() << endl;

	std::cout << "Done with track " << phtrk_iter->first << std::endl;
    }
  
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
  //cout << "PHTruthSiliconAssociation::EndRun(const int runnumber) Ending Run for Run " << runnumber << endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int PHTruthSiliconAssociation::End(PHCompositeNode *topNode)
{
  //cout << "PHTruthSiliconAssociation::End(PHCompositeNode *topNode) This is the End..." << endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int PHTruthSiliconAssociation::Reset(PHCompositeNode *topNode)
{
  //cout << "PHTruthSiliconAssociation::Reset(PHCompositeNode *topNode) being Reset" << endl;
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

  //  _vertex_map = findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMap");
  //if (!_vertex_map)
  //{
  //cerr << PHWHERE << " ERROR: Can't find SvtxVertexMap." << endl;
  //return Fun4AllReturnCodes::ABORTEVENT;
  //}

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
