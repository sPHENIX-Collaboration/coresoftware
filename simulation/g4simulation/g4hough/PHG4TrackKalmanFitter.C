////////////////////////////////////////////////////////////////////////////////
//
// This module is desgined to grab svtx tracks and put truth and cluster
// information into a TTree for GenFit testing
//
////////////////////////////////////////////////////////////////////////////////
//
// Darren McGlinchey
// 1 Apr 2016
//
////////////////////////////////////////////////////////////////////////////////


#include "PHG4TrackKalmanFitter.h"

#include <phool/phool.h>
#include <phool/getClass.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4VtxPoint.h>
#include <fun4all/PHTFileServer.h>
#include <fun4all/Fun4AllServer.h>

#include <g4hough/SvtxVertexMap.h>
#include <g4hough/SvtxVertex.h>
#include <g4hough/SvtxTrackMap.h>
#include <g4hough/SvtxTrack.h>
#include <g4hough/SvtxClusterMap.h>
#include <g4hough/SvtxCluster.h>
#include <g4hough/SvtxHitMap.h>
#include <g4hough/SvtxHit.h>

#include <g4eval/SvtxEvalStack.h>
#include <g4eval/SvtxTrackEval.h>
#include <g4eval/SvtxClusterEval.h>
#include <g4eval/SvtxTruthEval.h>
#include <g4eval/SvtxVertexEval.h>
#include <g4eval/SvtxHitEval.h>

#include <TTree.h>

#include <iostream>

using namespace std;

//----------------------------------------------------------------------------//
//-- Constructor:
//--  simple initialization
//----------------------------------------------------------------------------//
PHG4TrackKalmanFitter::PHG4TrackKalmanFitter(const string &name):
  SubsysReco( name ),
  _flags( NONE ),
  _tracks( NULL ),
  _svtxevalstack( NULL )
{
  //initialize
  _event = 0;
  _outfile = "PHG4TrackKalmanFitter.root";
}

//----------------------------------------------------------------------------//
//-- Init():
//--   Intialize all histograms, trees, and ntuples
//----------------------------------------------------------------------------//
int PHG4TrackKalmanFitter::Init(PHCompositeNode *topNode)
{
  cout << PHWHERE << " Openning file " << _outfile << endl;
  PHTFileServer::get().open( _outfile, "RECREATE");


  // create TTree
  _tracks = new TTree("tracks", "Svtx Tracks");
  _tracks->Branch("event", &event, "event/I");
  _tracks->Branch("gtrackID", &gtrackID, "gtrackID/I");
  _tracks->Branch("gflavor", &gflavor, "gflavor/I");
  _tracks->Branch("gpx", &gpx, "gpx/F");
  _tracks->Branch("gpy", &gpy, "gpy/F");
  _tracks->Branch("gpz", &gpz, "gpz/F");
  _tracks->Branch("gvx", &gvx, "gvx/F");
  _tracks->Branch("gvy", &gvy, "gvy/F");
  _tracks->Branch("gvz", &gvz, "gvz/F");
  _tracks->Branch("trackID", &trackID, "trackID/I");
  _tracks->Branch("charge", &charge, "charge/I");
  _tracks->Branch("nhits", &nhits, "nhits/I");
  _tracks->Branch("px", &px, "px/F");
  _tracks->Branch("py", &py, "py/F");
  _tracks->Branch("pz", &pz, "pz/F");
  _tracks->Branch("dca2d", &dca2d, "dca2d/F");
  _tracks->Branch("clusterID", &clusterID, "clusterID[nhits]/I");
  _tracks->Branch("layer", &layer, "layer[nhits]/I");
  _tracks->Branch("x", &x, "x[nhits]/F");
  _tracks->Branch("y", &y, "y[nhits]/F");
  _tracks->Branch("z", &z, "z[nhits]/F");
  _tracks->Branch("size_dphi", &size_dphi, "size_dphi[nhits]/F");
  _tracks->Branch("size_dz", &size_dz, "size_dz[nhits]/F");


  return 0;
}

//----------------------------------------------------------------------------//
//-- process_event():
//--   Call user instructions for every event.
//--   This function contains the analysis structure.
//----------------------------------------------------------------------------//
int PHG4TrackKalmanFitter::process_event(PHCompositeNode *topNode)
{
  _event++;
  if (_event % 1000 == 0)
    cout << PHWHERE << "Events processed: " << _event << endl;

  GetNodes(topNode);

  if (!_svtxevalstack) {
    _svtxevalstack = new SvtxEvalStack(topNode);
    _svtxevalstack->set_strict(false);
    _svtxevalstack->set_verbosity(verbosity + 1);
  } else {
    _svtxevalstack->next_event(topNode);
  }

  fill_tree(topNode);

  return 0;
}

//----------------------------------------------------------------------------//
//-- End():
//--   End method, wrap everything up
//----------------------------------------------------------------------------//
int PHG4TrackKalmanFitter::End(PHCompositeNode *topNode)
{

  PHTFileServer::get().cd( _outfile );

  _tracks->Write();

  if (_svtxevalstack) delete _svtxevalstack;

  return 0;
}


//----------------------------------------------------------------------------//
//-- fill_tree():
//--   Fill the trees with truth, track fit, and cluster information
//----------------------------------------------------------------------------//
void PHG4TrackKalmanFitter::fill_tree(PHCompositeNode *topNode)
{
  // Make sure to reset all the TTree variables before trying to set them.
  reset_variables();

  // get evaluators
  SvtxTrackEval*     trackeval = _svtxevalstack->get_track_eval();
  SvtxTruthEval*     trutheval = _svtxevalstack->get_truth_eval();

  if (_truth_container)
  {

    PHG4TruthInfoContainer::ConstRange range =
      _truth_container->GetPrimaryParticleRange();

    for (PHG4TruthInfoContainer::ConstIterator iter = range.first;
         iter != range.second;
         ++iter)
    {

      PHG4Particle* g4particle = iter->second;

      gtrackID = g4particle->get_track_id();
      gflavor  = g4particle->get_pid();

      PHG4VtxPoint* vtx = trutheval->get_vertex(g4particle);
      gvx      = vtx->get_x();
      gvy      = vtx->get_y();
      gvz      = vtx->get_z();

			gpx			 = g4particle->get_px();
			gpy			 = g4particle->get_py();
			gpz			 = g4particle->get_pz();

      SvtxTrack* track = trackeval->best_track_from(g4particle);
      if (track)
      {
        trackID   = track->get_id();
        charge    = track->get_charge();
        nhits     = track->size_clusters();
        px        = track->get_px();
        py        = track->get_py();
        pz        = track->get_pz();
        dca2d     = track->get_dca2d();


        int iclus = 0;
        for (SvtxTrack::ConstClusterIter iter = track->begin_clusters();
             iter != track->end_clusters();
             ++iter)
        {
          unsigned int cluster_id = *iter;
          SvtxCluster* cluster = _clustermap->get(cluster_id);
          unsigned int l = cluster->get_layer();

          clusterID[iclus] = (int)cluster_id;
          layer[iclus] = (int)l;
          x[iclus] = cluster->get_x();
          y[iclus] = cluster->get_y();
          z[iclus] = cluster->get_z();
          size_dphi[iclus] = cluster->get_phi_size();
          size_dz[iclus] = cluster->get_z_size();

          ++iclus;
        }

      } // if(track)
    } // for( iter)
  } //if (_truth_container)

  _tracks->Fill();
  return;

}

//----------------------------------------------------------------------------//
//-- reset_variables():
//--   Reset all the tree variables to their default values.
//--   Needs to be called at the start of every event
//----------------------------------------------------------------------------//
void PHG4TrackKalmanFitter::reset_variables()
{
  event = -9999;

  //-- truth
  gtrackID = -9999;
  gflavor = -9999;
  gpx = -9999;
  gpy = -9999;
  gpz = -9999;
  gvx = -9999;
  gvy = -9999;
  gvz = -9999;

  //-- reco
  trackID = -9999;
  charge = -9999;
  nhits = -9999;
  px = -9999;
  py = -9999;
  pz = -9999;
  dca2d = -9999;

  //-- clusters
  for (int i = 0; i < 7; i++)
  {
    clusterID[i] = -9999;
    layer[i] = -9999;
    x[i] = -9999;
    y[i] = -9999;
    z[i] = -9999;
    size_dphi[i] = -9999;
    size_dz[i] = -9999;
  }

}

//----------------------------------------------------------------------------//
//-- GetNodes():
//--   Get all the all the required nodes off the node tree
//----------------------------------------------------------------------------//
void PHG4TrackKalmanFitter::GetNodes(PHCompositeNode * topNode)
{
  //DST objects
  //Truth container
  _truth_container = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  if (!_truth_container && _event < 2)
  {
    cout << PHWHERE
         << " PHG4TruthInfoContainer node not found on node tree"
         << endl;
  }

  //Svtx Clusters
  _clustermap = findNode::getClass<SvtxClusterMap>(topNode, "SvtxClusterMap");
  if (!_clustermap && _event < 2)
  {
    cout << PHWHERE
         << " SvtxClusterMap node not found on node tree"
         << endl;
  }

}




