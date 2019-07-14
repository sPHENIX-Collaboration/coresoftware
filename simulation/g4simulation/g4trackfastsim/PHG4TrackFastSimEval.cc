/*!
 *  \file		PHG4TrackFastSimEval.cc
 *  \brief		Evaluation module for PHG4TrackFastSim output
 *  \details	input: PHG4TruthInfoContainer, SvtxTrackMap with SvtxTrack_FastSim inside
 *  \author		Haiwang Yu <yuhw@nmsu.edu>
 */

#include "PHG4TrackFastSimEval.h"

#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrack_FastSim.h>
#include <trackbase_historic/SvtxVertexMap.h>

#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4VtxPoint.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/PHTFileServer.h>
#include <fun4all/SubsysReco.h>  // for SubsysReco

#include <phool/getClass.h>
#include <phool/phool.h>

#include <TH2.h>
#include <TTree.h>
#include <TVector3.h>

#include <cassert>
#include <cmath>
#include <iostream>
#include <map>      // for _Rb_tree_const_ite...
#include <utility>  // for pair

#define LogDebug(exp) std::cout << "DEBUG: " << __FILE__ << ": " << __LINE__ << ": " << exp << "\n"
#define LogError(exp) std::cout << "ERROR: " << __FILE__ << ": " << __LINE__ << ": " << exp << "\n"
#define LogWarning(exp) std::cout << "WARNING: " << __FILE__ << ": " << __LINE__ << ": " << exp << "\n"

using namespace std;

//----------------------------------------------------------------------------//
//-- Constructor:
//--  simple initialization
//----------------------------------------------------------------------------//
PHG4TrackFastSimEval::PHG4TrackFastSimEval(const string &name, const string &filename, const string &trackmapname)
  : SubsysReco(name)
  , _outfile_name(filename)
  , _trackmapname(trackmapname)
  , _event(0)
  , _flags(NONE)
  , _eval_tree_tracks(nullptr)
  , _eval_tree_vertex(nullptr)
  , _h2d_Delta_mom_vs_truth_mom(nullptr)
  , _h2d_Delta_mom_vs_truth_eta(nullptr)
  , _truth_container(nullptr)
  , _trackmap(nullptr)
  , _vertexmap(nullptr)
{
  reset_variables();
}

//----------------------------------------------------------------------------//
//-- Init():
//--   Intialize all histograms, trees, and ntuples
//----------------------------------------------------------------------------//
int PHG4TrackFastSimEval::Init(PHCompositeNode *topNode)
{
  cout << PHWHERE << " Openning file " << _outfile_name << endl;
  PHTFileServer::get().open(_outfile_name, "RECREATE");

  // create TTree
  _eval_tree_tracks = new TTree("tracks", "FastSim Eval => tracks");
  _eval_tree_tracks->Branch("event", &event, "event/I");
  _eval_tree_tracks->Branch("gtrackID", &gtrackID, "gtrackID/I");
  _eval_tree_tracks->Branch("gflavor", &gflavor, "gflavor/I");
  _eval_tree_tracks->Branch("gpx", &gpx, "gpx/F");
  _eval_tree_tracks->Branch("gpy", &gpy, "gpy/F");
  _eval_tree_tracks->Branch("gpz", &gpz, "gpz/F");
  _eval_tree_tracks->Branch("gvx", &gvx, "gvx/F");
  _eval_tree_tracks->Branch("gvy", &gvy, "gvy/F");
  _eval_tree_tracks->Branch("gvz", &gvz, "gvz/F");
  _eval_tree_tracks->Branch("gvt", &gvt, "gvt/F");
  _eval_tree_tracks->Branch("trackID", &trackID, "trackID/I");
  _eval_tree_tracks->Branch("charge", &charge, "charge/I");
  _eval_tree_tracks->Branch("nhits", &nhits, "nhits/I");
  _eval_tree_tracks->Branch("px", &px, "px/F");
  _eval_tree_tracks->Branch("py", &py, "py/F");
  _eval_tree_tracks->Branch("pz", &pz, "pz/F");
  _eval_tree_tracks->Branch("pcax", &pcax, "pcax/F");
  _eval_tree_tracks->Branch("pcay", &pcay, "pcay/F");
  _eval_tree_tracks->Branch("pcaz", &pcaz, "pcaz/F");
  _eval_tree_tracks->Branch("dca2d", &dca2d, "dca2d/F");

  _h2d_Delta_mom_vs_truth_eta = new TH2D("_h2d_Delta_mom_vs_truth_eta",
                                         "#frac{#Delta p}{truth p} vs. truth #eta", 54, -4.5, +4.5, 1000, -1,
                                         1);

  _h2d_Delta_mom_vs_truth_mom = new TH2D("_h2d_Delta_mom_vs_truth_mom",
                                         "#frac{#Delta p}{truth p} vs. truth p", 41, -0.5, 40.5, 1000, -1,
                                         1);

  // create TTree - vertex
  _eval_tree_vertex = new TTree("vertex", "FastSim Eval => vertces");
  _eval_tree_vertex->Branch("gvx", &gvx, "gvx/F");
  _eval_tree_vertex->Branch("gvy", &gvy, "gvy/F");
  _eval_tree_vertex->Branch("gvz", &gvz, "gvz/F");
  _eval_tree_vertex->Branch("gvt", &gvt, "gvt/F");
  _eval_tree_vertex->Branch("vx", &vx, "vx/F");
  _eval_tree_vertex->Branch("vy", &vy, "vy/F");
  _eval_tree_vertex->Branch("vz", &vz, "vz/F");
  _eval_tree_vertex->Branch("deltavx", &deltavx, "deltavx/F");
  _eval_tree_vertex->Branch("deltavy", &deltavy, "deltavy/F");
  _eval_tree_vertex->Branch("deltavz", &deltavz, "deltavz/F");
  _eval_tree_vertex->Branch("gID", &gtrackID, "gID/I");
  _eval_tree_vertex->Branch("ID", &trackID, "ID/I");
  _eval_tree_vertex->Branch("ntracks", &ntracks, "ntracks/I");
  _eval_tree_vertex->Branch("n_from_truth", &n_from_truth, "n_from_truth/I");

  return Fun4AllReturnCodes::EVENT_OK;
}

//----------------------------------------------------------------------------//
//-- process_event():
//--   Call user instructions for every event.
//--   This function contains the analysis structure.
//----------------------------------------------------------------------------//
int PHG4TrackFastSimEval::process_event(PHCompositeNode *topNode)
{
  _event++;
  if (Verbosity() >= 2 and _event % 1000 == 0)
    cout << PHWHERE << "Events processed: " << _event << endl;

  //std::cout << "Opening nodes" << std::endl;
  GetNodes(topNode);

  //std::cout << "Filling trees" << std::endl;
  fill_track_tree(topNode);
  fill_vertex_tree(topNode);
  //std::cout << "DONE" << std::endl;

  return Fun4AllReturnCodes::EVENT_OK;
}

//----------------------------------------------------------------------------//
//-- End():
//--   End method, wrap everything up
//----------------------------------------------------------------------------//
int PHG4TrackFastSimEval::End(PHCompositeNode *topNode)
{
  PHTFileServer::get().cd(_outfile_name);

  _eval_tree_tracks->Write();
  _eval_tree_vertex->Write();

  _h2d_Delta_mom_vs_truth_eta->Write();
  _h2d_Delta_mom_vs_truth_mom->Write();

  //PHTFileServer::get().close();

  return Fun4AllReturnCodes::EVENT_OK;
}

//----------------------------------------------------------------------------//
//-- fill_tree():
//--   Fill the trees with truth, track fit, and cluster information
//----------------------------------------------------------------------------//
void PHG4TrackFastSimEval::fill_track_tree(PHCompositeNode *topNode)
{
  // Make sure to reset all the TTree variables before trying to set them.

  if (!_truth_container)
  {
    LogError("_truth_container not found!");
    return;
  }

  if (!_trackmap)
  {
    LogError("_trackmap not found!");
    return;
  }

  PHG4TruthInfoContainer::ConstRange range =
      _truth_container->GetPrimaryParticleRange();
  //std::cout << "A2" << std::endl;
  for (PHG4TruthInfoContainer::ConstIterator truth_itr = range.first;
       truth_itr != range.second; ++truth_itr)
  {
    reset_variables();
    //std::cout << "A1" << std::endl;
    event = _event;

    PHG4Particle *g4particle = truth_itr->second;
    if (!g4particle)
    {
      LogDebug("");
      continue;
    }
    //std::cout << "B1" << std::endl;

    SvtxTrack_FastSim *track = nullptr;

    //std::cout << "TRACKmap size " << _trackmap->size() << std::endl;
    for (SvtxTrackMap::ConstIter track_itr = _trackmap->begin();
         track_itr != _trackmap->end();
         track_itr++)
    {
      //std::cout << "TRACK * " << track_itr->first << std::endl;
      SvtxTrack_FastSim *temp = dynamic_cast<SvtxTrack_FastSim *>(track_itr->second);
      if (!temp)
      {
        std::cout << "ERROR CASTING PARTICLE!" << std::endl;
        continue;
      }
      //std::cout << " PARTICLE!" << std::endl;

      if ((temp->get_truth_track_id() - g4particle->get_track_id()) == 0)
      {
        track = temp;
      }
    }

    //std::cout << "B2" << std::endl;
    gtrackID = g4particle->get_track_id();
    gflavor = g4particle->get_pid();

    gpx = g4particle->get_px();
    gpy = g4particle->get_py();
    gpz = g4particle->get_pz();

    gvx = NAN;
    gvy = NAN;
    gvz = NAN;
    gvt = NAN;
    PHG4VtxPoint *vtx = _truth_container->GetVtx(g4particle->get_vtx_id());
    if (vtx)
    {
      gvx = vtx->get_x();
      gvy = vtx->get_y();
      gvz = vtx->get_z();
      gvt = vtx->get_t();
    }

    if (track)
    {
      //std::cout << "C1" << std::endl;
      trackID = track->get_id();
      charge = track->get_charge();
      nhits = track->size_clusters();

      px = track->get_px();
      py = track->get_py();
      pz = track->get_pz();
      pcax = track->get_x();
      pcay = track->get_y();
      pcaz = track->get_z();
      dca2d = track->get_dca2d();

      TVector3 truth_mom(gpx, gpy, gpz);
      TVector3 reco_mom(px, py, pz);
      //std::cout << "C2" << std::endl;

      _h2d_Delta_mom_vs_truth_mom->Fill(truth_mom.Mag(), (reco_mom.Mag() - truth_mom.Mag()) / truth_mom.Mag());
      _h2d_Delta_mom_vs_truth_eta->Fill(truth_mom.Eta(), (reco_mom.Mag() - truth_mom.Mag()) / truth_mom.Mag());
    }
    //std::cout << "B3" << std::endl;

    _eval_tree_tracks->Fill();
  }
  //std::cout << "A3" << std::endl;

  return;
}

//----------------------------------------------------------------------------//
//-- fill_tree():
//--   Fill the trees with truth, track fit, and cluster information
//----------------------------------------------------------------------------//
void PHG4TrackFastSimEval::fill_vertex_tree(PHCompositeNode *topNode)
{
  if (!_truth_container)
  {
    LogError("_truth_container not found!");
    return;
  }

  if (!_trackmap)
  {
    LogError("_trackmap not found!");
    return;
  }

  if (!_vertexmap)
  {
    return;
  }

  for (SvtxVertexMap::Iter iter = _vertexmap->begin();
       iter != _vertexmap->end();
       ++iter)
  {
    SvtxVertex *vertex = iter->second;

    // Make sure to reset all the TTree variables before trying to set them.
    reset_variables();
    //std::cout << "A1" << std::endl;
    event = _event;

    if (!vertex)
    {
      LogDebug("");
      continue;
    }

    //std::cout << "C1" << std::endl;
    trackID = vertex->get_id();
    ntracks = vertex->size_tracks();

    vx = vertex->get_x();
    vy = vertex->get_y();
    vz = vertex->get_z();
    deltavx = sqrt(vertex->get_error(1, 1));
    deltavy = sqrt(vertex->get_error(2, 2));
    deltavz = sqrt(vertex->get_error(3, 3));

    // best matched vertex
    PHG4VtxPoint *best_vtx = nullptr;
    int best_n_match = -1;
    map<PHG4VtxPoint *, int> vertex_match_map;
    for (auto iter = vertex->begin_tracks(); iter != vertex->end_tracks(); ++iter)
    {
      const auto &trackID = *iter;
      const auto trackIter = _trackmap->find(trackID);

      if (trackIter == _trackmap->end()) continue;

      SvtxTrack_FastSim *temp = dynamic_cast<SvtxTrack_FastSim *>(trackIter->second);

      if (!temp) continue;

      const auto g4trackID = temp->get_truth_track_id();
      const PHG4Particle *g4particle = _truth_container->GetParticle(g4trackID);
      assert(g4particle);
      PHG4VtxPoint *vtx = _truth_container->GetVtx(g4particle->get_vtx_id());

      int n_match = ++vertex_match_map[vtx];

      if (n_match > best_n_match)
      {
        best_n_match = n_match;
        best_vtx = vtx;
      }
    }
    if (best_vtx)
    {
      gvx = best_vtx->get_x();
      gvy = best_vtx->get_y();
      gvz = best_vtx->get_z();
      gvt = best_vtx->get_t();

      n_from_truth = best_n_match;
      gtrackID = best_vtx->get_id();
    }
    _eval_tree_vertex->Fill();
  }
  //std::cout << "B3" << std::endl;

  return;
}

//----------------------------------------------------------------------------//
//-- reset_variables():
//--   Reset all the tree variables to their default values.
//--   Needs to be called at the start of every event
//----------------------------------------------------------------------------//
void PHG4TrackFastSimEval::reset_variables()
{
  event = -9999;

  //-- truth
  gtrackID = -9999;
  gflavor = -9999;
  gpx = NAN;
  gpy = NAN;
  gpz = NAN;

  gvx = NAN;
  gvy = NAN;
  gvz = NAN;
  gvt = NAN;

  //-- reco
  trackID = -9999;
  charge = -9999;
  nhits = -9999;
  px = NAN;
  py = NAN;
  pz = NAN;
  pcax = NAN;
  pcay = NAN;
  pcaz = NAN;
  dca2d = NAN;

  vx = NAN;
  vy = NAN;
  vz = NAN;
  deltavx = NAN;
  deltavy = NAN;
  deltavz = NAN;
  ntracks = -9999;
  n_from_truth = -9999;
}

//----------------------------------------------------------------------------//
//-- GetNodes():
//--   Get all the all the required nodes off the node tree
//----------------------------------------------------------------------------//
int PHG4TrackFastSimEval::GetNodes(PHCompositeNode *topNode)
{
  //DST objects
  //Truth container
  _truth_container = findNode::getClass<PHG4TruthInfoContainer>(topNode,
                                                                "G4TruthInfo");
  if (!_truth_container && _event < 2)
  {
    cout << PHWHERE << " PHG4TruthInfoContainer node not found on node tree"
         << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  _trackmap = findNode::getClass<SvtxTrackMap>(topNode,
                                               _trackmapname);
  //std::cout << _trackmapname.c_str() << std::endl;
  if (!_trackmap)
  {
    cout << PHWHERE << "SvtxTrackMap node with name "
         << _trackmapname
         << " not found on node tree"
         << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  _vertexmap = findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMap");
  if (!_vertexmap && _event < 2)
  {
    cout << PHWHERE << "SvtxTrackMap node with name SvtxVertexMap not found on node tree. Will not build the vertex eval tree"
         << endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}
