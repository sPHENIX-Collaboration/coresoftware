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
#include <trackbase_historic/SvtxVertex.h>         // for SvtxVertex
#include <trackbase_historic/SvtxVertexMap.h>

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4VtxPoint.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/PHTFileServer.h>
#include <fun4all/SubsysReco.h>  // for SubsysReco

#include <phool/getClass.h>
#include <phool/phool.h>

#include <TH2.h>
#include <TSystem.h>
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
  , m_TruthInfoContainer(nullptr)
  , m_TrackMap(nullptr)
  , m_VertexMap(nullptr)
  , m_TracksEvalTree(nullptr)
  , m_VertexEvalTree(nullptr)
  , m_H2D_DeltaMomVsTruthMom(nullptr)
  , m_H2D_DeltaMomVsTruthEta(nullptr)
  , m_EventCounter(0)
  , _outfile_name(filename)
  , _trackmapname(trackmapname)
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
  m_TracksEvalTree = new TTree("tracks", "FastSim Eval => tracks");
  m_TracksEvalTree->Branch("event", &m_TTree_Event, "event/I");
  m_TracksEvalTree->Branch("gtrackID", &gtrackID, "gtrackID/I");
  m_TracksEvalTree->Branch("gflavor", &gflavor, "gflavor/I");
  m_TracksEvalTree->Branch("gpx", &gpx, "gpx/F");
  m_TracksEvalTree->Branch("gpy", &gpy, "gpy/F");
  m_TracksEvalTree->Branch("gpz", &gpz, "gpz/F");
  m_TracksEvalTree->Branch("gvx", &gvx, "gvx/F");
  m_TracksEvalTree->Branch("gvy", &gvy, "gvy/F");
  m_TracksEvalTree->Branch("gvz", &gvz, "gvz/F");
  m_TracksEvalTree->Branch("gvt", &gvt, "gvt/F");
  m_TracksEvalTree->Branch("trackID", &trackID, "trackID/I");
  m_TracksEvalTree->Branch("charge", &charge, "charge/I");
  m_TracksEvalTree->Branch("nhits", &nhits, "nhits/I");
  m_TracksEvalTree->Branch("px", &px, "px/F");
  m_TracksEvalTree->Branch("py", &py, "py/F");
  m_TracksEvalTree->Branch("pz", &pz, "pz/F");
  m_TracksEvalTree->Branch("pcax", &pcax, "pcax/F");
  m_TracksEvalTree->Branch("pcay", &pcay, "pcay/F");
  m_TracksEvalTree->Branch("pcaz", &pcaz, "pcaz/F");
  m_TracksEvalTree->Branch("dca2d", &dca2d, "dca2d/F");
  const string xyz[3] = {"x", "y", "z"};
  for (map<string, int>::const_iterator iter = m_ProjectionNameMap.begin(); iter != m_ProjectionNameMap.end(); ++iter)
  {
    char bname[100];
    char bdef[100];
    for (int i = 0; i < 3; i++)
    {
      sprintf(bname, "%s_%s", iter->first.c_str(), xyz[i].c_str());
      sprintf(bdef, "%s/F", bname);
      m_TracksEvalTree->Branch(bname, &ref[i][iter->second], bdef);
    }
    for (int i = 0; i < 3; i++)
    {
      sprintf(bname, "%s_p%s", iter->first.c_str(), xyz[i].c_str());
      sprintf(bdef, "%s/F", bname);
      m_TracksEvalTree->Branch(bname, &ref_p[i][iter->second], bdef);
    }
    for (int i = 0; i < 3; i++)
    {
      sprintf(bname, "%s_proj_%s", iter->first.c_str(), xyz[i].c_str());
      sprintf(bdef, "%s/F", bname);
      m_TracksEvalTree->Branch(bname, &proj[i][iter->second], bdef);
    }
    for (int i = 0; i < 3; i++)
    {
      sprintf(bname, "%s_proj_p%s", iter->first.c_str(), xyz[i].c_str());
      sprintf(bdef, "%s/F", bname);
      m_TracksEvalTree->Branch(bname, &proj_p[i][iter->second], bdef);
    }
  }

  m_H2D_DeltaMomVsTruthEta = new TH2D("DeltaMomVsTruthEta",
                                         "#frac{#Delta p}{truth p} vs. truth #eta", 54, -4.5, +4.5, 1000, -1,
                                         1);

  m_H2D_DeltaMomVsTruthMom = new TH2D("DeltaMomVsTruthMom",
                                         "#frac{#Delta p}{truth p} vs. truth p", 41, -0.5, 40.5, 1000, -1,
                                         1);

  // create TTree - vertex
  m_VertexEvalTree = new TTree("vertex", "FastSim Eval => vertces");
  m_VertexEvalTree->Branch("event", &m_TTree_Event, "event/I");
  m_VertexEvalTree->Branch("gvx", &gvx, "gvx/F");
  m_VertexEvalTree->Branch("gvy", &gvy, "gvy/F");
  m_VertexEvalTree->Branch("gvz", &gvz, "gvz/F");
  m_VertexEvalTree->Branch("gvt", &gvt, "gvt/F");
  m_VertexEvalTree->Branch("vx", &vx, "vx/F");
  m_VertexEvalTree->Branch("vy", &vy, "vy/F");
  m_VertexEvalTree->Branch("vz", &vz, "vz/F");
  m_VertexEvalTree->Branch("deltavx", &deltavx, "deltavx/F");
  m_VertexEvalTree->Branch("deltavy", &deltavy, "deltavy/F");
  m_VertexEvalTree->Branch("deltavz", &deltavz, "deltavz/F");
  m_VertexEvalTree->Branch("gID", &gtrackID, "gID/I");
  m_VertexEvalTree->Branch("ID", &trackID, "ID/I");
  m_VertexEvalTree->Branch("ntracks", &ntracks, "ntracks/I");
  m_VertexEvalTree->Branch("n_from_truth", &n_from_truth, "n_from_truth/I");

  return Fun4AllReturnCodes::EVENT_OK;
}

//----------------------------------------------------------------------------//
//-- process_event():
//--   Call user instructions for every event.
//--   This function contains the analysis structure.
//----------------------------------------------------------------------------//
int PHG4TrackFastSimEval::process_event(PHCompositeNode *topNode)
{
  m_EventCounter++;
  if (Verbosity() >= 2 and m_EventCounter % 1000 == 0)
    cout << PHWHERE << "Events processed: " << m_EventCounter << endl;

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

  m_TracksEvalTree->Write();
  m_VertexEvalTree->Write();

  m_H2D_DeltaMomVsTruthEta->Write();
  m_H2D_DeltaMomVsTruthMom->Write();

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

  if (!m_TruthInfoContainer)
  {
    LogError("m_TruthInfoContainer not found!");
    return;
  }

  if (!m_TrackMap)
  {
    LogError("m_TrackMap not found!");
    return;
  }

  PHG4TruthInfoContainer::ConstRange range =
      m_TruthInfoContainer->GetPrimaryParticleRange();
  //std::cout << "A2" << std::endl;
  for (PHG4TruthInfoContainer::ConstIterator truth_itr = range.first;
       truth_itr != range.second; ++truth_itr)
  {
    reset_variables();
    //std::cout << "A1" << std::endl;
    m_TTree_Event = m_EventCounter;

    PHG4Particle *g4particle = truth_itr->second;
    if (!g4particle)
    {
      LogDebug("");
      continue;
    }
    //std::cout << "B1" << std::endl;

    SvtxTrack_FastSim *track = nullptr;

    //std::cout << "TRACKmap size " << m_TrackMap->size() << std::endl;
    for (SvtxTrackMap::ConstIter track_itr = m_TrackMap->begin();
         track_itr != m_TrackMap->end();
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
    PHG4VtxPoint *vtx = m_TruthInfoContainer->GetVtx(g4particle->get_vtx_id());
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

      m_H2D_DeltaMomVsTruthMom->Fill(truth_mom.Mag(), (reco_mom.Mag() - truth_mom.Mag()) / truth_mom.Mag());
      m_H2D_DeltaMomVsTruthEta->Fill(truth_mom.Eta(), (reco_mom.Mag() - truth_mom.Mag()) / truth_mom.Mag());
      // find projections
      for (int k = 0; k < 3; k++)
      {
        for (int j = 0; j < nproj; j++)
        {
          proj[k][j] = -9999;
          proj_p[k][j] = -9999;
          ref[k][j] = -9999;
          ref_p[k][j] = -9999;
        }
      }
      for (SvtxTrack::ConstStateIter trkstates = track->begin_states();
           trkstates != track->end_states();
           ++trkstates)
      {
        //	cout << "checking " << trkstates->second->get_name() << endl;
        map<string, int>::const_iterator iter = m_ProjectionNameMap.find(trkstates->second->get_name());
        if (iter != m_ProjectionNameMap.end())
        {
          //	  cout << "found " << trkstates->second->get_name() << endl;
          // setting the projection (xyz and pxpypz)
          proj[0][iter->second] = trkstates->second->get_x();
          proj[1][iter->second] = trkstates->second->get_y();
          proj[2][iter->second] = trkstates->second->get_z();
          proj_p[0][iter->second] = trkstates->second->get_px();
          proj_p[1][iter->second] = trkstates->second->get_py();
          proj_p[2][iter->second] = trkstates->second->get_pz();

          string nodename = "G4HIT_" + trkstates->second->get_name();
          PHG4HitContainer *hits = findNode::getClass<PHG4HitContainer>(topNode, nodename);
          if (!hits)
          {
            cout << "could not find " << nodename << endl;
            continue;
          }
          //	  cout << "number of hits: " << hits->size() << endl;
          PHG4HitContainer::ConstRange hit_range = hits->getHits();
          for (PHG4HitContainer::ConstIterator hit_iter = hit_range.first; hit_iter != hit_range.second; hit_iter++)
          {
            //	    cout << "checking hit id " << hit_iter->second->get_trkid() << " against " << track->get_truth_track_id() << endl;
            if (hit_iter->second->get_trkid() - track->get_truth_track_id() == 0)
            {
              //	      cout << "found hit with id " << hit_iter->second->get_trkid() << endl;
              if (iter->second >= nproj)
              {
                cout << "bad index: " << iter->second << endl;
                gSystem->Exit(1);
              }
              ref[0][iter->second] = hit_iter->second->get_x(0);
              ref[1][iter->second] = hit_iter->second->get_y(0);
              ref[2][iter->second] = hit_iter->second->get_z(0);
              ref_p[0][iter->second] = hit_iter->second->get_px(0);
              ref_p[1][iter->second] = hit_iter->second->get_py(0);
              ref_p[2][iter->second] = hit_iter->second->get_pz(0);
            }
          }
        }
      }
    }
    //std::cout << "B3" << std::endl;

    m_TracksEvalTree->Fill();
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
  if (!m_TruthInfoContainer)
  {
    LogError("m_TruthInfoContainer not found!");
    return;
  }

  if (!m_TrackMap)
  {
    LogError("m_TrackMap not found!");
    return;
  }

  if (!m_VertexMap)
  {
    return;
  }

  for (SvtxVertexMap::Iter iter = m_VertexMap->begin();
       iter != m_VertexMap->end();
       ++iter)
  {
    SvtxVertex *vertex = iter->second;

    // Make sure to reset all the TTree variables before trying to set them.
    reset_variables();
    //std::cout << "A1" << std::endl;
    m_TTree_Event = m_EventCounter;

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
      const auto trackIter = m_TrackMap->find(trackID);

      if (trackIter == m_TrackMap->end()) continue;

      SvtxTrack_FastSim *temp = dynamic_cast<SvtxTrack_FastSim *>(trackIter->second);

      if (!temp) continue;

      const auto g4trackID = temp->get_truth_track_id();
      const PHG4Particle *g4particle = m_TruthInfoContainer->GetParticle(g4trackID);
      assert(g4particle);
      PHG4VtxPoint *vtx = m_TruthInfoContainer->GetVtx(g4particle->get_vtx_id());

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
    m_VertexEvalTree->Fill();
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
  m_TTree_Event = -9999;

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
  m_TruthInfoContainer = findNode::getClass<PHG4TruthInfoContainer>(topNode,
                                                                "G4TruthInfo");
  if (!m_TruthInfoContainer && m_EventCounter < 2)
  {
    cout << PHWHERE << " PHG4TruthInfoContainer node not found on node tree"
         << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_TrackMap = findNode::getClass<SvtxTrackMap>(topNode,
                                               _trackmapname);
  //std::cout << _trackmapname << std::endl;
  if (!m_TrackMap)
  {
    cout << PHWHERE << "SvtxTrackMap node with name "
         << _trackmapname
         << " not found on node tree"
         << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_VertexMap = findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMap");
  if (!m_VertexMap && Verbosity())
  {
    cout << PHWHERE << "SvtxTrackMap node with name SvtxVertexMap not found on node tree. Will not build the vertex eval tree"
         << endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

void PHG4TrackFastSimEval::AddProjection(const string &name)
{
  unsigned int size =  m_ProjectionNameMap.size();
  if (size >= nproj)
  {
    cout << "Too many projections in evaluator, maximum number is "
	 << nproj << endl;
    cout << "increase nproj in PHG4TrackFastSimEval and recompile" << endl;
    gSystem->Exit(1);
  }
  m_ProjectionNameMap.insert(make_pair(name, m_ProjectionNameMap.size()));
  return;
}
