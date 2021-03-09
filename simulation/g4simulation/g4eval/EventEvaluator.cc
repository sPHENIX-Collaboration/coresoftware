#include "EventEvaluator.h"

#include "CaloEvalStack.h"
#include "CaloRawClusterEval.h"
#include "CaloRawTowerEval.h"
#include "CaloTruthEval.h"

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4Shower.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4VtxPoint.h>

#include <g4vertex/GlobalVertex.h>
#include <g4vertex/GlobalVertexMap.h>

#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrack_FastSim.h>

#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawClusterUtility.h>
#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>

#include <phool/getClass.h>
#include <phool/phool.h>

#include <TFile.h>
#include <TNtuple.h>
#include <TTree.h>

#include <CLHEP/Vector/ThreeVector.h>

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <set>
#include <utility>

using namespace std;

EventEvaluator::EventEvaluator(const string& name, const string& filename)
  : SubsysReco(name)
  , _do_FHCAL(true)
  , _do_FEMC(true)
  , _do_DRCALO(true)
  , _do_HITS(true)
  , _do_TRACKS(true)
  , _do_CLUSTERS(true)
  , _do_VERTEX(true)
  , _do_PROJECTIONS(true)
  , _do_MCPARTICLES(true)
  , _ievent(0)
  , _nHitsLayers(0)
  , _hits_layerID(0)
  , _hits_x(0)
  , _hits_y(0)
  , _hits_z(0)
  , _hits_t(0)

  , _nTowers_FHCAL(0)
  , _tower_FHCAL_E(0)
  , _tower_FHCAL_iEta(0)
  , _tower_FHCAL_iPhi(0)
  , _tower_FHCAL_trueID(0)

  , _nTowers_DRCALO(0)
  , _tower_DRCALO_E(0)
  , _tower_DRCALO_iEta(0)
  , _tower_DRCALO_iPhi(0)
  , _tower_DRCALO_trueID(0)

  , _nTowers_FEMC(0)
  , _tower_FEMC_E(0)
  , _tower_FEMC_iEta(0)
  , _tower_FEMC_iPhi(0)
  , _tower_FEMC_trueID(0)

  , _nclusters_FHCAL(0)
  , _cluster_FHCAL_E(0)
  , _cluster_FHCAL_Eta(0)
  , _cluster_FHCAL_Phi(0)
  , _cluster_FHCAL_NTower(0)
  , _cluster_FHCAL_trueID(0)

  , _nclusters_FEMC(0)
  , _cluster_FEMC_E(0)
  , _cluster_FEMC_Eta(0)
  , _cluster_FEMC_Phi(0)
  , _cluster_FEMC_NTower(0)
  , _cluster_FEMC_trueID(0)

  , _vertex_x(0)
  , _vertex_y(0)
  , _vertex_z(0)

  , _nTracks(0)
  , _track_ID(0)
  , _track_px(0)
  , _track_py(0)
  , _track_pz(0)
  , _track_trueID(0)
  , _nProjections(0)
  , _track_ProjTrackID(0)
  , _track_ProjLayer(0)
  , _track_TLP_x(0)
  , _track_TLP_y(0)
  , _track_TLP_z(0)
  , _track_TLP_t(0)
  , _track_TLP_true_x(0)
  , _track_TLP_true_y(0)
  , _track_TLP_true_z(0)
  , _track_TLP_true_t(0)

  , _nMCPart(0)
  , _mcpart_ID(0)
  , _mcpart_ID_parent(0)
  , _mcpart_PDG(0)
  , _mcpart_E(0)
  , _mcpart_px(0)
  , _mcpart_py(0)
  , _mcpart_pz(0)

  , _reco_e_threshold(0.0)
  , _caloevalstackFHCAL(nullptr)
  , _caloevalstackDRCALO(nullptr)
  , _caloevalstackFEMC(nullptr)
  , _strict(false)
  , _event_tree(nullptr)
  , _filename(filename)
  , _tfile(nullptr)
{
  if(_do_HITS){
    _hits_layerID = new int[_maxNHits];
    _hits_x = new float[_maxNHits];
    _hits_y = new float[_maxNHits];
    _hits_z = new float[_maxNHits];
    _hits_t = new float[_maxNHits];
  }
  if(_do_FHCAL){
    _tower_FHCAL_E = new float[_maxNTowers];
    _tower_FHCAL_iEta  = new int[_maxNTowers];
    _tower_FHCAL_iPhi = new int[_maxNTowers];
    _tower_FHCAL_trueID = new int[_maxNTowers];
    if(_do_CLUSTERS){
      _cluster_FHCAL_E = new float[_maxNclusters];
      _cluster_FHCAL_Eta = new float[_maxNclusters];
      _cluster_FHCAL_Phi = new float[_maxNclusters];
      _cluster_FHCAL_NTower = new int[_maxNclusters];
      _cluster_FHCAL_trueID = new int[_maxNclusters];
    }
  }
  if(_do_DRCALO){
    _tower_DRCALO_E = new float[_maxNTowersDR];
    _tower_DRCALO_iEta  = new int[_maxNTowersDR];
    _tower_DRCALO_iPhi = new int[_maxNTowersDR];
    _tower_DRCALO_trueID = new int[_maxNTowersDR];
  }
  if(_do_FEMC){
    _tower_FEMC_E = new float[_maxNTowers];
    _tower_FEMC_iEta = new int[_maxNTowers];
    _tower_FEMC_iPhi = new int[_maxNTowers];
    _tower_FEMC_trueID = new int[_maxNTowers];
    if(_do_CLUSTERS){
      _cluster_FEMC_E = new float[_maxNclusters];
      _cluster_FEMC_Eta = new float[_maxNclusters];
      _cluster_FEMC_Phi = new float[_maxNclusters];
      _cluster_FEMC_NTower = new int[_maxNclusters];
      _cluster_FEMC_trueID = new int[_maxNclusters];
    }
  }
  if(_do_TRACKS){
    _track_ID = new float[_maxNTracks];
    _track_trueID = new float[_maxNTracks];
    _track_px = new float[_maxNTracks];
    _track_py= new float[_maxNTracks];
    _track_pz = new float[_maxNTracks];
  }
  if(_do_PROJECTIONS){
    _track_ProjTrackID = new float[_maxNProjections];
    _track_ProjLayer = new int[_maxNProjections];
    _track_TLP_x = new float[_maxNProjections];
    _track_TLP_y = new float[_maxNProjections];
    _track_TLP_z = new float[_maxNProjections];
    _track_TLP_t = new float[_maxNProjections];
    _track_TLP_true_x = new float[_maxNProjections];
    _track_TLP_true_y = new float[_maxNProjections];
    _track_TLP_true_z = new float[_maxNProjections];
    _track_TLP_true_t = new float[_maxNProjections];
  }
  if(_do_MCPARTICLES){
    _mcpart_ID = new float[_maxNMCPart];
    _mcpart_ID_parent = new float[_maxNMCPart];
    _mcpart_PDG = new float[_maxNMCPart];
    _mcpart_E = new float[_maxNMCPart];
    _mcpart_px = new float[_maxNMCPart];
    _mcpart_py = new float[_maxNMCPart];
    _mcpart_pz = new float[_maxNMCPart];
  }
}

int EventEvaluator::Init(PHCompositeNode* topNode)
{
  _ievent = 0;

  _tfile = new TFile(_filename.c_str(), "RECREATE");

  _event_tree = new TTree("event_tree", "event_tree");
  // tracks and hits
  if(_do_HITS){
    _event_tree->Branch("nHits", &_nHitsLayers, "nHits/I");
    _event_tree->Branch("hits_layerID", _hits_layerID, "hits_layerID[nHits]/I");
    _event_tree->Branch("hits_x", _hits_x, "hits_x[nHits]/F");
    _event_tree->Branch("hits_y", _hits_y, "hits_y[nHits]/F");
    _event_tree->Branch("hits_z", _hits_z, "hits_z[nHits]/F");
    _event_tree->Branch("hits_t", _hits_t, "hits_t[nHits]/F");
  }
  if(_do_TRACKS){
    _event_tree->Branch("nTracks", &_nTracks, "nTracks/I");
    _event_tree->Branch("tracks_ID", _track_ID, "tracks_ID[nTracks]/F");
    _event_tree->Branch("tracks_px", _track_px, "tracks_px[nTracks]/F");
    _event_tree->Branch("tracks_py", _track_py, "tracks_py[nTracks]/F");
    _event_tree->Branch("tracks_pz", _track_pz, "tracks_pz[nTracks]/F");
    _event_tree->Branch("tracks_trueID", _track_trueID, "tracks_trueID[nTracks]/F");
  }
  if(_do_PROJECTIONS){
    _event_tree->Branch("nProjections", &_nProjections, "nProjections/I");
    _event_tree->Branch("track_ProjTrackID", _track_ProjTrackID, "track_ProjTrackID[nProjections]/F");
    _event_tree->Branch("track_ProjLayer", _track_ProjLayer, "track_ProjLayer[nProjections]/I");
    _event_tree->Branch("track_TLP_x", _track_TLP_x, "track_TLP_x[nProjections]/F");
    _event_tree->Branch("track_TLP_y", _track_TLP_y, "track_TLP_y[nProjections]/F");
    _event_tree->Branch("track_TLP_z", _track_TLP_z, "track_TLP_z[nProjections]/F");
    _event_tree->Branch("track_TLP_t", _track_TLP_t, "track_TLP_t[nProjections]/F");
    _event_tree->Branch("track_TLP_true_x", _track_TLP_true_x, "track_TLP_true_x[nProjections]/F");
    _event_tree->Branch("track_TLP_true_y", _track_TLP_true_y, "track_TLP_true_y[nProjections]/F");
    _event_tree->Branch("track_TLP_true_z", _track_TLP_true_z, "track_TLP_true_z[nProjections]/F");
    _event_tree->Branch("track_TLP_true_t", _track_TLP_true_t, "track_TLP_true_t[nProjections]/F");
  }
  if(_do_FHCAL){
    // towers HCAL
    _event_tree->Branch("tower_FHCAL_N", &_nTowers_FHCAL,"tower_FHCAL_N/I");
    _event_tree->Branch("tower_FHCAL_E", _tower_FHCAL_E, "tower_FHCAL_E[tower_FHCAL_N]/F");
    _event_tree->Branch("tower_FHCAL_iEta",_tower_FHCAL_iEta, "tower_FHCAL_iEta[tower_FHCAL_N]/I");
    _event_tree->Branch("tower_FHCAL_iPhi",_tower_FHCAL_iPhi, "tower_FHCAL_iPhi[tower_FHCAL_N]/I");
    _event_tree->Branch("tower_FHCAL_trueID", _tower_FHCAL_trueID, "tower_FHCAL_trueID[tower_FHCAL_N]/I");
    if(_do_CLUSTERS){
      // clusters HCAL
      _event_tree->Branch("cluster_FHCAL_N", &_nclusters_FHCAL,"cluster_FHCAL_N/I");
      _event_tree->Branch("cluster_FHCAL_E", _cluster_FHCAL_E, "cluster_FHCAL_E[cluster_FHCAL_N]/F");
      _event_tree->Branch("cluster_FHCAL_Eta",_cluster_FHCAL_Eta, "cluster_FHCAL_Eta[cluster_FHCAL_N]/F");
      _event_tree->Branch("cluster_FHCAL_Phi",_cluster_FHCAL_Phi, "cluster_FHCAL_Phi[cluster_FHCAL_N]/F");
      _event_tree->Branch("cluster_FHCAL_NTower",_cluster_FHCAL_NTower, "cluster_FHCAL_NTower[cluster_FHCAL_N]/I");
      _event_tree->Branch("cluster_FHCAL_trueID", _cluster_FHCAL_trueID, "cluster_FHCAL_trueID[cluster_FHCAL_N]/I");
    }
  }
  if(_do_DRCALO){
    // towers DRCALO
    _event_tree->Branch("tower_DRCALO_N", &_nTowers_DRCALO,"tower_DRCALO_N/I");
    _event_tree->Branch("tower_DRCALO_E", _tower_DRCALO_E, "tower_DRCALO_E[tower_DRCALO_N]/F");
    _event_tree->Branch("tower_DRCALO_iEta",_tower_DRCALO_iEta, "tower_DRCALO_iEta[tower_DRCALO_N]/I");
    _event_tree->Branch("tower_DRCALO_iPhi",_tower_DRCALO_iPhi, "tower_DRCALO_iPhi[tower_DRCALO_N]/I");
    _event_tree->Branch("tower_DRCALO_trueID", _tower_DRCALO_trueID, "tower_DRCALO_trueID[tower_DRCALO_N]/I");
  }
  if(_do_FEMC){
    // towers EMC
    _event_tree->Branch("tower_FEMC_N",  &_nTowers_FEMC,"tower_FEMC_N/I");
    _event_tree->Branch("tower_FEMC_E",  _tower_FEMC_E, "tower_FEMC_E[tower_FEMC_N]/F");
    _event_tree->Branch("tower_FEMC_iEta", _tower_FEMC_iEta, "tower_FEMC_iEta[tower_FEMC_N]/I");
    _event_tree->Branch("tower_FEMC_iPhi", _tower_FEMC_iPhi, "tower_FEMC_iPhi[tower_FEMC_N]/I");
    _event_tree->Branch("tower_FEMC_trueID",  _tower_FEMC_trueID, "tower_FEMC_trueID[tower_FEMC_N]/I");
    if(_do_CLUSTERS){
      // clusters EMC
      _event_tree->Branch("cluster_FEMC_N",  &_nclusters_FEMC,"cluster_FEMC_N/I");
      _event_tree->Branch("cluster_FEMC_E",  _cluster_FEMC_E, "cluster_FEMC_E[cluster_FEMC_N]/F");
      _event_tree->Branch("cluster_FEMC_Eta", _cluster_FEMC_Eta, "cluster_FEMC_Eta[cluster_FEMC_N]/F");
      _event_tree->Branch("cluster_FEMC_Phi", _cluster_FEMC_Phi, "cluster_FEMC_Phi[cluster_FEMC_N]/F");
      _event_tree->Branch("cluster_FEMC_NTower", _cluster_FEMC_NTower, "cluster_FEMC_NTower[cluster_FEMC_N]/I");
      _event_tree->Branch("cluster_FEMC_trueID",  _cluster_FEMC_trueID, "cluster_FEMC_trueID[cluster_FEMC_N]/I");
    }
  }
  if(_do_VERTEX){
    // vertex
    _event_tree->Branch("vertex_x", &_vertex_x, "vertex_x/I");
    _event_tree->Branch("vertex_y", &_vertex_y, "vertex_y/I");
    _event_tree->Branch("vertex_z", &_vertex_z, "vertex_z/I");
  }
  if(_do_MCPARTICLES){
    // MC particles
    _event_tree->Branch("nMCPart", &_nMCPart, "nMCPart/I");
    _event_tree->Branch("mcpart_ID", _mcpart_ID, "mcpart_ID[nMCPart]/F");
    _event_tree->Branch("mcpart_ID_parent", _mcpart_ID_parent, "mcpart_ID_parent[nMCPart]/F");
    _event_tree->Branch("mcpart_PDG", _mcpart_PDG, "mcpart_PDG[nMCPart]/F");
    _event_tree->Branch("mcpart_E", _mcpart_E, "mcpart_E[nMCPart]/F");
    _event_tree->Branch("mcpart_px", _mcpart_px, "mcpart_px[nMCPart]/F");
    _event_tree->Branch("mcpart_py", _mcpart_py, "mcpart_py[nMCPart]/F");
    _event_tree->Branch("mcpart_pz", _mcpart_pz, "mcpart_pz[nMCPart]/F");
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int EventEvaluator::process_event(PHCompositeNode* topNode)
{
  if(Verbosity() > 0)cout << "entered process_event" << endl;
  if(_do_FHCAL){
    if (!_caloevalstackFHCAL)
    {
      _caloevalstackFHCAL = new CaloEvalStack(topNode, "FHCAL");
      _caloevalstackFHCAL->set_strict(_strict);
      _caloevalstackFHCAL->set_verbosity(Verbosity() + 1);
    }
    else
    {
      _caloevalstackFHCAL->next_event(topNode);
    }
  }
  if(_do_DRCALO){
    if (!_caloevalstackDRCALO)
    {
      _caloevalstackDRCALO = new CaloEvalStack(topNode, "DRCALO");
      _caloevalstackDRCALO->set_strict(_strict);
      _caloevalstackDRCALO->set_verbosity(Verbosity() + 1);
    }
    else
    {
      _caloevalstackDRCALO->next_event(topNode);
    }
  }
  if(_do_FEMC){
    if (!_caloevalstackFEMC)
    {
      _caloevalstackFEMC = new CaloEvalStack(topNode, "FEMC");
      _caloevalstackFEMC->set_strict(_strict);
      _caloevalstackFEMC->set_verbosity(Verbosity() + 1);
    }
    else
    {
      _caloevalstackFEMC->next_event(topNode);
    }
  }
  if (Verbosity() > 0) cout << "loaded evalstack" << endl;

  // print what is coming into the code
  printInputInfo(topNode);
  if (Verbosity() > 0) cout << "past printinputinfo" << endl;

  // fill the Evaluator Tree
  fillOutputNtuples(topNode);

  // Print out the ancestry information for this event
  printOutputInfo(topNode);

  ++_ievent;

  return Fun4AllReturnCodes::EVENT_OK;
}

void EventEvaluator::fillOutputNtuples(PHCompositeNode* topNode)
{
  if (Verbosity() > 2) cout << "EventEvaluator::fillOutputNtuples() entered" << endl;


  //----------------------
  // fill the Event Tree
  //----------------------

  //----------------------
  //    VERTEX
  //----------------------
  GlobalVertexMap* vertexmap = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");
  if(_do_VERTEX){
    if (vertexmap)
    {
      if (!vertexmap->empty())
      {
        if (Verbosity() > 0) cout << "saving vertex" << endl;
        GlobalVertex* vertex = (vertexmap->begin()->second);

        _vertex_x = vertex->get_x();
        _vertex_y = vertex->get_y();
        _vertex_z = vertex->get_z();
      }
    }
  }
  //----------------------
  //    HITS
  //----------------------
  if(_do_HITS){
    if (Verbosity() > 0) cout << "saving hits" << endl;
    _nHitsLayers = 0;
    for (int iIndex = 0; iIndex < 60; ++iIndex)
    {
      if (GetProjectionNameFromIndex(iIndex).find("NOTHING")!= std::string::npos ) continue;
      string nodename = "G4HIT_" + GetProjectionNameFromIndex(iIndex);
      PHG4HitContainer* hits = findNode::getClass<PHG4HitContainer>(topNode, nodename);
      if (hits)
      {
        if (Verbosity() > 1) cout << __PRETTY_FUNCTION__ << " number of hits: " << hits->size() << endl;
        PHG4HitContainer::ConstRange hit_range = hits->getHits();
        for (PHG4HitContainer::ConstIterator hit_iter = hit_range.first; hit_iter != hit_range.second; hit_iter++)
        {
          // if(Verbosity() > 0) cout << __PRETTY_FUNCTION__ << " checking hit id " << hit_iter->second->get_trkid() << " against " << track->get_truth_track_id() << endl;
          // if (hit_iter->second->get_trkid() - track->get_truth_track_id() == 0)
          // {
          if (Verbosity() > 1) cout << __PRETTY_FUNCTION__ << " found hit with id " << hit_iter->second->get_trkid() << endl;
          _hits_x[_nHitsLayers] = hit_iter->second->get_x(0);
          _hits_y[_nHitsLayers] = hit_iter->second->get_y(0);
          _hits_z[_nHitsLayers] = hit_iter->second->get_z(0);
          _hits_t[_nHitsLayers] = hit_iter->second->get_t(0);
          _hits_layerID[_nHitsLayers] = iIndex;
          _nHitsLayers++;

          // }
        }
        if (Verbosity() > 0) cout << "saved\t" << _nHitsLayers << "\thits for " << GetProjectionNameFromIndex(iIndex) << endl;
      }
      else
      {
        if (Verbosity() > 0) cout << __PRETTY_FUNCTION__ << " could not find " << nodename << endl;
        continue;
      }
    }
  }
  //----------------------
  //    TOWERS FHCAL
  //----------------------
  if(_do_FHCAL){
    CaloRawTowerEval* towerevalFHCAL = _caloevalstackFHCAL->get_rawtower_eval();
    _nTowers_FHCAL = 0;
    string towernodeFHCAL = "TOWER_CALIB_FHCAL";
    RawTowerContainer* towersFHCAL = findNode::getClass<RawTowerContainer>(topNode, towernodeFHCAL.c_str());
    if (towersFHCAL)
    {
      if (Verbosity() > 0) cout << "saving HCAL towers" << endl;
      string towergeomnodeFHCAL = "TOWERGEOM_FHCAL";
      RawTowerGeomContainer* towergeomFHCAL = findNode::getClass<RawTowerGeomContainer>(topNode, towergeomnodeFHCAL.c_str());
      if (towergeomFHCAL)
      {
        RawTowerContainer::ConstRange begin_end = towersFHCAL->getTowers();
        RawTowerContainer::ConstIterator rtiter;
        for (rtiter = begin_end.first; rtiter != begin_end.second; ++rtiter)
        {
          RawTower* tower = rtiter->second;
          if (tower)
          {
            // min energy cut
            if (tower->get_energy() < _reco_e_threshold) continue;

            _tower_FHCAL_iEta[_nTowers_FHCAL] = tower->get_bineta();
            _tower_FHCAL_iPhi[_nTowers_FHCAL] = tower->get_binphi();
            _tower_FHCAL_E[_nTowers_FHCAL] = tower->get_energy();

            PHG4Particle* primary = towerevalFHCAL->max_truth_primary_particle_by_energy(tower);
            if (primary)
            {
              _tower_FHCAL_trueID[_nTowers_FHCAL] = primary->get_track_id();
              // gflavor = primary->get_pid();
              // efromtruth = towerevalFHCAL->get_energy_contribution(tower, primary);
            } else {
              _tower_FHCAL_trueID[_nTowers_FHCAL] = -10;
            }
            _nTowers_FHCAL++;
          }
        }
      }
      else
      {
        if (Verbosity() > 0) cout << PHWHERE << " ERROR: Can't find " << towergeomnodeFHCAL << endl;
        // return;
      }
      if (Verbosity() > 0) cout << "saved\t" << _nTowers_FHCAL << "\tFHCAL towers" << endl;
    }
    else
    {
      if (Verbosity() > 0) cout << PHWHERE << " ERROR: Can't find " << towernodeFHCAL << endl;
      // return;
    }
  }
  //----------------------
  //    TOWERS DRCALO
  //----------------------
  if(_do_DRCALO){
    CaloRawTowerEval* towerevalDRCALO = _caloevalstackDRCALO->get_rawtower_eval();
    _nTowers_DRCALO = 0;
    string towernodeDRCALO = "TOWER_CALIB_DRCALO";
    RawTowerContainer* towersDRCALO = findNode::getClass<RawTowerContainer>(topNode, towernodeDRCALO.c_str());
    if (towersDRCALO)
    {
      if (Verbosity() > 0) cout << "saving DRCALO towers" << endl;
      string towergeomnodeDRCALO = "TOWERGEOM_DRCALO";
      RawTowerGeomContainer* towergeomDRCALO = findNode::getClass<RawTowerGeomContainer>(topNode, towergeomnodeDRCALO.c_str());
      if (towergeomDRCALO)
      {
        RawTowerContainer::ConstRange begin_end = towersDRCALO->getTowers();
        RawTowerContainer::ConstIterator rtiter;
        for (rtiter = begin_end.first; rtiter != begin_end.second; ++rtiter)
        {
          RawTower* tower = rtiter->second;
          if (tower)
          {
            // min energy cut
            if (tower->get_energy() < _reco_e_threshold) continue;

            _tower_DRCALO_iEta[_nTowers_DRCALO] = tower->get_bineta();
            _tower_DRCALO_iPhi[_nTowers_DRCALO] = tower->get_binphi();
            _tower_DRCALO_E[_nTowers_DRCALO] = tower->get_energy();

            PHG4Particle* primary = towerevalDRCALO->max_truth_primary_particle_by_energy(tower);
            if (primary)
            {
              _tower_DRCALO_trueID[_nTowers_DRCALO] = primary->get_track_id();
              // gflavor = primary->get_pid();
              // efromtruth = towerevalDRCALO->get_energy_contribution(tower, primary);
            } else {
              _tower_DRCALO_trueID[_nTowers_DRCALO] = -10;
            }
            _nTowers_DRCALO++;
          }
        }
      }
      else
      {
        if (Verbosity() > 0) cout << PHWHERE << " ERROR: Can't find " << towergeomnodeDRCALO << endl;
        // return;
      }
      if (Verbosity() > 0) cout << "saved\t" << _nTowers_DRCALO << "\tDRCALO towers" << endl;
    }
    else
    {
      if (Verbosity() > 0) cout << PHWHERE << " ERROR: Can't find " << towernodeDRCALO << endl;
      // return;
    }
  }
  //----------------------
  //    TOWERS
  //----------------------
  if(_do_FEMC){
    CaloRawTowerEval* towerevalFEMC = _caloevalstackFEMC->get_rawtower_eval();
    _nTowers_FEMC = 0;
    string towernodeFEMC = "TOWER_CALIB_FEMC";
    RawTowerContainer* towersFEMC = findNode::getClass<RawTowerContainer>(topNode, towernodeFEMC.c_str());
    if (towersFEMC)
    {
      if (Verbosity() > 0) cout << "saving EMC towers" << endl;
      string towergeomnodeFEMC = "TOWERGEOM_FEMC";
      RawTowerGeomContainer* towergeom = findNode::getClass<RawTowerGeomContainer>(topNode, towergeomnodeFEMC.c_str());
      if (towergeom)
      {
        RawTowerContainer::ConstRange begin_end = towersFEMC->getTowers();
        RawTowerContainer::ConstIterator rtiter;
        for (rtiter = begin_end.first; rtiter != begin_end.second; ++rtiter)
        {
          RawTower* tower = rtiter->second;
          if (tower)
          {
            // min energy cut
            if (tower->get_energy() < _reco_e_threshold) continue;

            _tower_FEMC_iEta[_nTowers_FEMC] = tower->get_bineta();
            _tower_FEMC_iPhi[_nTowers_FEMC] = tower->get_binphi();
            _tower_FEMC_E[_nTowers_FEMC] = tower->get_energy();

            PHG4Particle* primary = towerevalFEMC->max_truth_primary_particle_by_energy(tower);
            if (primary)
            {
              _tower_FEMC_trueID[_nTowers_FEMC] = primary->get_track_id();
              // gflavor = primary->get_pid();
              // efromtruth = towerevalFEMC->get_energy_contribution(tower, primary);
            } else {
              _tower_FEMC_trueID[_nTowers_FEMC] = -10;
            }
            _nTowers_FEMC++;
          }
        }
      }
      else
      {
        if (Verbosity() > 0) cout << PHWHERE << " ERROR: Can't find " << towergeomnodeFEMC << endl;
        // return;
      }
      if (Verbosity() > 0) cout << "saved\t" << _nTowers_FEMC << "\tFEMC towers" << endl;
    }
    else
    {
      if (Verbosity() > 0) cout << PHWHERE << " ERROR: Can't find " << towernodeFEMC << endl;
      // return;
    }
  }
  //------------------------
  // CLUSTERS FHCAL
  //------------------------
  if(_do_FHCAL && _do_CLUSTERS){
    CaloRawClusterEval* clusterevalFHCAL = _caloevalstackFHCAL->get_rawcluster_eval();
    _nclusters_FHCAL = 0;
    if (Verbosity() > 1) cout << "CaloEvaluator::filling gcluster ntuple..." << endl;

    string clusternodeFHCAL = "CLUSTER_FHCAL";
    RawClusterContainer* clustersFHCAL = findNode::getClass<RawClusterContainer>(topNode, clusternodeFHCAL.c_str());
    if (clustersFHCAL)
    {
      // for every cluster
      for (const auto& iterator : clustersFHCAL->getClustersMap())
      {
        RawCluster* cluster = iterator.second;

        if (cluster->get_energy() < _reco_e_threshold) continue;

        _cluster_FHCAL_E[_nclusters_FHCAL] = cluster->get_energy();
        _cluster_FHCAL_NTower[_nclusters_FHCAL] = cluster->getNTowers();
        _cluster_FHCAL_Phi[_nclusters_FHCAL] = cluster->get_phi();

        // require vertex for cluster eta calculation
        if (vertexmap)
        {
          if (!vertexmap->empty())
          {
            GlobalVertex* vertex = (vertexmap->begin()->second);
            _cluster_FHCAL_Eta[_nclusters_FHCAL] = RawClusterUtility::GetPseudorapidity(*cluster, CLHEP::Hep3Vector(vertex->get_x(), vertex->get_y(), vertex->get_z()));
          }
          else
            _cluster_FHCAL_Eta[_nclusters_FHCAL] = -10000;
        }
        else
          _cluster_FHCAL_Eta[_nclusters_FHCAL] = -10000;

        PHG4Particle* primary = clusterevalFHCAL->max_truth_primary_particle_by_energy(cluster);

        if (primary){
          _cluster_FHCAL_trueID[_nclusters_FHCAL]       = primary->get_track_id();
        } else {
          _cluster_FHCAL_trueID[_nclusters_FHCAL]       = -10;
        }

        _nclusters_FHCAL++;
      }
    }
    else
    {
      cerr << PHWHERE << " ERROR: Can't find " << clusternodeFHCAL << endl;
      // return;
    }
  }
  //------------------------
  // CLUSTERS FEMC
  //------------------------
  if(_do_FEMC && _do_CLUSTERS){
    CaloRawClusterEval* clusterevalFEMC = _caloevalstackFEMC->get_rawcluster_eval();
    _nclusters_FEMC = 0;
    if (Verbosity() > 1) cout << "CaloEvaluator::filling gcluster ntuple..." << endl;

    string clusternodeFEMC = "CLUSTER_FEMC";
    RawClusterContainer* clustersFEMC = findNode::getClass<RawClusterContainer>(topNode, clusternodeFEMC.c_str());
    if (clustersFEMC)
    {
      // for every cluster
      for (const auto& iterator : clustersFEMC->getClustersMap())
      {
        RawCluster* cluster = iterator.second;

        if (cluster->get_energy() < _reco_e_threshold) continue;

        _cluster_FEMC_E[_nclusters_FEMC] = cluster->get_energy();
        _cluster_FEMC_NTower[_nclusters_FEMC] = cluster->getNTowers();
        _cluster_FEMC_Phi[_nclusters_FEMC] = cluster->get_phi();

        // require vertex for cluster eta calculation
        if (vertexmap)
        {
          if (!vertexmap->empty())
          {
            GlobalVertex* vertex = (vertexmap->begin()->second);
            _cluster_FEMC_Eta[_nclusters_FEMC] = RawClusterUtility::GetPseudorapidity(*cluster, CLHEP::Hep3Vector(vertex->get_x(), vertex->get_y(), vertex->get_z()));
          }
          else
            _cluster_FEMC_Eta[_nclusters_FEMC] = -10000;
        }
        else
          _cluster_FEMC_Eta[_nclusters_FEMC] = -10000;

        PHG4Particle* primary = clusterevalFEMC->max_truth_primary_particle_by_energy(cluster);

        if (primary){
          _cluster_FEMC_trueID[_nclusters_FEMC]       = primary->get_track_id(); 
        } else {
          _cluster_FEMC_trueID[_nclusters_FEMC] = -10;
        }

        _nclusters_FEMC++;
      }
    }
    else
    {
      cerr << PHWHERE << " ERROR: Can't find " << clusternodeFEMC << endl;
      // return;
    }
  }
  //------------------------
  // TRACKS
  //------------------------
  if(_do_TRACKS){
    _nTracks = 0;
    _nProjections = 0;
    SvtxTrackMap* trackmap = findNode::getClass<SvtxTrackMap>(topNode, "TrackMap");
    if (trackmap)
    {
      if (Verbosity() > 0) cout << "saving tracks" << endl;
      for (SvtxTrackMap::ConstIter track_itr = trackmap->begin(); track_itr != trackmap->end(); track_itr++)
      {
        SvtxTrack_FastSim* track = dynamic_cast<SvtxTrack_FastSim*>(track_itr->second);
        if (track)
        {
          _track_ID[_nTracks] = track->get_id();
          _track_px[_nTracks] = track->get_px();
          _track_py[_nTracks] = track->get_py();
          _track_pz[_nTracks] = track->get_pz();
          _track_trueID[_nTracks] = track->get_truth_track_id();
          if(_do_PROJECTIONS){
            // find projections
            for (SvtxTrack::ConstStateIter trkstates = track->begin_states(); trkstates != track->end_states(); ++trkstates)
            {
              if (Verbosity() > 1) cout << __PRETTY_FUNCTION__ << " processing " << trkstates->second->get_name() << endl;
              string trackStateName = trkstates->second->get_name();
              if (Verbosity() > 1) cout << __PRETTY_FUNCTION__ << " found " << trkstates->second->get_name() << endl;
              int trackStateIndex = GetProjectionIndex(trackStateName);
              if (trackStateIndex > -1)
              {
                // save true projection info to given branch
                _track_TLP_true_x[_nProjections] = trkstates->second->get_pos(0);
                _track_TLP_true_y[_nProjections] = trkstates->second->get_pos(1);
                _track_TLP_true_z[_nProjections] = trkstates->second->get_pos(2);
                _track_TLP_true_t[_nProjections] = trkstates->first;
                _track_ProjLayer[_nProjections] = trackStateIndex;
                _track_ProjTrackID[_nProjections] = _nTracks;

                string nodename = "G4HIT_" + trkstates->second->get_name();
                PHG4HitContainer* hits = findNode::getClass<PHG4HitContainer>(topNode, nodename);
                if (hits)
                {
                  if (Verbosity() > 1) cout << __PRETTY_FUNCTION__ << " number of hits: " << hits->size() << endl;
                  PHG4HitContainer::ConstRange hit_range = hits->getHits();
                  for (PHG4HitContainer::ConstIterator hit_iter = hit_range.first; hit_iter != hit_range.second; hit_iter++)
                  {
                    if (Verbosity() > 1) cout << __PRETTY_FUNCTION__ << " checking hit id " << hit_iter->second->get_trkid() << " against " << track->get_truth_track_id() << endl;
                    if (hit_iter->second->get_trkid() - track->get_truth_track_id() == 0)
                    {
                      if (Verbosity() > 1) cout << __PRETTY_FUNCTION__ << " found hit with id " << hit_iter->second->get_trkid() << endl;
                      // save reco projection info to given branch
                      _track_TLP_x[_nProjections] = hit_iter->second->get_x(0);
                      _track_TLP_y[_nProjections] = hit_iter->second->get_y(0);
                      _track_TLP_z[_nProjections] = hit_iter->second->get_z(0);
                      _track_TLP_t[_nProjections] = hit_iter->second->get_t(0);
                    }
                  }
                }
                else
                {
                  if (Verbosity() > 1) cout << __PRETTY_FUNCTION__ << " could not find " << nodename << endl;
                  continue;
                }
                _nProjections++;
              }
            }
          }
          _nTracks++;
        }
        else
        {
          if (Verbosity() > 0)  //Verbosity()
          {
            cout << "PHG4TrackFastSimEval::fill_track_tree - ignore track that is not a SvtxTrack_FastSim:";
            track_itr->second->identify();
          }
          continue;
        }
      }
      if (Verbosity() > 0) cout << "saved\t" << _nTracks << "\ttracks" << endl;
    }
    else
    {
      if (Verbosity() > 0) cout << PHWHERE << "SvtxTrackMap node with name TrackMap not found on node tree" << endl;
      return;
    }
  }
  //------------------------
  // MC PARTICLES
  //------------------------
  _nMCPart = 0;
  if(_do_MCPARTICLES){
    PHG4TruthInfoContainer* truthinfocontainer = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
    if (truthinfocontainer)
    {
      if (Verbosity() > 0) cout << "saving MC particles" << endl;
      //GetParticleRange for all particles
      //GetPrimaryParticleRange for primary particles
      PHG4TruthInfoContainer::ConstRange range = truthinfocontainer->GetPrimaryParticleRange();
      for (PHG4TruthInfoContainer::ConstIterator truth_itr = range.first; truth_itr != range.second; ++truth_itr)
      {
        PHG4Particle* g4particle = truth_itr->second;
        if (!g4particle) continue;
        // in case of all MC particles, make restrictions on the secondary selection
        // if(g4particle->get_track_id()<0 && g4particle->get_e()<0.5) continue;
        // primary (g4particle->get_parent_id() == 0) selection via:
        // if(gtrackID < 0) continue;

        _mcpart_ID[_nMCPart] = g4particle->get_track_id();
        _mcpart_ID_parent[_nMCPart] = g4particle->get_parent_id();
        _mcpart_PDG[_nMCPart] = g4particle->get_pid();
        _mcpart_E[_nMCPart] = g4particle->get_e();
        _mcpart_px[_nMCPart] = g4particle->get_px();
        _mcpart_py[_nMCPart] = g4particle->get_py();
        _mcpart_pz[_nMCPart] = g4particle->get_pz();
        _nMCPart++;
      }
      if (Verbosity() > 0) cout << "saved\t" << _nMCPart << "\tMC particles" << endl;
    }
    else
    {
      if (Verbosity() > 0) cout << PHWHERE << " PHG4TruthInfoContainer node not found on node tree" << endl;
      return;
    }
  }

  _event_tree->Fill();
  resetBuffer();
  return;
}

void EventEvaluator::printOutputInfo(PHCompositeNode* topNode)
{
  if (Verbosity() > 2) cout << "EventEvaluator::printOutputInfo() entered" << endl;

  CaloRawClusterEval* clustereval = _caloevalstackFHCAL->get_rawcluster_eval();
  CaloTruthEval* trutheval = _caloevalstackFHCAL->get_truth_eval();

  //==========================================
  // print out some useful stuff for debugging
  //==========================================

  if (Verbosity() > 1)
  {
    // event information
    cout << endl;
    cout << PHWHERE << "   NEW OUTPUT FOR EVENT " << _ievent << endl;
    cout << endl;

    // need things off of the DST...
    PHG4TruthInfoContainer* truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
    if (!truthinfo)
    {
      cout << PHWHERE << " ERROR: Can't find G4TruthInfo" << endl;
      return;
    }

    // need things off of the DST...
    GlobalVertexMap* vertexmap = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");

    PHG4VtxPoint* gvertex = truthinfo->GetPrimaryVtx(truthinfo->GetPrimaryVertexIndex());
    float gvx = gvertex->get_x();
    float gvy = gvertex->get_y();
    float gvz = gvertex->get_z();

    float vx = NAN;
    float vy = NAN;
    float vz = NAN;
    if (vertexmap)
    {
      if (!vertexmap->empty())
      {
        GlobalVertex* vertex = (vertexmap->begin()->second);

        vx = vertex->get_x();
        vy = vertex->get_y();
        vz = vertex->get_z();
      }
    }

    cout << "vtrue = (" << gvx << "," << gvy << "," << gvz << ") => vreco = (" << vx << "," << vy << "," << vz << ")" << endl;

    PHG4TruthInfoContainer::ConstRange range = truthinfo->GetPrimaryParticleRange();
    for (PHG4TruthInfoContainer::ConstIterator iter = range.first;
         iter != range.second;
         ++iter)
    {
      PHG4Particle* primary = iter->second;

      cout << endl;

      cout << "===Primary PHG4Particle=========================================" << endl;
      cout << " particle id = " << primary->get_track_id() << endl;
      cout << " flavor = " << primary->get_pid() << endl;
      cout << " (px,py,pz,e) = (";

      float gpx = primary->get_px();
      float gpy = primary->get_py();
      float gpz = primary->get_pz();
      float ge = primary->get_e();

      cout.width(5);
      cout << gpx;
      cout << ",";
      cout.width(5);
      cout << gpy;
      cout << ",";
      cout.width(5);
      cout << gpz;
      cout << ",";
      cout.width(5);
      cout << ge;
      cout << ")" << endl;

      float gpt = sqrt(gpx * gpx + gpy * gpy);
      float geta = NAN;
      if (gpt != 0.0) geta = asinh(gpz / gpt);
      float gphi = atan2(gpy, gpx);

      cout << "(eta,phi,e,pt) = (";
      cout.width(5);
      cout << geta;
      cout << ",";
      cout.width(5);
      cout << gphi;
      cout << ",";
      cout.width(5);
      cout << ge;
      cout << ",";
      cout.width(5);
      cout << gpt;
      cout << ")" << endl;

      PHG4VtxPoint* vtx = trutheval->get_vertex(primary);

      cout << " vtrue = (";
      cout.width(5);
      cout << vtx->get_x();
      cout << ",";
      cout.width(5);
      cout << vtx->get_y();
      cout << ",";
      cout.width(5);
      cout << vtx->get_z();
      cout << ")" << endl;

      cout << " embed = " << trutheval->get_embed(primary) << endl;
      cout << " edep = " << trutheval->get_shower_energy_deposit(primary) << endl;

      std::set<RawCluster*> clusters = clustereval->all_clusters_from(primary);
      for (std::set<RawCluster*>::iterator clusiter = clusters.begin();
           clusiter != clusters.end();
           ++clusiter)
      {
        RawCluster* cluster = (*clusiter);

        float ntowers = cluster->getNTowers();
        float x = cluster->get_x();
        float y = cluster->get_y();
        float z = cluster->get_z();
        float phi = cluster->get_phi();
        float e = cluster->get_energy();

        float efromtruth = clustereval->get_energy_contribution(cluster, primary);

        cout << " => #" << cluster->get_id() << " (x,y,z,phi,e) = (";
        cout.width(5);
        cout << x;
        cout << ",";
        cout.width(5);
        cout << y;
        cout << ",";
        cout.width(5);
        cout << z;
        cout << ",";
        cout.width(5);
        cout << phi;
        cout << ",";
        cout.width(5);
        cout << e;
        cout << "), ntowers = " << ntowers << ", efromtruth = " << efromtruth << endl;
      }
    }
    cout << endl;
  }

  return;
}

void EventEvaluator::printInputInfo(PHCompositeNode* topNode)
{
  if (Verbosity() > 2) cout << "EventEvaluator::printInputInfo() entered" << endl;

  // print out the truth container

  if (Verbosity() > 1)
  {
    cout << endl;
    cout << PHWHERE << "   NEW INPUT FOR EVENT " << _ievent << endl;
    cout << endl;

    // need things off of the DST...
    PHG4TruthInfoContainer* truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
    if (!truthinfo)
    {
      cout << PHWHERE << " ERROR: Can't find G4TruthInfo" << endl;
      return;
    }

    cout << Name() << ": PHG4TruthInfoContainer contents: " << endl;

    PHG4TruthInfoContainer::Range truthrange = truthinfo->GetParticleRange();
    for (PHG4TruthInfoContainer::Iterator truthiter = truthrange.first;
         truthiter != truthrange.second;
         ++truthiter)
    {
      PHG4Particle* particle = truthiter->second;

      cout << truthiter->first << " => pid: " << particle->get_pid()
           << " pt: " << sqrt(pow(particle->get_px(), 2) + pow(particle->get_py(), 2)) << endl;
    }
  }

  return;
}

int EventEvaluator::End(PHCompositeNode* topNode)
{
  _tfile->cd();

  _event_tree->Write();

  _tfile->Close();

  delete _tfile;

  if (Verbosity() > 0)
  {
    cout << "========================= " << Name() << "::End() ============================" << endl;
    cout << " " << _ievent << " events of output written to: " << _filename << endl;
    cout << "===========================================================================" << endl;
  }

  if (_caloevalstackFHCAL) delete _caloevalstackFHCAL;
  if (_caloevalstackFEMC) delete _caloevalstackFEMC;

  return Fun4AllReturnCodes::EVENT_OK;
}

int EventEvaluator::GetProjectionIndex(std::string projname)
{
  if (projname.find("FTTL_0") != std::string::npos)
    return 0;
  else if (projname.find("FTTL_1") != std::string::npos)
    return 1;
  else if (projname.find("FTTL_2") != std::string::npos)
    return 2;
  else if (projname.find("ETTL_0") != std::string::npos)
    return 3;
  else if (projname.find("ETTL_1") != std::string::npos)
    return 4;
  else if (projname.find("FHCAL_0") != std::string::npos)
    return 5;
  else if (projname.find("FEMC_0") != std::string::npos)
    return 6;
  else if (projname.find("CTTL_0") != std::string::npos)
    return 7;
  else if (projname.find("CTTL_1") != std::string::npos)
    return 8;

  else if (projname.find("LBLVTX_CENTRAL_10") != std::string::npos)
    return 10;
  else if (projname.find("LBLVTX_CENTRAL_11") != std::string::npos)
    return 11;
  else if (projname.find("LBLVTX_CENTRAL_12") != std::string::npos)
    return 12;
  else if (projname.find("LBLVTX_CENTRAL_13") != std::string::npos)
    return 13;
  else if (projname.find("LBLVTX_CENTRAL_14") != std::string::npos)
    return 14;
  else if (projname.find("LBLVTX_CENTRAL_15") != std::string::npos)
    return 15;
              
  else if (projname.find("LBLVTX_FORWARD_20") != std::string::npos)
    return 20;
  else if (projname.find("LBLVTX_FORWARD_21") != std::string::npos)
    return 21;
  else if (projname.find("LBLVTX_FORWARD_22") != std::string::npos)
    return 22;
  else if (projname.find("LBLVTX_FORWARD_23") != std::string::npos)
    return 23;
  else if (projname.find("LBLVTX_FORWARD_24") != std::string::npos)
    return 24;
            
  else if (projname.find("LBLVTX_BACKWARD_30") != std::string::npos)
    return 30;
  else if (projname.find("LBLVTX_BACKWARD_31") != std::string::npos)
    return 31;
  else if (projname.find("LBLVTX_BACKWARD_32") != std::string::npos)
    return 32;
  else if (projname.find("LBLVTX_BACKWARD_33") != std::string::npos)
    return 33;
  else if (projname.find("LBLVTX_BACKWARD_34") != std::string::npos)
    return 34;
            
  else if (projname.find("BARREL_0") != std::string::npos)
    return 40;
  else if (projname.find("BARREL_1") != std::string::npos)
    return 41;
  else if (projname.find("BARREL_2") != std::string::npos)
    return 42;
  else if (projname.find("BARREL_3") != std::string::npos)
    return 43;
  else if (projname.find("BARREL_4") != std::string::npos)
    return 44;
  else if (projname.find("BARREL_5") != std::string::npos)
    return 45;

  else if (projname.find("FST_0") != std::string::npos)
    return 50;
  else if (projname.find("FST_1") != std::string::npos)
    return 51;
  else if (projname.find("FST_2") != std::string::npos)
    return 52;
  else if (projname.find("FST_3") != std::string::npos)
    return 53;
  else if (projname.find("FST_4") != std::string::npos)
    return 54;
  else if (projname.find("FST_5") != std::string::npos)
    return 55;

  else
    return -1;
  return -1;
}

std::string EventEvaluator::GetProjectionNameFromIndex(int projindex)
{
  switch (projindex)
  {
  case 0:
    return "FTTL_0";
  case 1:
    return "FTTL_1";
  case 2:
    return "FTTL_2";
  case 3:
    return "ETTL_0";
  case 4:
    return "ETTL_1";
  case 5:
    return "FHCAL";
  case 6:
    return "FEMC";
  case 7:
    return "CTTL_0";
  case 8:
    return "CTTL_1";


  case 10:
    return "LBLVTX_CENTRAL_10";
  case 11:
    return "LBLVTX_CENTRAL_11";
  case 12:
    return "LBLVTX_CENTRAL_12";
  case 13:
    return "LBLVTX_CENTRAL_13";
  case 14:
    return "LBLVTX_CENTRAL_14";
  case 15:
    return "LBLVTX_CENTRAL_15";

  case 20:
    return "LBLVTX_FORWARD_20";
  case 21:
    return "LBLVTX_FORWARD_21";
  case 22:
    return "LBLVTX_FORWARD_22";
  case 23: 
    return "LBLVTX_FORWARD_23";
  case 24:
    return "LBLVTX_FORWARD_24";
    
  case 30:
    return "LBLVTX_BACKWARD_30";
  case 31: 
    return "LBLVTX_BACKWARD_31";
  case 32:
    return "LBLVTX_BACKWARD_32";
  case 33:
    return "LBLVTX_BACKWARD_33";
  case 34:  
    return "LBLVTX_BACKWARD_34";
            
  case 40:
    return "BARREL_0";
  case 41:
    return "BARREL_1";
  case 42:
    return "BARREL_2";
  case 43:
    return "BARREL_3";
  case 44: 
    return "BARREL_4";
  case 45:
    return "BARREL_5";

  case 50:
    return "FST_0";
  case 51:
    return "FST_1";
  case 52:
    return "FST_2";
  case 53:
    return "FST_3";
  case 54:
    return "FST_4";
  case 55:
    return "FST_5";
  default:
    return "NOTHING";
  }
}

void EventEvaluator::resetBuffer()
{
  _vertex_x = 0;
  _vertex_y = 0;
  _vertex_z = 0;

  if(_do_HITS){
    _nHitsLayers = 0;
    for (Int_t ihit = 0; ihit < _maxNHits; ihit++)
    {
      _hits_layerID[ihit] = 0;
      _hits_x[ihit] = 0;
      _hits_y[ihit] = 0;
      _hits_z[ihit] = 0;
      _hits_t[ihit] = 0;
    }
  }
  if(_do_FHCAL){
    _nTowers_FHCAL = 0;
    for (Int_t itow = 0; itow < _maxNTowers; itow++)
    {
      _tower_FHCAL_E[itow] = 0;
      _tower_FHCAL_iEta[itow] = 0;
      _tower_FHCAL_iPhi[itow] = 0;
      _tower_FHCAL_trueID[itow] = 0;
    }
    if(_do_CLUSTERS){
      _nclusters_FHCAL = 0;
      for (Int_t itow = 0; itow < _maxNclusters; itow++)
      {
        _cluster_FHCAL_E[itow] = 0;
        _cluster_FHCAL_Eta[itow] = 0;
        _cluster_FHCAL_Phi[itow] = 0;
        _cluster_FHCAL_NTower[itow] = 0;
        _cluster_FHCAL_trueID[itow] = 0;
      }
    }
  }
  if(_do_FEMC){
    _nTowers_FEMC = 0;
    for (Int_t itow = 0; itow < _maxNTowers; itow++)
    {
      _tower_FEMC_E[itow] = 0;
      _tower_FEMC_iEta[itow] = 0;
      _tower_FEMC_iPhi[itow] = 0;
      _tower_FEMC_trueID[itow] = 0;
    }
    if(_do_CLUSTERS){
      _nclusters_FEMC = 0;
      for (Int_t itow = 0; itow < _maxNclusters; itow++)
      {
        _cluster_FEMC_E[itow] = 0;
        _cluster_FEMC_Eta[itow] = 0;
        _cluster_FEMC_Phi[itow] = 0;
        _cluster_FEMC_NTower[itow] = 0;
        _cluster_FEMC_trueID[itow] = 0;
      }
    }
  }

  _nTracks = 0;
  for (Int_t itrk = 0; itrk < _maxNTracks; itrk++)
  {
    _track_ID[itrk] = 0;
    _track_trueID[itrk] = 0;
    _track_px[itrk] = 0;
    _track_py[itrk] = 0;
    _track_pz[itrk] = 0;
  }
  _nProjections = 0;
  for (Int_t iproj = 0; iproj < _maxNProjections; iproj++)
  {
    _track_ProjLayer[iproj] = -1;
    _track_ProjTrackID[iproj] = 0;
    _track_TLP_x[iproj] = 0;
    _track_TLP_y[iproj] = 0;
    _track_TLP_z[iproj] = 0;
    _track_TLP_t[iproj] = 0;
    _track_TLP_true_x[iproj] = 0;
    _track_TLP_true_y[iproj] = 0;
    _track_TLP_true_z[iproj] = 0;
    _track_TLP_true_t[iproj] = 0;
  }
  _nMCPart = 0;
  for (Int_t imcpart = 0; imcpart < _maxNMCPart; imcpart++)
  {
    _mcpart_ID[imcpart] = 0;
    _mcpart_ID_parent[imcpart] = 0;
    _mcpart_PDG[imcpart] = 0;
    _mcpart_E[imcpart] = 0;
    _mcpart_px[imcpart] = 0;
    _mcpart_py[imcpart] = 0;
    _mcpart_pz[imcpart] = 0;
  }
}
