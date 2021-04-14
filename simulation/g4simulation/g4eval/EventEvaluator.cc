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
#include <trackbase_historic/SvtxVertex.h>  // for SvtxVertex
#include <trackbase_historic/SvtxVertexMap.h>

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
  , _do_FHCAL(false)
  , _do_HCALIN(false)
  , _do_HCALOUT(false)
  , _do_EHCAL(false)
  , _do_FEMC(false)
  , _do_CEMC(false)
  , _do_EEMC(false)
  , _do_DRCALO(false)
  , _do_HITS(false)
  , _do_TRACKS(false)
  , _do_CLUSTERS(false)
  , _do_VERTEX(false)
  , _do_PROJECTIONS(false)
  , _do_MCPARTICLES(false)
  , _ievent(0)
  , _nHitsLayers(0)
  , _hits_layerID(0)
  , _hits_trueID(0)
  , _hits_x(0)
  , _hits_y(0)
  , _hits_z(0)
  , _hits_t(0)

  , _nTowers_FHCAL(0)
  , _tower_FHCAL_E(0)
  , _tower_FHCAL_iEta(0)
  , _tower_FHCAL_iPhi(0)
  , _tower_FHCAL_trueID(0)

  , _nTowers_HCALIN(0)
  , _tower_HCALIN_E(0)
  , _tower_HCALIN_iEta(0)
  , _tower_HCALIN_iPhi(0)
  , _tower_HCALIN_trueID(0)

  , _nTowers_HCALOUT(0)
  , _tower_HCALOUT_E(0)
  , _tower_HCALOUT_iEta(0)
  , _tower_HCALOUT_iPhi(0)
  , _tower_HCALOUT_trueID(0)

  , _nTowers_EHCAL(0)
  , _tower_EHCAL_E(0)
  , _tower_EHCAL_iEta(0)
  , _tower_EHCAL_iPhi(0)
  , _tower_EHCAL_trueID(0)
  
  , _nTowers_DRCALO(0)
  , _tower_DRCALO_E(0)
  , _tower_DRCALO_NScint(0)
  , _tower_DRCALO_NCerenkov(0)
  , _tower_DRCALO_iEta(0)
  , _tower_DRCALO_iPhi(0)
  , _tower_DRCALO_trueID(0)

  , _nTowers_FEMC(0)
  , _tower_FEMC_E(0)
  , _tower_FEMC_iEta(0)
  , _tower_FEMC_iPhi(0)
  , _tower_FEMC_trueID(0)

  , _nTowers_CEMC(0)
  , _tower_CEMC_E(0)
  , _tower_CEMC_iEta(0)
  , _tower_CEMC_iPhi(0)
  , _tower_CEMC_trueID(0)

  , _nTowers_EEMC(0)
  , _tower_EEMC_E(0)
  , _tower_EEMC_iEta(0)
  , _tower_EEMC_iPhi(0)
  , _tower_EEMC_trueID(0)

  , _nclusters_FHCAL(0)
  , _cluster_FHCAL_E(0)
  , _cluster_FHCAL_Eta(0)
  , _cluster_FHCAL_Phi(0)
  , _cluster_FHCAL_NTower(0)
  , _cluster_FHCAL_trueID(0)

  , _nclusters_HCALIN(0)
  , _cluster_HCALIN_E(0)
  , _cluster_HCALIN_Eta(0)
  , _cluster_HCALIN_Phi(0)
  , _cluster_HCALIN_NTower(0)
  , _cluster_HCALIN_trueID(0)

  , _nclusters_HCALOUT(0)
  , _cluster_HCALOUT_E(0)
  , _cluster_HCALOUT_Eta(0)
  , _cluster_HCALOUT_Phi(0)
  , _cluster_HCALOUT_NTower(0)
  , _cluster_HCALOUT_trueID(0)

  , _nclusters_EHCAL(0)
  , _cluster_EHCAL_E(0)
  , _cluster_EHCAL_Eta(0)
  , _cluster_EHCAL_Phi(0)
  , _cluster_EHCAL_NTower(0)
  , _cluster_EHCAL_trueID(0)

  , _nclusters_FEMC(0)
  , _cluster_FEMC_E(0)
  , _cluster_FEMC_Eta(0)
  , _cluster_FEMC_Phi(0)
  , _cluster_FEMC_NTower(0)
  , _cluster_FEMC_trueID(0)

  , _nclusters_CEMC(0)
  , _cluster_CEMC_E(0)
  , _cluster_CEMC_Eta(0)
  , _cluster_CEMC_Phi(0)
  , _cluster_CEMC_NTower(0)
  , _cluster_CEMC_trueID(0)

  , _nclusters_EEMC(0)
  , _cluster_EEMC_E(0)
  , _cluster_EEMC_Eta(0)
  , _cluster_EEMC_Phi(0)
  , _cluster_EEMC_NTower(0)
  , _cluster_EEMC_trueID(0)
  
  , _vertex_x(0)
  , _vertex_y(0)
  , _vertex_z(0)
  , _vertex_NCont(0)
  , _vertex_true_x(0)
  , _vertex_true_y(0)
  , _vertex_true_z(0)
  
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
  , _depth_MCstack(3)
  , _caloevalstackFHCAL(nullptr)
  , _caloevalstackHCALIN(nullptr)
  , _caloevalstackHCALOUT(nullptr)
  , _caloevalstackEHCAL(nullptr)
  , _caloevalstackDRCALO(nullptr)
  , _caloevalstackFEMC(nullptr)
  , _caloevalstackCEMC(nullptr)
  , _caloevalstackEEMC(nullptr)
  , _strict(false)
  , _event_tree(nullptr)
  , _filename(filename)
  , _tfile(nullptr)
{
  _hits_layerID = new int[_maxNHits];
  _hits_trueID = new int[_maxNHits];
  _hits_x = new float[_maxNHits];
  _hits_y = new float[_maxNHits];
  _hits_z = new float[_maxNHits];
  _hits_t = new float[_maxNHits];

  _tower_FHCAL_E = new float[_maxNTowers];
  _tower_FHCAL_iEta  = new int[_maxNTowers];
  _tower_FHCAL_iPhi = new int[_maxNTowers];
  _tower_FHCAL_trueID = new int[_maxNTowers];
  _cluster_FHCAL_E = new float[_maxNclusters];
  _cluster_FHCAL_Eta = new float[_maxNclusters];
  _cluster_FHCAL_Phi = new float[_maxNclusters];
  _cluster_FHCAL_NTower = new int[_maxNclusters];
  _cluster_FHCAL_trueID = new int[_maxNclusters];

  _tower_HCALIN_E = new float[_maxNTowers];
  _tower_HCALIN_iEta  = new int[_maxNTowers];
  _tower_HCALIN_iPhi = new int[_maxNTowers];
  _tower_HCALIN_trueID = new int[_maxNTowers];
  _cluster_HCALIN_E = new float[_maxNclusters];
  _cluster_HCALIN_Eta = new float[_maxNclusters];
  _cluster_HCALIN_Phi = new float[_maxNclusters];
  _cluster_HCALIN_NTower = new int[_maxNclusters];
  _cluster_HCALIN_trueID = new int[_maxNclusters];

  _tower_HCALOUT_E = new float[_maxNTowers];
  _tower_HCALOUT_iEta  = new int[_maxNTowers];
  _tower_HCALOUT_iPhi = new int[_maxNTowers];
  _tower_HCALOUT_trueID = new int[_maxNTowers];
  _cluster_HCALOUT_E = new float[_maxNclusters];
  _cluster_HCALOUT_Eta = new float[_maxNclusters];
  _cluster_HCALOUT_Phi = new float[_maxNclusters];
  _cluster_HCALOUT_NTower = new int[_maxNclusters];
  _cluster_HCALOUT_trueID = new int[_maxNclusters];

  _tower_EHCAL_E = new float[_maxNTowers];
  _tower_EHCAL_iEta  = new int[_maxNTowers];
  _tower_EHCAL_iPhi = new int[_maxNTowers];
  _tower_EHCAL_trueID = new int[_maxNTowers];
  _cluster_EHCAL_E = new float[_maxNclusters];
  _cluster_EHCAL_Eta = new float[_maxNclusters];
  _cluster_EHCAL_Phi = new float[_maxNclusters];
  _cluster_EHCAL_NTower = new int[_maxNclusters];
  _cluster_EHCAL_trueID = new int[_maxNclusters];
  
  _tower_DRCALO_E = new float[_maxNTowersDR];
  _tower_DRCALO_NScint  = new int[_maxNTowersDR];
  _tower_DRCALO_NCerenkov  = new int[_maxNTowersDR];
  _tower_DRCALO_iEta  = new int[_maxNTowersDR];
  _tower_DRCALO_iPhi = new int[_maxNTowersDR];
  _tower_DRCALO_trueID = new int[_maxNTowersDR];

  _tower_FEMC_E = new float[_maxNTowers];
  _tower_FEMC_iEta = new int[_maxNTowers];
  _tower_FEMC_iPhi = new int[_maxNTowers];
  _tower_FEMC_trueID = new int[_maxNTowers];
  _cluster_FEMC_E = new float[_maxNclusters];
  _cluster_FEMC_Eta = new float[_maxNclusters];
  _cluster_FEMC_Phi = new float[_maxNclusters];
  _cluster_FEMC_NTower = new int[_maxNclusters];
  _cluster_FEMC_trueID = new int[_maxNclusters];

  _tower_CEMC_E = new float[_maxNTowers];
  _tower_CEMC_iEta = new int[_maxNTowers];
  _tower_CEMC_iPhi = new int[_maxNTowers];
  _tower_CEMC_trueID = new int[_maxNTowers];
  _cluster_CEMC_E = new float[_maxNclusters];
  _cluster_CEMC_Eta = new float[_maxNclusters];
  _cluster_CEMC_Phi = new float[_maxNclusters];
  _cluster_CEMC_NTower = new int[_maxNclusters];
  _cluster_CEMC_trueID = new int[_maxNclusters];

  _tower_EEMC_E = new float[_maxNTowers];
  _tower_EEMC_iEta = new int[_maxNTowers];
  _tower_EEMC_iPhi = new int[_maxNTowers];
  _tower_EEMC_trueID = new int[_maxNTowers];
  _cluster_EEMC_E = new float[_maxNclusters];
  _cluster_EEMC_Eta = new float[_maxNclusters];
  _cluster_EEMC_Phi = new float[_maxNclusters];
  _cluster_EEMC_NTower = new int[_maxNclusters];
  _cluster_EEMC_trueID = new int[_maxNclusters];
  
  _track_ID = new float[_maxNTracks];
  _track_trueID = new float[_maxNTracks];
  _track_px = new float[_maxNTracks];
  _track_py= new float[_maxNTracks];
  _track_pz = new float[_maxNTracks];
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

  _mcpart_ID = new float[_maxNMCPart];
  _mcpart_ID_parent = new float[_maxNMCPart];
  _mcpart_PDG = new float[_maxNMCPart];
  _mcpart_E = new float[_maxNMCPart];
  _mcpart_px = new float[_maxNMCPart];
  _mcpart_py = new float[_maxNMCPart];
  _mcpart_pz = new float[_maxNMCPart];
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
    _event_tree->Branch("hits_trueID", _hits_trueID, "hits_trueID[nHits]/I");
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
    // towers FHCAL
    _event_tree->Branch("tower_FHCAL_N", &_nTowers_FHCAL,"tower_FHCAL_N/I");
    _event_tree->Branch("tower_FHCAL_E", _tower_FHCAL_E, "tower_FHCAL_E[tower_FHCAL_N]/F");
    _event_tree->Branch("tower_FHCAL_iEta",_tower_FHCAL_iEta, "tower_FHCAL_iEta[tower_FHCAL_N]/I");
    _event_tree->Branch("tower_FHCAL_iPhi",_tower_FHCAL_iPhi, "tower_FHCAL_iPhi[tower_FHCAL_N]/I");
    _event_tree->Branch("tower_FHCAL_trueID", _tower_FHCAL_trueID, "tower_FHCAL_trueID[tower_FHCAL_N]/I");
    if(_do_CLUSTERS){
      // clusters FHCAL
      _event_tree->Branch("cluster_FHCAL_N", &_nclusters_FHCAL,"cluster_FHCAL_N/I");
      _event_tree->Branch("cluster_FHCAL_E", _cluster_FHCAL_E, "cluster_FHCAL_E[cluster_FHCAL_N]/F");
      _event_tree->Branch("cluster_FHCAL_Eta",_cluster_FHCAL_Eta, "cluster_FHCAL_Eta[cluster_FHCAL_N]/F");
      _event_tree->Branch("cluster_FHCAL_Phi",_cluster_FHCAL_Phi, "cluster_FHCAL_Phi[cluster_FHCAL_N]/F");
      _event_tree->Branch("cluster_FHCAL_NTower",_cluster_FHCAL_NTower, "cluster_FHCAL_NTower[cluster_FHCAL_N]/I");
      _event_tree->Branch("cluster_FHCAL_trueID", _cluster_FHCAL_trueID, "cluster_FHCAL_trueID[cluster_FHCAL_N]/I");
    }
  }
  if(_do_HCALIN){
    // towers HCAL-in
    _event_tree->Branch("tower_HCALIN_N", &_nTowers_HCALIN,"tower_HCALIN_N/I");
    _event_tree->Branch("tower_HCALIN_E", _tower_HCALIN_E, "tower_HCALIN_E[tower_HCALIN_N]/F");
    _event_tree->Branch("tower_HCALIN_iEta",_tower_HCALIN_iEta, "tower_HCALIN_iEta[tower_HCALIN_N]/I");
    _event_tree->Branch("tower_HCALIN_iPhi",_tower_HCALIN_iPhi, "tower_HCALIN_iPhi[tower_HCALIN_N]/I");
    _event_tree->Branch("tower_HCALIN_trueID", _tower_HCALIN_trueID, "tower_HCALIN_trueID[tower_HCALIN_N]/I");
    if(_do_CLUSTERS){
      // clusters HCAL-in
      _event_tree->Branch("cluster_HCALIN_N", &_nclusters_HCALIN,"cluster_HCALIN_N/I");
      _event_tree->Branch("cluster_HCALIN_E", _cluster_HCALIN_E, "cluster_HCALIN_E[cluster_HCALIN_N]/F");
      _event_tree->Branch("cluster_HCALIN_Eta",_cluster_HCALIN_Eta, "cluster_HCALIN_Eta[cluster_HCALIN_N]/F");
      _event_tree->Branch("cluster_HCALIN_Phi",_cluster_HCALIN_Phi, "cluster_HCALIN_Phi[cluster_HCALIN_N]/F");
      _event_tree->Branch("cluster_HCALIN_NTower",_cluster_HCALIN_NTower, "cluster_HCALIN_NTower[cluster_HCALIN_N]/I");
      _event_tree->Branch("cluster_HCALIN_trueID", _cluster_HCALIN_trueID, "cluster_HCALIN_trueID[cluster_HCALIN_N]/I");
    }
  }
  if(_do_HCALOUT){
    // towers HCAL-out
    _event_tree->Branch("tower_HCALOUT_N", &_nTowers_HCALOUT,"tower_HCALOUT_N/I");
    _event_tree->Branch("tower_HCALOUT_E", _tower_HCALOUT_E, "tower_HCALOUT_E[tower_HCALOUT_N]/F");
    _event_tree->Branch("tower_HCALOUT_iEta",_tower_HCALOUT_iEta, "tower_HCALOUT_iEta[tower_HCALOUT_N]/I");
    _event_tree->Branch("tower_HCALOUT_iPhi",_tower_HCALOUT_iPhi, "tower_HCALOUT_iPhi[tower_HCALOUT_N]/I");
    _event_tree->Branch("tower_HCALOUT_trueID", _tower_HCALOUT_trueID, "tower_HCALOUT_trueID[tower_HCALOUT_N]/I");
    if(_do_CLUSTERS){
      // clusters HCAL-out
      _event_tree->Branch("cluster_HCALOUT_N", &_nclusters_HCALOUT,"cluster_HCALOUT_N/I");
      _event_tree->Branch("cluster_HCALOUT_E", _cluster_HCALOUT_E, "cluster_HCALOUT_E[cluster_HCALOUT_N]/F");
      _event_tree->Branch("cluster_HCALOUT_Eta",_cluster_HCALOUT_Eta, "cluster_HCALOUT_Eta[cluster_HCALOUT_N]/F");
      _event_tree->Branch("cluster_HCALOUT_Phi",_cluster_HCALOUT_Phi, "cluster_HCALOUT_Phi[cluster_HCALOUT_N]/F");
      _event_tree->Branch("cluster_HCALOUT_NTower",_cluster_HCALOUT_NTower, "cluster_HCALOUT_NTower[cluster_HCALOUT_N]/I");
      _event_tree->Branch("cluster_HCALOUT_trueID", _cluster_HCALOUT_trueID, "cluster_HCALOUT_trueID[cluster_HCALOUT_N]/I");
    }
  }
  if(_do_EHCAL){
    // towers EHCAL
    _event_tree->Branch("tower_EHCAL_N", &_nTowers_EHCAL,"tower_EHCAL_N/I");
    _event_tree->Branch("tower_EHCAL_E", _tower_EHCAL_E, "tower_EHCAL_E[tower_EHCAL_N]/F");
    _event_tree->Branch("tower_EHCAL_iEta",_tower_EHCAL_iEta, "tower_EHCAL_iEta[tower_EHCAL_N]/I");
    _event_tree->Branch("tower_EHCAL_iPhi",_tower_EHCAL_iPhi, "tower_EHCAL_iPhi[tower_EHCAL_N]/I");
    _event_tree->Branch("tower_EHCAL_trueID", _tower_EHCAL_trueID, "tower_EHCAL_trueID[tower_EHCAL_N]/I");
    if(_do_CLUSTERS){
      // clusters EHCAL
      _event_tree->Branch("cluster_EHCAL_N", &_nclusters_EHCAL,"cluster_EHCAL_N/I");
      _event_tree->Branch("cluster_EHCAL_E", _cluster_EHCAL_E, "cluster_EHCAL_E[cluster_EHCAL_N]/F");
      _event_tree->Branch("cluster_EHCAL_Eta",_cluster_EHCAL_Eta, "cluster_EHCAL_Eta[cluster_EHCAL_N]/F");
      _event_tree->Branch("cluster_EHCAL_Phi",_cluster_EHCAL_Phi, "cluster_EHCAL_Phi[cluster_EHCAL_N]/F");
      _event_tree->Branch("cluster_EHCAL_NTower",_cluster_EHCAL_NTower, "cluster_EHCAL_NTower[cluster_EHCAL_N]/I");
      _event_tree->Branch("cluster_EHCAL_trueID", _cluster_EHCAL_trueID, "cluster_EHCAL_trueID[cluster_EHCAL_N]/I");
    }
  }
  if(_do_DRCALO){
    // towers DRCALO
    _event_tree->Branch("tower_DRCALO_N", &_nTowers_DRCALO,"tower_DRCALO_N/I");
    _event_tree->Branch("tower_DRCALO_E", _tower_DRCALO_E, "tower_DRCALO_E[tower_DRCALO_N]/F");
    _event_tree->Branch("tower_DRCALO_NScint",_tower_DRCALO_NScint, "tower_DRCALO_NScint[tower_DRCALO_N]/I");
    _event_tree->Branch("tower_DRCALO_NCerenkov",_tower_DRCALO_NCerenkov, "tower_DRCALO_NCerenkov[tower_DRCALO_N]/I");
    _event_tree->Branch("tower_DRCALO_iEta",_tower_DRCALO_iEta, "tower_DRCALO_iEta[tower_DRCALO_N]/I");
    _event_tree->Branch("tower_DRCALO_iPhi",_tower_DRCALO_iPhi, "tower_DRCALO_iPhi[tower_DRCALO_N]/I");
    _event_tree->Branch("tower_DRCALO_trueID", _tower_DRCALO_trueID, "tower_DRCALO_trueID[tower_DRCALO_N]/I");
  }
  if(_do_FEMC){
    // towers FEMC
    _event_tree->Branch("tower_FEMC_N",  &_nTowers_FEMC,"tower_FEMC_N/I");
    _event_tree->Branch("tower_FEMC_E",  _tower_FEMC_E, "tower_FEMC_E[tower_FEMC_N]/F");
    _event_tree->Branch("tower_FEMC_iEta", _tower_FEMC_iEta, "tower_FEMC_iEta[tower_FEMC_N]/I");
    _event_tree->Branch("tower_FEMC_iPhi", _tower_FEMC_iPhi, "tower_FEMC_iPhi[tower_FEMC_N]/I");
    _event_tree->Branch("tower_FEMC_trueID",  _tower_FEMC_trueID, "tower_FEMC_trueID[tower_FEMC_N]/I");
    if(_do_CLUSTERS){
      // clusters FEMC
      _event_tree->Branch("cluster_FEMC_N",  &_nclusters_FEMC,"cluster_FEMC_N/I");
      _event_tree->Branch("cluster_FEMC_E",  _cluster_FEMC_E, "cluster_FEMC_E[cluster_FEMC_N]/F");
      _event_tree->Branch("cluster_FEMC_Eta", _cluster_FEMC_Eta, "cluster_FEMC_Eta[cluster_FEMC_N]/F");
      _event_tree->Branch("cluster_FEMC_Phi", _cluster_FEMC_Phi, "cluster_FEMC_Phi[cluster_FEMC_N]/F");
      _event_tree->Branch("cluster_FEMC_NTower", _cluster_FEMC_NTower, "cluster_FEMC_NTower[cluster_FEMC_N]/I");
      _event_tree->Branch("cluster_FEMC_trueID",  _cluster_FEMC_trueID, "cluster_FEMC_trueID[cluster_FEMC_N]/I");
    }
  }
  if(_do_CEMC){
    // towers CEMC
    _event_tree->Branch("tower_CEMC_N",  &_nTowers_CEMC,"tower_CEMC_N/I");
    _event_tree->Branch("tower_CEMC_E",  _tower_CEMC_E, "tower_CEMC_E[tower_CEMC_N]/F");
    _event_tree->Branch("tower_CEMC_iEta", _tower_CEMC_iEta, "tower_CEMC_iEta[tower_CEMC_N]/I");
    _event_tree->Branch("tower_CEMC_iPhi", _tower_CEMC_iPhi, "tower_CEMC_iPhi[tower_CEMC_N]/I");
    _event_tree->Branch("tower_CEMC_trueID",  _tower_CEMC_trueID, "tower_CEMC_trueID[tower_CEMC_N]/I");
    if(_do_CLUSTERS){
      // clusters CEMC
      _event_tree->Branch("cluster_CEMC_N",  &_nclusters_CEMC,"cluster_CEMC_N/I");
      _event_tree->Branch("cluster_CEMC_E",  _cluster_CEMC_E, "cluster_CEMC_E[cluster_CEMC_N]/F");
      _event_tree->Branch("cluster_CEMC_Eta", _cluster_CEMC_Eta, "cluster_CEMC_Eta[cluster_CEMC_N]/F");
      _event_tree->Branch("cluster_CEMC_Phi", _cluster_CEMC_Phi, "cluster_CEMC_Phi[cluster_CEMC_N]/F");
      _event_tree->Branch("cluster_CEMC_NTower", _cluster_CEMC_NTower, "cluster_CEMC_NTower[cluster_CEMC_N]/I");
      _event_tree->Branch("cluster_CEMC_trueID",  _cluster_CEMC_trueID, "cluster_CEMC_trueID[cluster_CEMC_N]/I");
    }
  }
  if(_do_EEMC){
    // towers EEMC
    _event_tree->Branch("tower_EEMC_N",  &_nTowers_EEMC,"tower_EEMC_N/I");
    _event_tree->Branch("tower_EEMC_E",  _tower_EEMC_E, "tower_EEMC_E[tower_EEMC_N]/F");
    _event_tree->Branch("tower_EEMC_iEta", _tower_EEMC_iEta, "tower_EEMC_iEta[tower_EEMC_N]/I");
    _event_tree->Branch("tower_EEMC_iPhi", _tower_EEMC_iPhi, "tower_EEMC_iPhi[tower_EEMC_N]/I");
    _event_tree->Branch("tower_EEMC_trueID",  _tower_EEMC_trueID, "tower_EEMC_trueID[tower_EEMC_N]/I");
    if(_do_CLUSTERS){
      // clusters EEMC
      _event_tree->Branch("cluster_EEMC_N",  &_nclusters_EEMC,"cluster_EEMC_N/I");
      _event_tree->Branch("cluster_EEMC_E",  _cluster_EEMC_E, "cluster_EEMC_E[cluster_EEMC_N]/F");
      _event_tree->Branch("cluster_EEMC_Eta", _cluster_EEMC_Eta, "cluster_EEMC_Eta[cluster_EEMC_N]/F");
      _event_tree->Branch("cluster_EEMC_Phi", _cluster_EEMC_Phi, "cluster_EEMC_Phi[cluster_EEMC_N]/F");
      _event_tree->Branch("cluster_EEMC_NTower", _cluster_EEMC_NTower, "cluster_EEMC_NTower[cluster_EEMC_N]/I");
      _event_tree->Branch("cluster_EEMC_trueID",  _cluster_EEMC_trueID, "cluster_EEMC_trueID[cluster_EEMC_N]/I");
    }
  }
  if(_do_VERTEX){
    // vertex
    _event_tree->Branch("vertex_x", &_vertex_x, "vertex_x/F");
    _event_tree->Branch("vertex_y", &_vertex_y, "vertex_y/F");
    _event_tree->Branch("vertex_z", &_vertex_z, "vertex_z/F");
    _event_tree->Branch("vertex_NCont", &_vertex_NCont, "vertex_NCont/I");
    _event_tree->Branch("vertex_true_x", &_vertex_true_x, "vertex_true_x/F");
    _event_tree->Branch("vertex_true_y", &_vertex_true_y, "vertex_true_y/F");
    _event_tree->Branch("vertex_true_z", &_vertex_true_z, "vertex_true_z/F");

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
  // if(Verbosity() > 0){
    cout << "entered process_event" << endl;
    // }
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
  if(_do_HCALIN){
    if (!_caloevalstackHCALIN)
    {
      _caloevalstackHCALIN = new CaloEvalStack(topNode, "HCALIN");
      _caloevalstackHCALIN->set_strict(_strict);
      _caloevalstackHCALIN->set_verbosity(Verbosity() + 1);
    }
    else
    {
      _caloevalstackHCALIN->next_event(topNode);
    }
  }
  if(_do_HCALOUT){
    if (!_caloevalstackHCALOUT)
    {
      _caloevalstackHCALOUT = new CaloEvalStack(topNode, "HCALOUT");
      _caloevalstackHCALOUT->set_strict(_strict);
      _caloevalstackHCALOUT->set_verbosity(Verbosity() + 1);
    }
    else
    {
      _caloevalstackHCALOUT->next_event(topNode);
    }
  }
  if(_do_EHCAL){
    if (!_caloevalstackEHCAL)
    {
      _caloevalstackEHCAL = new CaloEvalStack(topNode, "EHCAL");
      _caloevalstackEHCAL->set_strict(_strict);
      _caloevalstackEHCAL->set_verbosity(Verbosity() + 1);
    }
    else
    {
      _caloevalstackEHCAL->next_event(topNode);
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
  if(_do_CEMC){
    if (!_caloevalstackCEMC)
    {
      _caloevalstackCEMC = new CaloEvalStack(topNode, "CEMC");
      _caloevalstackCEMC->set_strict(_strict);
      _caloevalstackCEMC->set_verbosity(Verbosity() + 1);
    }
    else
    {
      _caloevalstackCEMC->next_event(topNode);
    }
  }
  if(_do_EEMC){
    if (!_caloevalstackEEMC)
    {
      _caloevalstackEEMC = new CaloEvalStack(topNode, "EEMC");
      _caloevalstackEEMC->set_strict(_strict);
      _caloevalstackEEMC->set_verbosity(Verbosity() + 1);
    }
    else
    {
      _caloevalstackEEMC->next_event(topNode);
    }
  }
  if (Verbosity() > 0) {cout << "loaded evalstack" << endl;}


  // fill the Evaluator Tree
  fillOutputNtuples(topNode);


  ++_ievent;

  return Fun4AllReturnCodes::EVENT_OK;
}

void EventEvaluator::fillOutputNtuples(PHCompositeNode* topNode)
{
  if (Verbosity() > 2){ cout << "EventEvaluator::fillOutputNtuples() entered" << endl;}


  //----------------------
  // fill the Event Tree
  //----------------------

  //----------------------
  //    VERTEX
  //----------------------
  SvtxVertexMap* vertexmap = findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMap");
  if(_do_VERTEX){
    if (vertexmap)
    {
      if (!vertexmap->empty())
      {
        if (Verbosity() > 0){ cout << "saving vertex" << endl;}
        SvtxVertex* vertex = (vertexmap->begin())->second;

        _vertex_x = vertex->get_x();
        _vertex_y = vertex->get_y();
        _vertex_z = vertex->get_z();
        _vertex_NCont = vertex->size_tracks();
      }
    }
  }
  //----------------------
  //    HITS
  //----------------------
  if(_do_HITS){
    if (Verbosity() > 0){cout << "saving hits" << endl;}
    _nHitsLayers = 0;
    PHG4TruthInfoContainer* truthinfocontainerHits = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
    for (int iIndex = 0; iIndex < 60; ++iIndex)
    {
      if (GetProjectionNameFromIndex(iIndex).find("NOTHING")!= std::string::npos ) continue;
      if (GetProjectionNameFromIndex(iIndex).find("FHCAL")!= std::string::npos ) continue;
      if (GetProjectionNameFromIndex(iIndex).find("EHCAL")!= std::string::npos ) continue;
      if (GetProjectionNameFromIndex(iIndex).find("FEMC")!= std::string::npos ) continue;
      if (GetProjectionNameFromIndex(iIndex).find("CEMC")!= std::string::npos ) continue;
      if (GetProjectionNameFromIndex(iIndex).find("EEMC")!= std::string::npos ) continue;
      if (GetProjectionNameFromIndex(iIndex).find("HCALIN")!= std::string::npos ) continue;
      if (GetProjectionNameFromIndex(iIndex).find("HCALOUT")!= std::string::npos ) continue;
      string nodename = "G4HIT_" + GetProjectionNameFromIndex(iIndex);
      PHG4HitContainer* hits = findNode::getClass<PHG4HitContainer>(topNode, nodename);
      if (hits)
      {
        if (Verbosity() > 1){ cout << __PRETTY_FUNCTION__ << " number of hits: " << hits->size() << endl;}
        PHG4HitContainer::ConstRange hit_range = hits->getHits();
        for (PHG4HitContainer::ConstIterator hit_iter = hit_range.first; hit_iter != hit_range.second; hit_iter++)
        {
          // if(Verbosity() > 0) cout << __PRETTY_FUNCTION__ << " checking hit id " << hit_iter->second->get_trkid() << " against " << track->get_truth_track_id() << endl;
          // if (hit_iter->second->get_trkid() - track->get_truth_track_id() == 0)
          // {
          if (Verbosity() > 1){ cout << __PRETTY_FUNCTION__ << " found hit with id " << hit_iter->second->get_trkid() << endl;}
          _hits_x[_nHitsLayers] = hit_iter->second->get_x(0);
          _hits_y[_nHitsLayers] = hit_iter->second->get_y(0);
          _hits_z[_nHitsLayers] = hit_iter->second->get_z(0);
          _hits_t[_nHitsLayers] = hit_iter->second->get_t(0);
          _hits_layerID[_nHitsLayers] = iIndex;
          if (truthinfocontainerHits)
          {
            PHG4Particle *particle = truthinfocontainerHits->GetParticle( hit_iter->second->get_trkid() );

            if(particle->get_parent_id()!=0){
              PHG4Particle* g4particleMother = truthinfocontainerHits->GetParticle( hit_iter->second->get_trkid() );
              int mcSteps = 0;
              while(g4particleMother->get_parent_id()!=0){
                g4particleMother = truthinfocontainerHits->GetParticle(g4particleMother->get_parent_id());
                mcSteps+=1;
              }
              if(mcSteps<=_depth_MCstack){
                _hits_trueID[_nHitsLayers] = hit_iter->second->get_trkid();
              } else {
                PHG4Particle* g4particleMother2 = truthinfocontainerHits->GetParticle( hit_iter->second->get_trkid() );
                int mcSteps2 = 0;
                while(g4particleMother2->get_parent_id()!=0 && (mcSteps2<(mcSteps-_depth_MCstack+1))){
                  g4particleMother2 = truthinfocontainerHits->GetParticle(g4particleMother2->get_parent_id());
                  mcSteps2+=1;
                }
                _hits_trueID[_nHitsLayers] = g4particleMother2->get_parent_id();
              }
            } else {
              _hits_trueID[_nHitsLayers] = hit_iter->second->get_trkid();
            }
          }
          _nHitsLayers++;

          // }
        }
        if (Verbosity() > 0){ cout << "saved\t" << _nHitsLayers << "\thits for " << GetProjectionNameFromIndex(iIndex) << endl;}
      }
      else
      {
        if (Verbosity() > 0){ cout << __PRETTY_FUNCTION__ << " could not find " << nodename << endl;}
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
      if (Verbosity() > 0){ cout << "saving HCAL towers" << endl;}
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
        if (Verbosity() > 0){ cout << PHWHERE << " ERROR: Can't find " << towergeomnodeFHCAL << endl;}
        // return;
      }
      if (Verbosity() > 0){ cout << "saved\t" << _nTowers_FHCAL << "\tFHCAL towers" << endl;}
    }
    else
    {
      if (Verbosity() > 0){ cout << PHWHERE << " ERROR: Can't find " << towernodeFHCAL << endl;}
      // return;
    }
  }
  //----------------------
  //    TOWERS HCALIN
  //----------------------
  if(_do_HCALIN){
    CaloRawTowerEval* towerevalHCALIN = _caloevalstackHCALIN->get_rawtower_eval();
    _nTowers_HCALIN = 0;
    string towernodeHCALIN = "TOWER_CALIB_HCALIN";
    RawTowerContainer* towersHCALIN = findNode::getClass<RawTowerContainer>(topNode, towernodeHCALIN.c_str());
    if (towersHCALIN)
    {
      if (Verbosity() > 0){ cout << "saving HCAL towers" << endl;}
      string towergeomnodeHCALIN = "TOWERGEOM_HCALIN";
      RawTowerGeomContainer* towergeomHCALIN = findNode::getClass<RawTowerGeomContainer>(topNode, towergeomnodeHCALIN.c_str());
      if (towergeomHCALIN)
      {
        RawTowerContainer::ConstRange begin_end = towersHCALIN->getTowers();
        RawTowerContainer::ConstIterator rtiter;
        for (rtiter = begin_end.first; rtiter != begin_end.second; ++rtiter)
        {
          RawTower* tower = rtiter->second;
          if (tower)
          {
            // min energy cut
            if (tower->get_energy() < _reco_e_threshold) continue;
            _tower_HCALIN_iEta[_nTowers_HCALIN] = tower->get_bineta();
            _tower_HCALIN_iPhi[_nTowers_HCALIN] = tower->get_binphi();
            _tower_HCALIN_E[_nTowers_HCALIN] = tower->get_energy();

            PHG4Particle* primary = towerevalHCALIN->max_truth_primary_particle_by_energy(tower);
            if (primary)
            {
              _tower_HCALIN_trueID[_nTowers_HCALIN] = primary->get_track_id();
              // gflavor = primary->get_pid();
              // efromtruth = towerevalHCALIN->get_energy_contribution(tower, primary);
            } else {
              _tower_HCALIN_trueID[_nTowers_HCALIN] = -10;
            }
            _nTowers_HCALIN++;
          }
        }
      }
      else
      {
        if (Verbosity() > 0){ cout << PHWHERE << " ERROR: Can't find " << towergeomnodeHCALIN << endl;}
        // return;
      }
      if (Verbosity() > 0){ cout << "saved\t" << _nTowers_HCALIN << "\tHCALIN towers" << endl;}
    }
    else
    {
      if (Verbosity() > 0){ cout << PHWHERE << " ERROR: Can't find " << towernodeHCALIN << endl;}
      // return;
    }
  }
  //----------------------
  //    TOWERS HCALOUT
  //----------------------
  if(_do_HCALOUT){
    CaloRawTowerEval* towerevalHCALOUT = _caloevalstackHCALOUT->get_rawtower_eval();
    _nTowers_HCALOUT = 0;
    string towernodeHCALOUT = "TOWER_CALIB_HCALOUT";
    RawTowerContainer* towersHCALOUT = findNode::getClass<RawTowerContainer>(topNode, towernodeHCALOUT.c_str());
    if (towersHCALOUT)
    {
      if (Verbosity() > 0){ cout << "saving HCAL towers" << endl;}
      string towergeomnodeHCALOUT = "TOWERGEOM_HCALOUT";
      RawTowerGeomContainer* towergeomHCALOUT = findNode::getClass<RawTowerGeomContainer>(topNode, towergeomnodeHCALOUT.c_str());
      if (towergeomHCALOUT)
      {
        RawTowerContainer::ConstRange begin_end = towersHCALOUT->getTowers();
        RawTowerContainer::ConstIterator rtiter;
        for (rtiter = begin_end.first; rtiter != begin_end.second; ++rtiter)
        {
          RawTower* tower = rtiter->second;
          if (tower)
          {
            // min energy cut
            if (tower->get_energy() < _reco_e_threshold) continue;
            _tower_HCALOUT_iEta[_nTowers_HCALOUT] = tower->get_bineta();
            _tower_HCALOUT_iPhi[_nTowers_HCALOUT] = tower->get_binphi();
            _tower_HCALOUT_E[_nTowers_HCALOUT] = tower->get_energy();

            PHG4Particle* primary = towerevalHCALOUT->max_truth_primary_particle_by_energy(tower);
            if (primary)
            {
              _tower_HCALOUT_trueID[_nTowers_HCALOUT] = primary->get_track_id();
              // gflavor = primary->get_pid();
              // efromtruth = towerevalHCALOUT->get_energy_contribution(tower, primary);
            } else {
              _tower_HCALOUT_trueID[_nTowers_HCALOUT] = -10;
            }
            _nTowers_HCALOUT++;
          }
        }
      }
      else
      {
        if (Verbosity() > 0){ cout << PHWHERE << " ERROR: Can't find " << towergeomnodeHCALOUT << endl;}
        // return;
      }
      if (Verbosity() > 0){ cout << "saved\t" << _nTowers_HCALOUT << "\tHCALOUT towers" << endl;}
    }
    else
    {
      if (Verbosity() > 0){ cout << PHWHERE << " ERROR: Can't find " << towernodeHCALOUT << endl;}
      // return;
    }
  }
  //----------------------
  //    TOWERS EHCAL
  //----------------------
  if(_do_EHCAL){
    CaloRawTowerEval* towerevalEHCAL = _caloevalstackEHCAL->get_rawtower_eval();
    _nTowers_EHCAL = 0;
    string towernodeEHCAL = "TOWER_CALIB_EHCAL";
    RawTowerContainer* towersEHCAL = findNode::getClass<RawTowerContainer>(topNode, towernodeEHCAL.c_str());
    if (towersEHCAL)
    {
      if (Verbosity() > 0){ cout << "saving HCAL towers" << endl;}
      string towergeomnodeEHCAL = "TOWERGEOM_EHCAL";
      RawTowerGeomContainer* towergeomEHCAL = findNode::getClass<RawTowerGeomContainer>(topNode, towergeomnodeEHCAL.c_str());
      if (towergeomEHCAL)
      {
        RawTowerContainer::ConstRange begin_end = towersEHCAL->getTowers();
        RawTowerContainer::ConstIterator rtiter;
        for (rtiter = begin_end.first; rtiter != begin_end.second; ++rtiter)
        {
          RawTower* tower = rtiter->second;
          if (tower)
          {
            // min energy cut
            if (tower->get_energy() < _reco_e_threshold) continue;
            _tower_EHCAL_iEta[_nTowers_EHCAL] = tower->get_bineta();
            _tower_EHCAL_iPhi[_nTowers_EHCAL] = tower->get_binphi();
            _tower_EHCAL_E[_nTowers_EHCAL] = tower->get_energy();

            PHG4Particle* primary = towerevalEHCAL->max_truth_primary_particle_by_energy(tower);
            if (primary)
            {
              _tower_EHCAL_trueID[_nTowers_EHCAL] = primary->get_track_id();
              // gflavor = primary->get_pid();
              // efromtruth = towerevalEHCAL->get_energy_contribution(tower, primary);
            } else {
              _tower_EHCAL_trueID[_nTowers_EHCAL] = -10;
            }
            _nTowers_EHCAL++;
          }
        }
      }
      else
      {
        if (Verbosity() > 0){ cout << PHWHERE << " ERROR: Can't find " << towergeomnodeEHCAL << endl;}
        // return;
      }
      if (Verbosity() > 0){ cout << "saved\t" << _nTowers_EHCAL << "\tEHCAL towers" << endl;}
    }
    else
    {
      if (Verbosity() > 0){ cout << PHWHERE << " ERROR: Can't find " << towernodeEHCAL << endl;}
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
      if (Verbosity() > 0){ cout << "saving DRCALO towers" << endl;}
      string towergeomnodeDRCALO = "TOWERGEOM_DRCALO";
      RawTowerGeomContainer* towergeomDRCALO = findNode::getClass<RawTowerGeomContainer>(topNode, towergeomnodeDRCALO.c_str());
      if (towergeomDRCALO)
      {
      if (Verbosity() > 0){ cout << "found DRCALO geom" << endl;}
        RawTowerContainer::ConstRange begin_end = towersDRCALO->getTowers();
        RawTowerContainer::ConstIterator rtiter;
        for (rtiter = begin_end.first; rtiter != begin_end.second; ++rtiter)
        {
          RawTower* tower = rtiter->second;
          if (tower)
          {
            // min energy cut
            // if (tower->get_energy() < _reco_e_threshold) continue;

            _tower_DRCALO_iEta[_nTowers_DRCALO] = tower->get_bineta();
            _tower_DRCALO_iPhi[_nTowers_DRCALO] = tower->get_binphi();
            _tower_DRCALO_E[_nTowers_DRCALO] = tower->get_energy();
            _tower_DRCALO_NScint[_nTowers_DRCALO] = tower->get_scint_gammas();
            _tower_DRCALO_NCerenkov[_nTowers_DRCALO] = tower->get_cerenkov_gammas();
            // cout << "sci gammas: " << tower->get_scint_gammas() << "\tcerenk gammas: " << tower->get_cerenkov_gammas() << endl;

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
      if (Verbosity() > 0){ cout << "finished DRCALO twr loop" << endl;}
      }
      else
      {
        if (Verbosity() > 0){ cout << PHWHERE << " ERROR: Can't find " << towergeomnodeDRCALO << endl;}
        // return;
      }
      if (Verbosity() > 0){ cout << "saved\t" << _nTowers_DRCALO << "\tDRCALO towers" << endl;}
    }
    else
    {
      if (Verbosity() > 0){ cout << PHWHERE << " ERROR: Can't find " << towernodeDRCALO << endl;}
      // return;
    }
  }
  //----------------------
  //    TOWERS FEMC
  //----------------------
  if(_do_FEMC){
    CaloRawTowerEval* towerevalFEMC = _caloevalstackFEMC->get_rawtower_eval();
    _nTowers_FEMC = 0;
    string towernodeFEMC = "TOWER_CALIB_FEMC";
    RawTowerContainer* towersFEMC = findNode::getClass<RawTowerContainer>(topNode, towernodeFEMC.c_str());
    if (towersFEMC)
    {
      if (Verbosity() > 0){ cout << "saving EMC towers" << endl;}
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
        if (Verbosity() > 0){ cout << PHWHERE << " ERROR: Can't find " << towergeomnodeFEMC << endl;}
        // return;
      }
      if (Verbosity() > 0){ cout << "saved\t" << _nTowers_FEMC << "\tFEMC towers" << endl;}
    }
    else
    {
      if (Verbosity() > 0){ cout << PHWHERE << " ERROR: Can't find " << towernodeFEMC << endl;}
      // return;
    }
  }
  //----------------------
  //    TOWERS CEMC
  //----------------------
  if(_do_CEMC){
    CaloRawTowerEval* towerevalCEMC = _caloevalstackCEMC->get_rawtower_eval();
    _nTowers_CEMC = 0;
    string towernodeCEMC = "TOWER_CALIB_CEMC";
    RawTowerContainer* towersCEMC = findNode::getClass<RawTowerContainer>(topNode, towernodeCEMC.c_str());
    if (towersCEMC)
    {
      if (Verbosity() > 0){ cout << "saving EMC towers" << endl;}
      string towergeomnodeCEMC = "TOWERGEOM_CEMC";
      RawTowerGeomContainer* towergeom = findNode::getClass<RawTowerGeomContainer>(topNode, towergeomnodeCEMC.c_str());
      if (towergeom)
      {
        RawTowerContainer::ConstRange begin_end = towersCEMC->getTowers();
        RawTowerContainer::ConstIterator rtiter;
        for (rtiter = begin_end.first; rtiter != begin_end.second; ++rtiter)
        {
          RawTower* tower = rtiter->second;
          if (tower)
          {
            // min energy cut
            if (tower->get_energy() < _reco_e_threshold) continue;

            _tower_CEMC_iEta[_nTowers_CEMC] = tower->get_bineta();
            _tower_CEMC_iPhi[_nTowers_CEMC] = tower->get_binphi();
            _tower_CEMC_E[_nTowers_CEMC] = tower->get_energy();

            PHG4Particle* primary = towerevalCEMC->max_truth_primary_particle_by_energy(tower);
            if (primary)
            {
              _tower_CEMC_trueID[_nTowers_CEMC] = primary->get_track_id();
              // gflavor = primary->get_pid();
              // efromtruth = towerevalCEMC->get_energy_contribution(tower, primary);
            } else {
              _tower_CEMC_trueID[_nTowers_CEMC] = -10;
            }
            _nTowers_CEMC++;
          }
        }
      }
      else
      {
        if (Verbosity() > 0){ cout << PHWHERE << " ERROR: Can't find " << towergeomnodeCEMC << endl;}
        // return;
      }
      if (Verbosity() > 0){ cout << "saved\t" << _nTowers_CEMC << "\tCEMC towers" << endl;}
    }
    else
    {
      if (Verbosity() > 0){ cout << PHWHERE << " ERROR: Can't find " << towernodeCEMC << endl;}
      // return;
    }
  }
  //----------------------
  //    TOWERS EEMC
  //----------------------
  if(_do_EEMC){
    CaloRawTowerEval* towerevalEEMC = _caloevalstackEEMC->get_rawtower_eval();
    _nTowers_EEMC = 0;
    string towernodeEEMC = "TOWER_CALIB_EEMC";
    RawTowerContainer* towersEEMC = findNode::getClass<RawTowerContainer>(topNode, towernodeEEMC.c_str());
    if (towersEEMC)
    {
      if (Verbosity() > 0){ cout << "saving EMC towers" << endl;}
      string towergeomnodeEEMC = "TOWERGEOM_EEMC";
      RawTowerGeomContainer* towergeom = findNode::getClass<RawTowerGeomContainer>(topNode, towergeomnodeEEMC.c_str());
      if (towergeom)
      {
        RawTowerContainer::ConstRange begin_end = towersEEMC->getTowers();
        RawTowerContainer::ConstIterator rtiter;
        for (rtiter = begin_end.first; rtiter != begin_end.second; ++rtiter)
        {
          RawTower* tower = rtiter->second;
          if (tower)
          {
            // min energy cut
            if (tower->get_energy() < _reco_e_threshold) continue;

            _tower_EEMC_iEta[_nTowers_EEMC] = tower->get_bineta();
            _tower_EEMC_iPhi[_nTowers_EEMC] = tower->get_binphi();
            _tower_EEMC_E[_nTowers_EEMC] = tower->get_energy();

            PHG4Particle* primary = towerevalEEMC->max_truth_primary_particle_by_energy(tower);
            if (primary)
            {
              _tower_EEMC_trueID[_nTowers_EEMC] = primary->get_track_id();
              // gflavor = primary->get_pid();
              // efromtruth = towerevalEEMC->get_energy_contribution(tower, primary);
            } else {
              _tower_EEMC_trueID[_nTowers_EEMC] = -10;
            }
            _nTowers_EEMC++;
          }
        }
      }
      else
      {
        if (Verbosity() > 0){ cout << PHWHERE << " ERROR: Can't find " << towergeomnodeEEMC << endl;}
        // return;
      }
      if (Verbosity() > 0){ cout << "saved\t" << _nTowers_EEMC << "\tEEMC towers" << endl;}
    }
    else
    {
      if (Verbosity() > 0){ cout << PHWHERE << " ERROR: Can't find " << towernodeEEMC << endl;}
      // return;
    }
  }
  //------------------------
  // CLUSTERS FHCAL
  //------------------------
  if(_do_FHCAL && _do_CLUSTERS){
    CaloRawClusterEval* clusterevalFHCAL = _caloevalstackFHCAL->get_rawcluster_eval();
    _nclusters_FHCAL = 0;
    if (Verbosity() > 1){ cout << "CaloEvaluator::filling gcluster ntuple..." << endl;}

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
            SvtxVertex* vertex = (vertexmap->begin()->second);
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
  // CLUSTERS HCALIN
  //------------------------
  if(_do_HCALIN && _do_CLUSTERS){
    CaloRawClusterEval* clusterevalHCALIN = _caloevalstackHCALIN->get_rawcluster_eval();
    _nclusters_HCALIN = 0;
    if (Verbosity() > 1){ cout << "CaloEvaluator::filling gcluster ntuple..." << endl;}

    string clusternodeHCALIN = "CLUSTER_HCALIN";
    RawClusterContainer* clustersHCALIN = findNode::getClass<RawClusterContainer>(topNode, clusternodeHCALIN.c_str());
    if (clustersHCALIN)
    {
      // for every cluster
      for (const auto& iterator : clustersHCALIN->getClustersMap())
      {
        RawCluster* cluster = iterator.second;

        if (cluster->get_energy() < _reco_e_threshold) continue;

        _cluster_HCALIN_E[_nclusters_HCALIN] = cluster->get_energy();
        _cluster_HCALIN_NTower[_nclusters_HCALIN] = cluster->getNTowers();
        _cluster_HCALIN_Phi[_nclusters_HCALIN] = cluster->get_phi();

        // require vertex for cluster eta calculation
        if (vertexmap)
        {
          if (!vertexmap->empty())
          {
            SvtxVertex* vertex = (vertexmap->begin()->second);
            _cluster_HCALIN_Eta[_nclusters_HCALIN] = RawClusterUtility::GetPseudorapidity(*cluster, CLHEP::Hep3Vector(vertex->get_x(), vertex->get_y(), vertex->get_z()));
          }
          else
            _cluster_HCALIN_Eta[_nclusters_HCALIN] = -10000;
        }
        else
          _cluster_HCALIN_Eta[_nclusters_HCALIN] = -10000;

        PHG4Particle* primary = clusterevalHCALIN->max_truth_primary_particle_by_energy(cluster);

        if (primary){
          _cluster_HCALIN_trueID[_nclusters_HCALIN]       = primary->get_track_id();
        } else {
          _cluster_HCALIN_trueID[_nclusters_HCALIN]       = -10;
        }

        _nclusters_HCALIN++;
      }
    }
    else
    {
      cerr << PHWHERE << " ERROR: Can't find " << clusternodeHCALIN << endl;
      // return;
    }
  }
  //------------------------
  // CLUSTERS HCALOUT
  //------------------------
  if(_do_HCALOUT && _do_CLUSTERS){
    CaloRawClusterEval* clusterevalHCALOUT = _caloevalstackHCALOUT->get_rawcluster_eval();
    _nclusters_HCALOUT = 0;
    if (Verbosity() > 1){ cout << "CaloEvaluator::filling gcluster ntuple..." << endl;}

    string clusternodeHCALOUT = "CLUSTER_HCALOUT";
    RawClusterContainer* clustersHCALOUT = findNode::getClass<RawClusterContainer>(topNode, clusternodeHCALOUT.c_str());
    if (clustersHCALOUT)
    {
      // for every cluster
      for (const auto& iterator : clustersHCALOUT->getClustersMap())
      {
        RawCluster* cluster = iterator.second;

        if (cluster->get_energy() < _reco_e_threshold) continue;

        _cluster_HCALOUT_E[_nclusters_HCALOUT] = cluster->get_energy();
        _cluster_HCALOUT_NTower[_nclusters_HCALOUT] = cluster->getNTowers();
        _cluster_HCALOUT_Phi[_nclusters_HCALOUT] = cluster->get_phi();

        // require vertex for cluster eta calculation
        if (vertexmap)
        {
          if (!vertexmap->empty())
          {
            SvtxVertex* vertex = (vertexmap->begin()->second);
            _cluster_HCALOUT_Eta[_nclusters_HCALOUT] = RawClusterUtility::GetPseudorapidity(*cluster, CLHEP::Hep3Vector(vertex->get_x(), vertex->get_y(), vertex->get_z()));
          }
          else
            _cluster_HCALOUT_Eta[_nclusters_HCALOUT] = -10000;
        }
        else
          _cluster_HCALOUT_Eta[_nclusters_HCALOUT] = -10000;

        PHG4Particle* primary = clusterevalHCALOUT->max_truth_primary_particle_by_energy(cluster);

        if (primary){
          _cluster_HCALOUT_trueID[_nclusters_HCALOUT]       = primary->get_track_id();
        } else {
          _cluster_HCALOUT_trueID[_nclusters_HCALOUT]       = -10;
        }

        _nclusters_HCALOUT++;
      }
    }
    else
    {
      cerr << PHWHERE << " ERROR: Can't find " << clusternodeHCALOUT << endl;
      // return;
    }
  }
  //------------------------
  // CLUSTERS EHCAL
  //------------------------
  if(_do_EHCAL && _do_CLUSTERS){
    CaloRawClusterEval* clusterevalEHCAL = _caloevalstackEHCAL->get_rawcluster_eval();
    _nclusters_EHCAL = 0;
    if (Verbosity() > 1){ cout << "CaloEvaluator::filling gcluster ntuple..." << endl;}

    string clusternodeEHCAL = "CLUSTER_EHCAL";
    RawClusterContainer* clustersEHCAL = findNode::getClass<RawClusterContainer>(topNode, clusternodeEHCAL.c_str());
    if (clustersEHCAL)
    {
      // for every cluster
      for (const auto& iterator : clustersEHCAL->getClustersMap())
      {
        RawCluster* cluster = iterator.second;

        if (cluster->get_energy() < _reco_e_threshold) continue;

        _cluster_EHCAL_E[_nclusters_EHCAL] = cluster->get_energy();
        _cluster_EHCAL_NTower[_nclusters_EHCAL] = cluster->getNTowers();
        _cluster_EHCAL_Phi[_nclusters_EHCAL] = cluster->get_phi();

        // require vertex for cluster eta calculation
        if (vertexmap)
        {
          if (!vertexmap->empty())
          {
            SvtxVertex* vertex = (vertexmap->begin()->second);
            _cluster_EHCAL_Eta[_nclusters_EHCAL] = RawClusterUtility::GetPseudorapidity(*cluster, CLHEP::Hep3Vector(vertex->get_x(), vertex->get_y(), vertex->get_z()));
          }
          else
            _cluster_EHCAL_Eta[_nclusters_EHCAL] = -10000;
        }
        else
          _cluster_EHCAL_Eta[_nclusters_EHCAL] = -10000;

        PHG4Particle* primary = clusterevalEHCAL->max_truth_primary_particle_by_energy(cluster);

        if (primary){
          _cluster_EHCAL_trueID[_nclusters_EHCAL]       = primary->get_track_id();
        } else {
          _cluster_EHCAL_trueID[_nclusters_EHCAL]       = -10;
        }

        _nclusters_EHCAL++;
      }
    }
    else
    {
      cerr << PHWHERE << " ERROR: Can't find " << clusternodeEHCAL << endl;
      // return;
    }
  }

  //------------------------
  // CLUSTERS FEMC
  //------------------------
  if(_do_FEMC && _do_CLUSTERS){
    CaloRawClusterEval* clusterevalFEMC = _caloevalstackFEMC->get_rawcluster_eval();
    _nclusters_FEMC = 0;
    if (Verbosity() > 1){ cout << "CaloEvaluator::filling gcluster ntuple..." << endl;}

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
            SvtxVertex* vertex = (vertexmap->begin()->second);
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
  // CLUSTERS CEMC
  //------------------------
  if(_do_CEMC && _do_CLUSTERS){
    CaloRawClusterEval* clusterevalCEMC = _caloevalstackCEMC->get_rawcluster_eval();
    _nclusters_CEMC = 0;
    if (Verbosity() > 1){ cout << "CaloEvaluator::filling gcluster ntuple..." << endl;}

    string clusternodeCEMC = "CLUSTER_CEMC";
    RawClusterContainer* clustersCEMC = findNode::getClass<RawClusterContainer>(topNode, clusternodeCEMC.c_str());
    if (clustersCEMC)
    {
      // for every cluster
      for (const auto& iterator : clustersCEMC->getClustersMap())
      {
        RawCluster* cluster = iterator.second;

        if (cluster->get_energy() < _reco_e_threshold) continue;

        _cluster_CEMC_E[_nclusters_CEMC] = cluster->get_energy();
        _cluster_CEMC_NTower[_nclusters_CEMC] = cluster->getNTowers();
        _cluster_CEMC_Phi[_nclusters_CEMC] = cluster->get_phi();

        // require vertex for cluster eta calculation
        if (vertexmap)
        {
          if (!vertexmap->empty())
          {
            SvtxVertex* vertex = (vertexmap->begin()->second);
            _cluster_CEMC_Eta[_nclusters_CEMC] = RawClusterUtility::GetPseudorapidity(*cluster, CLHEP::Hep3Vector(vertex->get_x(), vertex->get_y(), vertex->get_z()));
          }
          else
            _cluster_CEMC_Eta[_nclusters_CEMC] = -10000;
        }
        else
          _cluster_CEMC_Eta[_nclusters_CEMC] = -10000;

        PHG4Particle* primary = clusterevalCEMC->max_truth_primary_particle_by_energy(cluster);

        if (primary){
          _cluster_CEMC_trueID[_nclusters_CEMC]       = primary->get_track_id(); 
        } else {
          _cluster_CEMC_trueID[_nclusters_CEMC] = -10;
        }

        _nclusters_CEMC++;
      }
    }
    else
    {
      cerr << PHWHERE << " ERROR: Can't find " << clusternodeCEMC << endl;
      // return;
    }
  }
  //------------------------
  // CLUSTERS EEMC
  //------------------------
  if(_do_EEMC && _do_CLUSTERS){
    CaloRawClusterEval* clusterevalEEMC = _caloevalstackEEMC->get_rawcluster_eval();
    _nclusters_EEMC = 0;
    if (Verbosity() > 1){ cout << "CaloEvaluator::filling gcluster ntuple..." << endl;}

    string clusternodeEEMC = "CLUSTER_EEMC";
    RawClusterContainer* clustersEEMC = findNode::getClass<RawClusterContainer>(topNode, clusternodeEEMC.c_str());
    if (clustersEEMC)
    {
      // for every cluster
      for (const auto& iterator : clustersEEMC->getClustersMap())
      {
        RawCluster* cluster = iterator.second;

        if (cluster->get_energy() < _reco_e_threshold) continue;

        _cluster_EEMC_E[_nclusters_EEMC] = cluster->get_energy();
        _cluster_EEMC_NTower[_nclusters_EEMC] = cluster->getNTowers();
        _cluster_EEMC_Phi[_nclusters_EEMC] = cluster->get_phi();

        // require vertex for cluster eta calculation
        if (vertexmap)
        {
          if (!vertexmap->empty())
          {
            SvtxVertex* vertex = (vertexmap->begin()->second);
            _cluster_EEMC_Eta[_nclusters_EEMC] = RawClusterUtility::GetPseudorapidity(*cluster, CLHEP::Hep3Vector(vertex->get_x(), vertex->get_y(), vertex->get_z()));
          }
          else
            _cluster_EEMC_Eta[_nclusters_EEMC] = -10000;
        }
        else
          _cluster_EEMC_Eta[_nclusters_EEMC] = -10000;

        PHG4Particle* primary = clusterevalEEMC->max_truth_primary_particle_by_energy(cluster);

        if (primary){
          _cluster_EEMC_trueID[_nclusters_EEMC]       = primary->get_track_id(); 
        } else {
          _cluster_EEMC_trueID[_nclusters_EEMC] = -10;
        }

        _nclusters_EEMC++;
      }
    }
    else
    {
      cerr << PHWHERE << " ERROR: Can't find " << clusternodeEEMC << endl;
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
      if (Verbosity() > 0){ cout << "saving tracks" << endl;}
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
              if (Verbosity() > 1){ cout << __PRETTY_FUNCTION__ << " processing " << trkstates->second->get_name() << endl;}
              string trackStateName = trkstates->second->get_name();
              if (Verbosity() > 1){ cout << __PRETTY_FUNCTION__ << " found " << trkstates->second->get_name() << endl;}
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
                  if (Verbosity() > 1){ cout << __PRETTY_FUNCTION__ << " number of hits: " << hits->size() << endl;}
                  PHG4HitContainer::ConstRange hit_range = hits->getHits();
                  for (PHG4HitContainer::ConstIterator hit_iter = hit_range.first; hit_iter != hit_range.second; hit_iter++)
                  {
                    if (Verbosity() > 1){ cout << __PRETTY_FUNCTION__ << " checking hit id " << hit_iter->second->get_trkid() << " against " << track->get_truth_track_id() << endl;}
                    if (hit_iter->second->get_trkid() - track->get_truth_track_id() == 0)
                    {
                      if (Verbosity() > 1){ cout << __PRETTY_FUNCTION__ << " found hit with id " << hit_iter->second->get_trkid() << endl;}
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
                  if (Verbosity() > 1){ cout << __PRETTY_FUNCTION__ << " could not find " << nodename << endl;}
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
      if (Verbosity() > 0){ cout << "saved\t" << _nTracks << "\ttracks" << endl;}
    }
    else
    {
      if (Verbosity() > 0){ cout << PHWHERE << "SvtxTrackMap node with name TrackMap not found on node tree" << endl;}
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
      if (Verbosity() > 0){ cout << "saving MC particles" << endl;}
      //GetParticleRange for all particles
      //GetPrimaryParticleRange for primary particles
      PHG4TruthInfoContainer::ConstRange range = truthinfocontainer->GetParticleRange();
      for (PHG4TruthInfoContainer::ConstIterator truth_itr = range.first; truth_itr != range.second; ++truth_itr)
      {
        PHG4Particle* g4particle = truth_itr->second;
        if (!g4particle) continue;

        int mcSteps = 0;
        PHG4Particle* g4particleMother = truth_itr->second;
        if(g4particle->get_parent_id()!=0){
          while(g4particleMother->get_parent_id()!=0){
            g4particleMother = truthinfocontainer->GetParticle(g4particleMother->get_parent_id());
            mcSteps+=1;
          }
        }
        if(mcSteps>_depth_MCstack) continue;

        // evaluating true primary vertex
        if (_do_VERTEX && _nMCPart == 0){
            PHG4VtxPoint *vtx = truthinfocontainer->GetVtx(g4particle->get_vtx_id());
            if (vtx)
            {
              _vertex_true_x = vtx->get_x();
              _vertex_true_y = vtx->get_y();
              _vertex_true_z = vtx->get_z();
            }

          
        }
        
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
            // TVector3 projvec(_mcpart_px[0],_mcpart_py[0],_mcpart_pz[0]);
            // float projeta = projvec.Eta();
        _nMCPart++;
      }
      if (Verbosity() > 0){ cout << "saved\t" << _nMCPart << "\tMC particles" << endl;}
    }
    else
    {
      if (Verbosity() > 0){ cout << PHWHERE << " PHG4TruthInfoContainer node not found on node tree" << endl;}
      return;
    }
  }
  _event_tree->Fill();
  resetBuffer();
  if (Verbosity() > 0){ cout << "EventEvaluator buffer reset" << endl;}
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
  if (_caloevalstackHCALIN) delete _caloevalstackHCALIN;
  if (_caloevalstackHCALOUT) delete _caloevalstackHCALOUT;
  if (_caloevalstackEHCAL) delete _caloevalstackEHCAL;
  if (_caloevalstackDRCALO) delete _caloevalstackDRCALO;
  if (_caloevalstackFEMC) delete _caloevalstackFEMC;
  if (_caloevalstackCEMC) delete _caloevalstackCEMC;
  if (_caloevalstackEEMC) delete _caloevalstackEEMC;
  

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

  else if (projname.find("EHCAL_0") != std::string::npos)
    return 60;
  else if (projname.find("EEMC_0") != std::string::npos)
    return 61;
  else if (projname.find("HCALIN_0") != std::string::npos)
    return 62;
  else if (projname.find("HCALOUT_0") != std::string::npos)
    return 63;
  else if (projname.find("CEMC_0") != std::string::npos)
    return 64;
  
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
    
  case 60:
    return "EHCAL";
  case 61:
    return "EEMC";
  case 62:
    return "HCALIN";
  case 63:
    return "HCALOUT";
  case 64:
    return "CEMC";
   
  default:
    return "NOTHING";
  }
}

void EventEvaluator::resetBuffer()
{
  if(_do_VERTEX){
    _vertex_x = 0;
    _vertex_y = 0;
    _vertex_z = 0;
    _vertex_NCont = 0;
    _vertex_true_x = 0;
    _vertex_true_y = 0;
    _vertex_true_z = 0;
  }
  if(_do_HITS){
    _nHitsLayers = 0;
    for (Int_t ihit = 0; ihit < _maxNHits; ihit++)
    {
      _hits_layerID[ihit] = 0;
      _hits_trueID[ihit] = 0;
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
  if(_do_DRCALO){
    _nTowers_DRCALO = 0;
    for (Int_t itow = 0; itow < _maxNTowers; itow++)
    {
      _tower_DRCALO_E[itow] = 0;
      _tower_DRCALO_NScint[itow] = 0;
      _tower_DRCALO_NCerenkov[itow] = 0;
      _tower_DRCALO_iEta[itow] = 0;
      _tower_DRCALO_iPhi[itow] = 0;
      _tower_DRCALO_trueID[itow] = 0;
    }
  }
  if(_do_TRACKS){
    _nTracks = 0;
    for (Int_t itrk = 0; itrk < _maxNTracks; itrk++)
    {
      _track_ID[itrk] = 0;
      _track_trueID[itrk] = 0;
      _track_px[itrk] = 0;
      _track_py[itrk] = 0;
      _track_pz[itrk] = 0;
    }
    if(_do_PROJECTIONS){
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
    }
  }
  if(_do_MCPARTICLES){
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
}
