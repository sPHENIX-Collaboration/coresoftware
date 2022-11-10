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
#include <calobase/RawTowerv2.h>

#include <phhepmc/PHHepMCGenEvent.h>
#include <phhepmc/PHHepMCGenEventMap.h>
#include <phhepmc/PHGenIntegral.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>

#include <phool/getClass.h>
#include <phool/phool.h>
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/PHCompositeNode.h>

#include <TFile.h>
#include <TNtuple.h>
#include <TTree.h>

#include <CLHEP/Vector/ThreeVector.h>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <HepMC/GenEvent.h>
#include <HepMC/GenVertex.h>
#pragma GCC diagnostic pop

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <set>
#include <utility>

using namespace std;

EventEvaluator::EventEvaluator(const string& name, const string& filename)
  : SubsysReco(name)
  , _do_store_event_info(false)
  , _do_HCALIN(false)
  , _do_HCALOUT(false)
  , _do_CEMC(false)
  , _do_HITS(false)
  , _do_TRACKS(false)
  , _do_CLUSTERS(false)
  , _do_VERTEX(false)
  , _do_PROJECTIONS(false)
  , _do_MCPARTICLES(false)
  , _do_HEPMC(false)
  , _do_GEOMETRY(false)
  , _ievent(0)
  , _cross_section(0)
  , _event_weight(0)
  , _n_generator_accepted(0)
  , _nHitsLayers(0)
  , _hits_layerID(0)
  , _hits_trueID(0)
  , _hits_x(0)
  , _hits_y(0)
  , _hits_z(0)
  , _hits_t(0)

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
  
  , _nTowers_CEMC(0)
  , _tower_CEMC_E(0)
  , _tower_CEMC_iEta(0)
  , _tower_CEMC_iPhi(0)
  , _tower_CEMC_trueID(0)

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

  , _nclusters_CEMC(0)
  , _cluster_CEMC_E(0)
  , _cluster_CEMC_Eta(0)
  , _cluster_CEMC_Phi(0)
  , _cluster_CEMC_NTower(0)
  , _cluster_CEMC_trueID(0)

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
  , _track_dca(0)
  , _track_dca_2d(0)
  , _track_trueID(0)
  , _track_source(0)
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
  , _mcpart_BCID(0)

  , _nHepmcp(0)
  , _hepmcp_procid(0)
  , _hepmcp_x1(NAN)
  , _hepmcp_x2(NAN)

  //  , _hepmcp_ID_parent(0)
  , _hepmcp_status(0)
  , _hepmcp_PDG(0)
  , _hepmcp_E(0)
  , _hepmcp_px(0)
  , _hepmcp_py(0)
  , _hepmcp_pz(0)
  , _hepmcp_m1(0)
  , _hepmcp_m2(0)
  , _hepmcp_BCID(0)

  , _calo_ID(0)
  , _calo_towers_N(0)
  , _calo_towers_iEta(0)
  , _calo_towers_iPhi(0)
  , _calo_towers_Eta(0)
  , _calo_towers_Phi(0)
  , _calo_towers_x(0)
  , _calo_towers_y(0)
  , _calo_towers_z(0)
  , _geometry_done(0)

  , _reco_e_threshold(0.0)
  , _reco_e_threshold_BECAL(0.0)
  , _depth_MCstack(0)
  , _caloevalstackHCALIN(nullptr)
  , _caloevalstackHCALOUT(nullptr)
  , _caloevalstackCEMC(nullptr)
  , _strict(false)
  , _event_tree(nullptr)
  , _geometry_tree(nullptr)
  , _filename(filename)
  , _tfile(nullptr)
  , _tfile_geometry(nullptr)
{
  _hits_layerID = new int[_maxNHits];
  _hits_trueID = new int[_maxNHits];
  _hits_x = new float[_maxNHits];
  _hits_y = new float[_maxNHits];
  _hits_z = new float[_maxNHits];
  _hits_t = new float[_maxNHits];

  _tower_HCALIN_E = new float[_maxNTowersCentral];
  _tower_HCALIN_iEta  = new int[_maxNTowersCentral];
  _tower_HCALIN_iPhi = new int[_maxNTowersCentral];
  _tower_HCALIN_trueID = new int[_maxNTowersCentral];
  _cluster_HCALIN_E = new float[_maxNclustersCentral];
  _cluster_HCALIN_Eta = new float[_maxNclustersCentral];
  _cluster_HCALIN_Phi = new float[_maxNclustersCentral];
  _cluster_HCALIN_NTower = new int[_maxNclustersCentral];
  _cluster_HCALIN_trueID = new int[_maxNclustersCentral];

  _tower_HCALOUT_E = new float[_maxNTowersCentral];
  _tower_HCALOUT_iEta  = new int[_maxNTowersCentral];
  _tower_HCALOUT_iPhi = new int[_maxNTowersCentral];
  _tower_HCALOUT_trueID = new int[_maxNTowersCentral];
  _cluster_HCALOUT_E = new float[_maxNclustersCentral];
  _cluster_HCALOUT_Eta = new float[_maxNclustersCentral];
  _cluster_HCALOUT_Phi = new float[_maxNclustersCentral];
  _cluster_HCALOUT_NTower = new int[_maxNclustersCentral];
  _cluster_HCALOUT_trueID = new int[_maxNclustersCentral];

  _tower_CEMC_E = new float[_maxNTowersCentral];
  _tower_CEMC_iEta = new int[_maxNTowersCentral];
  _tower_CEMC_iPhi = new int[_maxNTowersCentral];
  _tower_CEMC_trueID = new int[_maxNTowersCentral];
  _cluster_CEMC_E = new float[_maxNclustersCentral];
  _cluster_CEMC_Eta = new float[_maxNclustersCentral];
  _cluster_CEMC_Phi = new float[_maxNclustersCentral];
  _cluster_CEMC_NTower = new int[_maxNclustersCentral];
  _cluster_CEMC_trueID = new int[_maxNclustersCentral];
  
  _track_ID = new float[_maxNTracks];
  _track_trueID = new float[_maxNTracks];
  _track_px = new float[_maxNTracks];
  _track_py = new float[_maxNTracks];
  _track_pz = new float[_maxNTracks];
  _track_dca = new float[_maxNTracks];
  _track_dca_2d = new float[_maxNTracks];
  _track_source = new unsigned short[_maxNTracks];
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

  _mcpart_ID = new int[_maxNMCPart];
  _mcpart_ID_parent = new int[_maxNMCPart];
  _mcpart_PDG = new int[_maxNMCPart];
  _mcpart_E = new float[_maxNMCPart];
  _mcpart_px = new float[_maxNMCPart];
  _mcpart_py = new float[_maxNMCPart];
  _mcpart_pz = new float[_maxNMCPart];
  _mcpart_BCID = new int[_maxNMCPart];

  _hepmcp_BCID = new int[_maxNHepmcp];
  //  _hepmcp_ID_parent = new float[_maxNHepmcp];
  _hepmcp_status = new int[_maxNHepmcp];
  _hepmcp_PDG = new int[_maxNHepmcp];
  _hepmcp_E = new float[_maxNHepmcp];
  _hepmcp_px = new float[_maxNHepmcp];
  _hepmcp_py = new float[_maxNHepmcp];
  _hepmcp_pz = new float[_maxNHepmcp];
  _hepmcp_m1 = new int[_maxNHepmcp];
  _hepmcp_m2 = new int[_maxNHepmcp];

  _calo_towers_iEta = new int[_maxNTowersCalo];
  _calo_towers_iPhi = new int[_maxNTowersCalo];
  _calo_towers_Eta = new float[_maxNTowersCalo];
  _calo_towers_Phi = new float[_maxNTowersCalo];
  _calo_towers_x = new float[_maxNTowersCalo];
  _calo_towers_y = new float[_maxNTowersCalo];
  _calo_towers_z = new float[_maxNTowersCalo];
  _geometry_done = new int[20];
  for(int igem=0;igem<20;igem++) _geometry_done[igem] = 0;

}

int EventEvaluator::Init(PHCompositeNode* /*topNode*/)
{
  _ievent = 0;

  _tfile = new TFile(_filename.c_str(), "RECREATE");

  _event_tree = new TTree("event_tree", "event_tree");
  if (_do_store_event_info)
  {
    // Event level info. This isn't the most efficient way to store this info, but it's straightforward
    // within the structure of the class, so the size is small compared to the rest of the output.
    _event_tree->Branch("cross_section", &_cross_section, "cross_section/F");
    _event_tree->Branch("event_weight", &_event_weight, "event_weight/F");
    _event_tree->Branch("n_generator_accepted", &_n_generator_accepted, "n_generator_accepted/I");
  }
  // tracks and hits
  if (_do_HITS)
  {
    _event_tree->Branch("nHits", &_nHitsLayers, "nHits/I");
    _event_tree->Branch("hits_layerID", _hits_layerID, "hits_layerID[nHits]/I");
    _event_tree->Branch("hits_trueID", _hits_trueID, "hits_trueID[nHits]/I");
    _event_tree->Branch("hits_x", _hits_x, "hits_x[nHits]/F");
    _event_tree->Branch("hits_y", _hits_y, "hits_y[nHits]/F");
    _event_tree->Branch("hits_z", _hits_z, "hits_z[nHits]/F");
    _event_tree->Branch("hits_t", _hits_t, "hits_t[nHits]/F");
  }
  if (_do_TRACKS)
  {
    _event_tree->Branch("nTracks", &_nTracks, "nTracks/I");
    _event_tree->Branch("tracks_ID", _track_ID, "tracks_ID[nTracks]/F");
    _event_tree->Branch("tracks_px", _track_px, "tracks_px[nTracks]/F");
    _event_tree->Branch("tracks_py", _track_py, "tracks_py[nTracks]/F");
    _event_tree->Branch("tracks_pz", _track_pz, "tracks_pz[nTracks]/F");
    _event_tree->Branch("tracks_dca", _track_dca, "tracks_dca[nTracks]/F");
    _event_tree->Branch("tracks_dca_2d", _track_dca_2d, "tracks_dca_2d[nTracks]/F");
    _event_tree->Branch("tracks_trueID", _track_trueID, "tracks_trueID[nTracks]/F");
    _event_tree->Branch("tracks_source", _track_source, "tracks_source[nTracks]/s");
  }
  if (_do_PROJECTIONS)
  {
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
  if (_do_HCALIN)
  {
    // towers HCAL-in
    _event_tree->Branch("tower_HCALIN_N", &_nTowers_HCALIN, "tower_HCALIN_N/I");
    _event_tree->Branch("tower_HCALIN_E", _tower_HCALIN_E, "tower_HCALIN_E[tower_HCALIN_N]/F");
    _event_tree->Branch("tower_HCALIN_iEta", _tower_HCALIN_iEta, "tower_HCALIN_iEta[tower_HCALIN_N]/I");
    _event_tree->Branch("tower_HCALIN_iPhi", _tower_HCALIN_iPhi, "tower_HCALIN_iPhi[tower_HCALIN_N]/I");
    _event_tree->Branch("tower_HCALIN_trueID", _tower_HCALIN_trueID, "tower_HCALIN_trueID[tower_HCALIN_N]/I");
    if (_do_CLUSTERS)
    {
      // clusters HCAL-in
      _event_tree->Branch("cluster_HCALIN_N", &_nclusters_HCALIN, "cluster_HCALIN_N/I");
      _event_tree->Branch("cluster_HCALIN_E", _cluster_HCALIN_E, "cluster_HCALIN_E[cluster_HCALIN_N]/F");
      _event_tree->Branch("cluster_HCALIN_Eta", _cluster_HCALIN_Eta, "cluster_HCALIN_Eta[cluster_HCALIN_N]/F");
      _event_tree->Branch("cluster_HCALIN_Phi", _cluster_HCALIN_Phi, "cluster_HCALIN_Phi[cluster_HCALIN_N]/F");
      _event_tree->Branch("cluster_HCALIN_NTower", _cluster_HCALIN_NTower, "cluster_HCALIN_NTower[cluster_HCALIN_N]/I");
      _event_tree->Branch("cluster_HCALIN_trueID", _cluster_HCALIN_trueID, "cluster_HCALIN_trueID[cluster_HCALIN_N]/I");
    }
  }
  if (_do_HCALOUT)
  {
    // towers HCAL-out
    _event_tree->Branch("tower_HCALOUT_N", &_nTowers_HCALOUT, "tower_HCALOUT_N/I");
    _event_tree->Branch("tower_HCALOUT_E", _tower_HCALOUT_E, "tower_HCALOUT_E[tower_HCALOUT_N]/F");
    _event_tree->Branch("tower_HCALOUT_iEta", _tower_HCALOUT_iEta, "tower_HCALOUT_iEta[tower_HCALOUT_N]/I");
    _event_tree->Branch("tower_HCALOUT_iPhi", _tower_HCALOUT_iPhi, "tower_HCALOUT_iPhi[tower_HCALOUT_N]/I");
    _event_tree->Branch("tower_HCALOUT_trueID", _tower_HCALOUT_trueID, "tower_HCALOUT_trueID[tower_HCALOUT_N]/I");
    if (_do_CLUSTERS)
    {
      // clusters HCAL-out
      _event_tree->Branch("cluster_HCALOUT_N", &_nclusters_HCALOUT, "cluster_HCALOUT_N/I");
      _event_tree->Branch("cluster_HCALOUT_E", _cluster_HCALOUT_E, "cluster_HCALOUT_E[cluster_HCALOUT_N]/F");
      _event_tree->Branch("cluster_HCALOUT_Eta", _cluster_HCALOUT_Eta, "cluster_HCALOUT_Eta[cluster_HCALOUT_N]/F");
      _event_tree->Branch("cluster_HCALOUT_Phi", _cluster_HCALOUT_Phi, "cluster_HCALOUT_Phi[cluster_HCALOUT_N]/F");
      _event_tree->Branch("cluster_HCALOUT_NTower", _cluster_HCALOUT_NTower, "cluster_HCALOUT_NTower[cluster_HCALOUT_N]/I");
      _event_tree->Branch("cluster_HCALOUT_trueID", _cluster_HCALOUT_trueID, "cluster_HCALOUT_trueID[cluster_HCALOUT_N]/I");
    }
  }
  if (_do_CEMC)
  {
    // towers CEMC
    _event_tree->Branch("tower_CEMC_N", &_nTowers_CEMC, "tower_CEMC_N/I");
    _event_tree->Branch("tower_CEMC_E", _tower_CEMC_E, "tower_CEMC_E[tower_CEMC_N]/F");
    _event_tree->Branch("tower_CEMC_iEta", _tower_CEMC_iEta, "tower_CEMC_iEta[tower_CEMC_N]/I");
    _event_tree->Branch("tower_CEMC_iPhi", _tower_CEMC_iPhi, "tower_CEMC_iPhi[tower_CEMC_N]/I");
    _event_tree->Branch("tower_CEMC_trueID", _tower_CEMC_trueID, "tower_CEMC_trueID[tower_CEMC_N]/I");
    if (_do_CLUSTERS)
    {
      // clusters CEMC
      _event_tree->Branch("cluster_CEMC_N", &_nclusters_CEMC, "cluster_CEMC_N/I");
      _event_tree->Branch("cluster_CEMC_E", _cluster_CEMC_E, "cluster_CEMC_E[cluster_CEMC_N]/F");
      _event_tree->Branch("cluster_CEMC_Eta", _cluster_CEMC_Eta, "cluster_CEMC_Eta[cluster_CEMC_N]/F");
      _event_tree->Branch("cluster_CEMC_Phi", _cluster_CEMC_Phi, "cluster_CEMC_Phi[cluster_CEMC_N]/F");
      _event_tree->Branch("cluster_CEMC_NTower", _cluster_CEMC_NTower, "cluster_CEMC_NTower[cluster_CEMC_N]/I");
      _event_tree->Branch("cluster_CEMC_trueID", _cluster_CEMC_trueID, "cluster_CEMC_trueID[cluster_CEMC_N]/I");
    }
  }
  if (_do_VERTEX)
  {
    // vertex
    _event_tree->Branch("vertex_x", &_vertex_x, "vertex_x/F");
    _event_tree->Branch("vertex_y", &_vertex_y, "vertex_y/F");
    _event_tree->Branch("vertex_z", &_vertex_z, "vertex_z/F");
    _event_tree->Branch("vertex_NCont", &_vertex_NCont, "vertex_NCont/I");
    _event_tree->Branch("vertex_true_x", &_vertex_true_x, "vertex_true_x/F");
    _event_tree->Branch("vertex_true_y", &_vertex_true_y, "vertex_true_y/F");
    _event_tree->Branch("vertex_true_z", &_vertex_true_z, "vertex_true_z/F");
  }
  if (_do_MCPARTICLES)
  {
    // MC particles
    _event_tree->Branch("nMCPart", &_nMCPart, "nMCPart/I");
    _event_tree->Branch("mcpart_ID", _mcpart_ID, "mcpart_ID[nMCPart]/I");
    _event_tree->Branch("mcpart_ID_parent", _mcpart_ID_parent, "mcpart_ID_parent[nMCPart]/I");
    _event_tree->Branch("mcpart_PDG", _mcpart_PDG, "mcpart_PDG[nMCPart]/I");
    _event_tree->Branch("mcpart_E", _mcpart_E, "mcpart_E[nMCPart]/F");
    _event_tree->Branch("mcpart_px", _mcpart_px, "mcpart_px[nMCPart]/F");
    _event_tree->Branch("mcpart_py", _mcpart_py, "mcpart_py[nMCPart]/F");
    _event_tree->Branch("mcpart_pz", _mcpart_pz, "mcpart_pz[nMCPart]/F");
    _event_tree->Branch("mcpart_BCID", _mcpart_BCID, "mcpart_BCID[nMCPart]/I");
  }
  if (_do_HEPMC)
  {
    // MC particles
    _event_tree->Branch("nHepmcp", &_nHepmcp, "nHepmcp/I");
    _event_tree->Branch("hepmcp_procid", &_hepmcp_procid, "hepmcp_procid/I");
    _event_tree->Branch("hepmcp_x1", &_hepmcp_x1, "hepmcp_x1/F");
    _event_tree->Branch("hepmcp_x2", &_hepmcp_x2, "hepmcp_x2/F");

    //    _event_tree->Branch("hepmcp_ID_parent", _hepmcp_ID_parent, "hepmcp_ID_parent[nHepmcp]/F");
    _event_tree->Branch("hepmcp_status", _hepmcp_status, "hepmcp_status[nHepmcp]/I");
    _event_tree->Branch("hepmcp_PDG", _hepmcp_PDG, "hepmcp_PDG[nHepmcp]/I");
    _event_tree->Branch("hepmcp_E", _hepmcp_E, "hepmcp_E[nHepmcp]/F");
    _event_tree->Branch("hepmcp_px", _hepmcp_px, "hepmcp_px[nHepmcp]/F");
    _event_tree->Branch("hepmcp_py", _hepmcp_py, "hepmcp_py[nHepmcp]/F");
    _event_tree->Branch("hepmcp_pz", _hepmcp_pz, "hepmcp_pz[nHepmcp]/F");
    _event_tree->Branch("hepmcp_BCID", _hepmcp_BCID, "hepmcp_BCID[nHepmcp]/I");
    _event_tree->Branch("hepmcp_m1", _hepmcp_m1, "hepmcp_m1[nHepmcp]/I");
    _event_tree->Branch("hepmcp_m2", _hepmcp_m2, "hepmcp_m2[nHepmcp]/I");
  }


  if(_do_GEOMETRY){
    _tfile_geometry = new TFile("geometry.root", "RECREATE");

    _geometry_tree = new TTree("geometry_tree", "geometry_tree");
    // tracks and hits
    _geometry_tree->Branch("calo", &_calo_ID, "nHits/I");
    _geometry_tree->Branch("calo_towers_N",  &_calo_towers_N,"calo_towers_N/I");
    _geometry_tree->Branch("calo_towers_iEta", _calo_towers_iEta, "calo_towers_iEta[calo_towers_N]/I");
    _geometry_tree->Branch("calo_towers_iPhi", _calo_towers_iPhi, "calo_towers_iPhi[calo_towers_N]/I");
    _geometry_tree->Branch("calo_towers_Eta", _calo_towers_Eta, "calo_towers_Eta[calo_towers_N]/F");
    _geometry_tree->Branch("calo_towers_Phi", _calo_towers_Phi, "calo_towers_Phi[calo_towers_N]/F");
    _geometry_tree->Branch("calo_towers_x", _calo_towers_x, "calo_towers_x[calo_towers_N]/F");
    _geometry_tree->Branch("calo_towers_y", _calo_towers_y, "calo_towers_y[calo_towers_N]/F");
    _geometry_tree->Branch("calo_towers_z", _calo_towers_z, "calo_towers_z[calo_towers_N]/F");
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int EventEvaluator::process_event(PHCompositeNode* topNode)
{
  if (Verbosity() > 0)
  {
    cout << "entered process_event" << endl;
  }
  if (_do_HCALIN)
  {
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
  if (_do_HCALOUT)
  {
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
  if (_do_CEMC)
  {
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
  
  if (Verbosity() > 0)
  {
    cout << "loaded evalstack" << endl;
  }

  // fill the Evaluator Tree
  fillOutputNtuples(topNode);

  ++_ievent;

  return Fun4AllReturnCodes::EVENT_OK;
}

void EventEvaluator::fillOutputNtuples(PHCompositeNode* topNode)
{
  if (Verbosity() > 2)
  {
    cout << "EventEvaluator::fillOutputNtuples() entered" << endl;
  }

  //----------------------
  // fill the Event Tree
  //----------------------

  //----------------------
  // Event level info
  //---------------------
  // Extract weight info from the stored HepMC event.
  if (_do_store_event_info) {
    PHHepMCGenEventMap* hepmceventmap = findNode::getClass<PHHepMCGenEventMap>(topNode, "PHHepMCGenEventMap");
    if (hepmceventmap)
    {
      if (Verbosity() > 0)
      {
        cout << "saving event level info" << endl;
      }

      for (PHHepMCGenEventMap::ConstIter eventIter = hepmceventmap->begin();
           eventIter != hepmceventmap->end();
           ++eventIter)
      {
        PHHepMCGenEvent* hepmcevent = eventIter->second;

        if (hepmcevent)
        {
          HepMC::GenEvent* truthevent = hepmcevent->getEvent();
          if (!truthevent)
          {
            cout << PHWHERE
                 << "no evt pointer under phhepmvgeneventmap found "
                 << endl;
            return;
          }

          auto xsec = truthevent->cross_section();
          if (xsec)
          {
            _cross_section = xsec->cross_section();
          }
          // Only fill the event weight if available.
          // The overall event weight will be stored in the last entry in the vector.
          auto weights = truthevent->weights();
          if (weights.size() > 0) {
              _event_weight = weights[weights.size() - 1];
          }
        }
      }
    }
    else
    {
      if (Verbosity() > 0)
      {
        cout << PHWHERE << " PHHepMCGenEventMap node (for event level info) not found on node tree" << endl;
      }
      return;
    }

    // Retrieve the number of generator accepted events
    // Following how this was implemented in PHPythia8
    PHNodeIterator iter(topNode);
    PHCompositeNode *sumNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "RUN"));
    if (!sumNode)
    {
      cout << PHWHERE << "RUN Node missing doing nothing" << endl;
      return;
    }
    auto * integralNode = findNode::getClass<PHGenIntegral>(sumNode, "PHGenIntegral");
    if (integralNode)
    {
      _n_generator_accepted = integralNode->get_N_Generator_Accepted_Event();
    }
    else
    {
      if (Verbosity() > 0)
      {
        cout << PHWHERE << " PHGenIntegral node (for n generator accepted) not found on node tree. Continuing" << endl;
      }
    }
  }
  //----------------------
  //    VERTEX
  //----------------------
  SvtxVertexMap* vertexmap = findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMap");
  if (_do_VERTEX)
  {
    if (vertexmap)
    {
      if (!vertexmap->empty())
      {
        if (Verbosity() > 0)
        {
          cout << "saving vertex" << endl;
        }
        SvtxVertex* vertex = (vertexmap->begin())->second;

        _vertex_x = vertex->get_x();
        _vertex_y = vertex->get_y();
        _vertex_z = vertex->get_z();
        _vertex_NCont = vertex->size_tracks();
      } else {
        _vertex_x = 0.;
        _vertex_y = 0.;
        _vertex_z = 0.;
        _vertex_NCont = -1;
      }
    }
  }
  //----------------------
  //    HITS
  //----------------------
  if (_do_HITS)
  {
    if (Verbosity() > 0)
    {
      cout << "saving hits" << endl;
    }
    _nHitsLayers = 0;
    PHG4TruthInfoContainer* truthinfocontainerHits = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
    for (int iIndex = 0; iIndex < 60; ++iIndex)
    {
      // you need to add your layer name here to be saved! This has to be done
      // as we do not want to save thousands of calorimeter hits!
      if ((GetProjectionNameFromIndex(iIndex).find("MVTX") != std::string::npos) ||
          (GetProjectionNameFromIndex(iIndex).find("INTT") != std::string::npos) 
      ){
        string nodename = "G4HIT_" + GetProjectionNameFromIndex(iIndex);
        PHG4HitContainer* hits = findNode::getClass<PHG4HitContainer>(topNode, nodename);
        if (hits)
        {
          if (Verbosity() > 1)
          {
            cout << __PRETTY_FUNCTION__ << " number of hits: " << hits->size() << endl;
          }
          PHG4HitContainer::ConstRange hit_range = hits->getHits();
          for (PHG4HitContainer::ConstIterator hit_iter = hit_range.first; hit_iter != hit_range.second; hit_iter++)
          {
            if (Verbosity() > 1)
            {
              cout << __PRETTY_FUNCTION__ << " found hit with id " << hit_iter->second->get_trkid() << endl;
            }
            if(_nHitsLayers > _maxNHits){
              cout << __PRETTY_FUNCTION__ << " exceededed maximum hit array size! Please check where these hits come from!" << endl;
              break;
            }
            _hits_x[_nHitsLayers] = hit_iter->second->get_x(0);
            _hits_y[_nHitsLayers] = hit_iter->second->get_y(0);
            _hits_z[_nHitsLayers] = hit_iter->second->get_z(0);
            _hits_t[_nHitsLayers] = hit_iter->second->get_t(0);
            _hits_layerID[_nHitsLayers] = iIndex;
            if (truthinfocontainerHits)
            {
              PHG4Particle* particle = truthinfocontainerHits->GetParticle(hit_iter->second->get_trkid());

              if (particle->get_parent_id() != 0)
              {
                PHG4Particle* g4particleMother = truthinfocontainerHits->GetParticle(hit_iter->second->get_trkid());
                int mcSteps = 0;
                while (g4particleMother->get_parent_id() != 0)
                {
                  g4particleMother = truthinfocontainerHits->GetParticle(g4particleMother->get_parent_id());
                  if (g4particleMother == NULL) break;
                  mcSteps += 1;
                }
                if (mcSteps <= _depth_MCstack)
                {
                  _hits_trueID[_nHitsLayers] = hit_iter->second->get_trkid();
                }
                else
                {
                  PHG4Particle* g4particleMother2 = truthinfocontainerHits->GetParticle(hit_iter->second->get_trkid());
                  int mcSteps2 = 0;
                  while (g4particleMother2->get_parent_id() != 0 && (mcSteps2 < (mcSteps - _depth_MCstack + 1)))
                  {
                    g4particleMother2 = truthinfocontainerHits->GetParticle(g4particleMother2->get_parent_id());
                    if (g4particleMother2 == NULL){
                      break;
                    } else {
                      _hits_trueID[_nHitsLayers] = g4particleMother2->get_parent_id();
                      mcSteps2 += 1;
                    }
                  }
                }
              }
              else
              {
                _hits_trueID[_nHitsLayers] = hit_iter->second->get_trkid();
              }
            }
            _nHitsLayers++;

          }
          if (Verbosity() > 0)
          {
            cout << "saved\t" << _nHitsLayers << "\thits for " << GetProjectionNameFromIndex(iIndex) << endl;
          }
        }
        else
        {
          if (Verbosity() > 0)
          {
            cout << __PRETTY_FUNCTION__ << " could not find " << nodename << endl;
          }
          continue;
        }
      }
    }
  }
  //----------------------
  //    TOWERS HCALIN
  //----------------------
  if (_do_HCALIN)
  {
    CaloRawTowerEval* towerevalHCALIN = _caloevalstackHCALIN->get_rawtower_eval();
    _nTowers_HCALIN = 0;
    string towernodeHCALIN = "TOWER_CALIB_HCALIN";
    RawTowerContainer* towersHCALIN = findNode::getClass<RawTowerContainer>(topNode, towernodeHCALIN.c_str());
    if (towersHCALIN)
    {
      if (Verbosity() > 0)
      {
        cout << "saving HCAL towers" << endl;
      }
      string towergeomnodeHCALIN = "TOWERGEOM_HCALIN";
      RawTowerGeomContainer* towergeomHCALIN = findNode::getClass<RawTowerGeomContainer>(topNode, towergeomnodeHCALIN.c_str());
      if (towergeomHCALIN)
      {
        if(_do_GEOMETRY && !_geometry_done[kHCALIN]){
          RawTowerGeomContainer::ConstRange all_towers = towergeomHCALIN->get_tower_geometries();
          for (RawTowerGeomContainer::ConstIterator it = all_towers.first;
              it != all_towers.second; ++it)
          {
            _calo_ID = kHCALIN;
            _calo_towers_iEta[_calo_towers_N] = it->second->get_bineta();
            _calo_towers_iPhi[_calo_towers_N] = it->second->get_binphi();
            _calo_towers_Eta[_calo_towers_N] = it->second->get_eta();
            _calo_towers_Phi[_calo_towers_N] = it->second->get_phi();
            _calo_towers_x[_calo_towers_N] = it->second->get_center_x();
            _calo_towers_y[_calo_towers_N] = it->second->get_center_y();
            _calo_towers_z[_calo_towers_N] = it->second->get_center_z();
            _calo_towers_N++;
          }
          _geometry_done[kHCALIN] = 1;
          _geometry_tree->Fill();
          resetGeometryArrays();
        }
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
            }
            else
            {
              _tower_HCALIN_trueID[_nTowers_HCALIN] = -10;
            }
            _nTowers_HCALIN++;
          }
        }
      }
      else
      {
        if (Verbosity() > 0)
        {
          cout << PHWHERE << " ERROR: Can't find " << towergeomnodeHCALIN << endl;
        }
        // return;
      }
      if (Verbosity() > 0)
      {
        cout << "saved\t" << _nTowers_HCALIN << "\tHCALIN towers" << endl;
      }
    }
    else
    {
      if (Verbosity() > 0)
      {
        cout << PHWHERE << " ERROR: Can't find " << towernodeHCALIN << endl;
      }
      // return;
    }
  }
  //----------------------
  //    TOWERS HCALOUT
  //----------------------
  if (_do_HCALOUT)
  {
    CaloRawTowerEval* towerevalHCALOUT = _caloevalstackHCALOUT->get_rawtower_eval();
    _nTowers_HCALOUT = 0;
    string towernodeHCALOUT = "TOWER_CALIB_HCALOUT";
    RawTowerContainer* towersHCALOUT = findNode::getClass<RawTowerContainer>(topNode, towernodeHCALOUT.c_str());
    if (towersHCALOUT)
    {
      if (Verbosity() > 0)
      {
        cout << "saving HCAL towers" << endl;
      }
      string towergeomnodeHCALOUT = "TOWERGEOM_HCALOUT";
      RawTowerGeomContainer* towergeomHCALOUT = findNode::getClass<RawTowerGeomContainer>(topNode, towergeomnodeHCALOUT.c_str());
      if (towergeomHCALOUT)
      {
        if(_do_GEOMETRY && !_geometry_done[kHCALOUT]){
          RawTowerGeomContainer::ConstRange all_towers = towergeomHCALOUT->get_tower_geometries();
          for (RawTowerGeomContainer::ConstIterator it = all_towers.first;
              it != all_towers.second; ++it)
          {
            _calo_ID = kHCALOUT;
            _calo_towers_iEta[_calo_towers_N] = it->second->get_bineta();
            _calo_towers_iPhi[_calo_towers_N] = it->second->get_binphi();
            _calo_towers_Eta[_calo_towers_N] = it->second->get_eta();
            _calo_towers_Phi[_calo_towers_N] = it->second->get_phi();
            _calo_towers_x[_calo_towers_N] = it->second->get_center_x();
            _calo_towers_y[_calo_towers_N] = it->second->get_center_y();
            _calo_towers_z[_calo_towers_N] = it->second->get_center_z();
            _calo_towers_N++;
          }
          _geometry_done[kHCALOUT] = 1;
          _geometry_tree->Fill();
          resetGeometryArrays();
        }
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
            }
            else
            {
              _tower_HCALOUT_trueID[_nTowers_HCALOUT] = -10;
            }
            _nTowers_HCALOUT++;
          }
        }
      }
      else
      {
        if (Verbosity() > 0)
        {
          cout << PHWHERE << " ERROR: Can't find " << towergeomnodeHCALOUT << endl;
        }
        // return;
      }
      if (Verbosity() > 0)
      {
        cout << "saved\t" << _nTowers_HCALOUT << "\tHCALOUT towers" << endl;
      }
    }
    else
    {
      if (Verbosity() > 0)
      {
        cout << PHWHERE << " ERROR: Can't find " << towernodeHCALOUT << endl;
      }
      // return;
    }
  }
  //----------------------
  //    TOWERS CEMC
  //----------------------
  if (_do_CEMC)
  {
    CaloRawTowerEval* towerevalCEMC = _caloevalstackCEMC->get_rawtower_eval();
    _nTowers_CEMC = 0;
    string towernodeCEMC = "TOWER_CALIB_CEMC";
    RawTowerContainer* towersCEMC = findNode::getClass<RawTowerContainer>(topNode, towernodeCEMC.c_str());
    if (towersCEMC)
    {
      if (Verbosity() > 0)
      {
        cout << "saving EMC towers" << endl;
      }
      string towergeomnodeCEMC = "TOWERGEOM_CEMC";
      RawTowerGeomContainer* towergeom = findNode::getClass<RawTowerGeomContainer>(topNode, towergeomnodeCEMC.c_str());
      if (towergeom)
      {
        if(_do_GEOMETRY && !_geometry_done[kCEMC]){
          RawTowerGeomContainer::ConstRange all_towers = towergeom->get_tower_geometries();
          for (RawTowerGeomContainer::ConstIterator it = all_towers.first;
              it != all_towers.second; ++it)
          {
            _calo_ID = kCEMC;
            _calo_towers_iEta[_calo_towers_N] = it->second->get_bineta();
            _calo_towers_iPhi[_calo_towers_N] = it->second->get_binphi();
            _calo_towers_Eta[_calo_towers_N] = it->second->get_eta();
            _calo_towers_Phi[_calo_towers_N] = it->second->get_phi();
            _calo_towers_x[_calo_towers_N] = it->second->get_center_x();
            _calo_towers_y[_calo_towers_N] = it->second->get_center_y();
            _calo_towers_z[_calo_towers_N] = it->second->get_center_z();
            _calo_towers_N++;
          }
          _geometry_done[kCEMC] = 1;
          _geometry_tree->Fill();
          resetGeometryArrays();
        }
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
            }
            else
            {
              _tower_CEMC_trueID[_nTowers_CEMC] = -10;
            }
            _nTowers_CEMC++;
          }
        }
      }
      else
      {
        if (Verbosity() > 0)
        {
          cout << PHWHERE << " ERROR: Can't find " << towergeomnodeCEMC << endl;
        }
        // return;
      }
      if (Verbosity() > 0)
      {
        cout << "saved\t" << _nTowers_CEMC << "\tCEMC towers" << endl;
      }
    }
    else
    {
      if (Verbosity() > 0)
      {
        cout << PHWHERE << " ERROR: Can't find " << towernodeCEMC << endl;
      }
      // return;
    }
  }
  
  //------------------------
  // CLUSTERS HCALIN
  //------------------------
  if (_do_HCALIN && _do_CLUSTERS)
  {
    CaloRawClusterEval* clusterevalHCALIN = _caloevalstackHCALIN->get_rawcluster_eval();
    _nclusters_HCALIN = 0;
    if (Verbosity() > 1)
    {
      cout << "CaloEvaluator::filling gcluster ntuple..." << endl;
    }

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
            _cluster_HCALIN_Eta[_nclusters_HCALIN] = RawClusterUtility::GetPseudorapidity(*cluster, CLHEP::Hep3Vector(0, 0, 0));;
        }
        else
          _cluster_HCALIN_Eta[_nclusters_HCALIN] = RawClusterUtility::GetPseudorapidity(*cluster, CLHEP::Hep3Vector(0, 0, 0));;

        PHG4Particle* primary = clusterevalHCALIN->max_truth_primary_particle_by_energy(cluster);

        if (primary)
        {
          _cluster_HCALIN_trueID[_nclusters_HCALIN] = primary->get_track_id();
        }
        else
        {
          _cluster_HCALIN_trueID[_nclusters_HCALIN] = -10;
        }

        _nclusters_HCALIN++;
      }
    }
    else
    {
      cerr << PHWHERE << " ERROR: Can't find " << clusternodeHCALIN << endl;
      // return;
    }
    if (Verbosity() > 0){ cout << "saved\t" << _nclusters_HCALIN << "\tHCALIN clusters" << endl;}
  }
  //------------------------
  // CLUSTERS HCALOUT
  //------------------------
  if (_do_HCALOUT && _do_CLUSTERS)
  {
    CaloRawClusterEval* clusterevalHCALOUT = _caloevalstackHCALOUT->get_rawcluster_eval();
    _nclusters_HCALOUT = 0;
    if (Verbosity() > 1)
    {
      cout << "CaloEvaluator::filling gcluster ntuple..." << endl;
    }

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
            _cluster_HCALOUT_Eta[_nclusters_HCALOUT] = RawClusterUtility::GetPseudorapidity(*cluster, CLHEP::Hep3Vector(0, 0, 0));;
        }
        else
          _cluster_HCALOUT_Eta[_nclusters_HCALOUT] = RawClusterUtility::GetPseudorapidity(*cluster, CLHEP::Hep3Vector(0, 0, 0));;

        PHG4Particle* primary = clusterevalHCALOUT->max_truth_primary_particle_by_energy(cluster);

        if (primary)
        {
          _cluster_HCALOUT_trueID[_nclusters_HCALOUT] = primary->get_track_id();
        }
        else
        {
          _cluster_HCALOUT_trueID[_nclusters_HCALOUT] = -10;
        }

        _nclusters_HCALOUT++;
      }
    }
    else
    {
      cerr << PHWHERE << " ERROR: Can't find " << clusternodeHCALOUT << endl;
      // return;
    }
    if (Verbosity() > 0){ cout << "saved\t" << _nclusters_HCALOUT << "\tHCALOUT clusters" << endl;}
  }
  //------------------------
  // CLUSTERS CEMC
  //------------------------
  if (_do_CEMC && _do_CLUSTERS)
  {
    CaloRawClusterEval* clusterevalCEMC = _caloevalstackCEMC->get_rawcluster_eval();
    _nclusters_CEMC = 0;
    if (Verbosity() > 1)
    {
      cout << "CaloEvaluator::filling gcluster ntuple..." << endl;
    }

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
            _cluster_CEMC_Eta[_nclusters_CEMC] = RawClusterUtility::GetPseudorapidity(*cluster, CLHEP::Hep3Vector(0, 0, 0));;
        }
        else
          _cluster_CEMC_Eta[_nclusters_CEMC] = RawClusterUtility::GetPseudorapidity(*cluster, CLHEP::Hep3Vector(0, 0, 0));;

        PHG4Particle* primary = clusterevalCEMC->max_truth_primary_particle_by_energy(cluster);

        if (primary)
        {
          _cluster_CEMC_trueID[_nclusters_CEMC] = primary->get_track_id();
        }
        else
        {
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
    if (Verbosity() > 0){ cout << "saved\t" << _nclusters_CEMC << "\tCEMC clusters" << endl;}
  }

  //------------------------
  // TRACKS
  //------------------------
  if (_do_TRACKS)
  {
    _nTracks = 0;
    _nProjections = 0;
    // Loop over track maps, identifiy each source.
    // Although this configuration is fixed here, it doesn't require multiple sources.
    // It will only store them if they're available.
    std::vector<std::pair<std::string, TrackSource_t>> trackMapInfo = {
        {"TrackMap", TrackSource_t::all},
        {"TrackMapInner", TrackSource_t::inner}};
    bool foundAtLeastOneTrackSource = false;
    for (const auto& trackMapInfo : trackMapInfo)
    {
      SvtxTrackMap* trackmap = findNode::getClass<SvtxTrackMap>(topNode, trackMapInfo.first);
      if (trackmap)
      {
        foundAtLeastOneTrackSource = true;
        int nTracksInASource = 0;
        if (Verbosity() > 0)
        {
          cout << "saving tracks for track map: " << trackMapInfo.first << endl;
        }
        for (SvtxTrackMap::ConstIter track_itr = trackmap->begin(); track_itr != trackmap->end(); track_itr++)
        {
          SvtxTrack_FastSim* track = dynamic_cast<SvtxTrack_FastSim*>(track_itr->second);
          if (track)
          {
            _track_ID[_nTracks] = track->get_id();
            _track_px[_nTracks] = track->get_px();
            _track_py[_nTracks] = track->get_py();
            _track_pz[_nTracks] = track->get_pz();
            // Ideally, would be dca3d_xy and dca3d_z, but these don't seem to be calculated properly in the
            // current (June 2021) simulations (they return NaN). So we take dca (seems to be ~ the 3d distance)
            // and dca_2d (seems to be ~ the distance in the transverse plane).
            // The names of the branches are based on the method names.
            _track_dca[_nTracks] = static_cast<float>(track->get_dca());
            _track_dca_2d[_nTracks] = static_cast<float>(track->get_dca2d());
            _track_trueID[_nTracks] = track->get_truth_track_id();
            _track_source[_nTracks] = static_cast<unsigned short>(trackMapInfo.second);
            if (_do_PROJECTIONS)
            {
              // find projections
              for (SvtxTrack::ConstStateIter trkstates = track->begin_states(); trkstates != track->end_states(); ++trkstates)
              {
                if (Verbosity() > 1)
                {
                  cout << __PRETTY_FUNCTION__ << " processing " << trkstates->second->get_name() << endl;
                }
                string trackStateName = trkstates->second->get_name();
                if (Verbosity() > 1)
                {
                  cout << __PRETTY_FUNCTION__ << " found " << trkstates->second->get_name() << endl;
                }
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
                    if (Verbosity() > 1)
                    {
                      cout << __PRETTY_FUNCTION__ << " number of hits: " << hits->size() << endl;
                    }
                    PHG4HitContainer::ConstRange hit_range = hits->getHits();
                    for (PHG4HitContainer::ConstIterator hit_iter = hit_range.first; hit_iter != hit_range.second; hit_iter++)
                    {
                      if (Verbosity() > 1)
                      {
                        cout << __PRETTY_FUNCTION__ << " checking hit id " << hit_iter->second->get_trkid() << " against " << track->get_truth_track_id() << endl;
                      }
                      if (hit_iter->second->get_trkid() - track->get_truth_track_id() == 0)
                      {
                        if (Verbosity() > 1)
                        {
                          cout << __PRETTY_FUNCTION__ << " found hit with id " << hit_iter->second->get_trkid() << endl;
                        }
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
                    if (Verbosity() > 1)
                    {
                      cout << __PRETTY_FUNCTION__ << " could not find " << nodename << endl;
                    }
                    continue;
                  }
                  _nProjections++;
                }
              }
            }
            _nTracks++;
            nTracksInASource++;
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
        if (Verbosity() > 0)
        {
          cout << "saved\t" << nTracksInASource << "\ttracks from track map " << trackMapInfo.first << ". Total saved tracks: " << _nTracks << endl;
        }
      }
      else
      {
        if (Verbosity() > 0)
        {
          cout << PHWHERE << "SvtxTrackMap node with name '" << trackMapInfo.first << "' not found on node tree" << endl;
        }
      }
    }
    if (foundAtLeastOneTrackSource == false) {
      cout << PHWHERE << "Requested tracks, but found no sources on node tree. Returning" << endl;
      return;
    }
  }
  //------------------------
  // MC PARTICLES
  //------------------------
  _nMCPart = 0;
  if (_do_MCPARTICLES)
  {
    PHG4TruthInfoContainer* truthinfocontainer = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
    if (truthinfocontainer)
    {
      if (Verbosity() > 0)
      {
        cout << "saving MC particles" << endl;
      }
      //GetParticleRange for all particles
      //GetPrimaryParticleRange for primary particles
      PHG4TruthInfoContainer::ConstRange range = truthinfocontainer->GetParticleRange();
      for (PHG4TruthInfoContainer::ConstIterator truth_itr = range.first; truth_itr != range.second; ++truth_itr)
      {
        PHG4Particle* g4particle = truth_itr->second;
        if (!g4particle) continue;

        int mcSteps = 0;
        PHG4Particle* g4particleMother = truth_itr->second;
        if (g4particle->get_parent_id() != 0)
        {
          while (g4particleMother->get_parent_id() != 0)
          {
            g4particleMother = truthinfocontainer->GetParticle(g4particleMother->get_parent_id());
            if (g4particleMother == NULL) break;
            mcSteps += 1;
          }
        }
        if (mcSteps > _depth_MCstack) continue;

        // evaluating true primary vertex
        if (_do_VERTEX && _nMCPart == 0)
        {
          PHG4VtxPoint* vtx = truthinfocontainer->GetVtx(g4particle->get_vtx_id());
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

        //using the e threshold also for the truth particles gets rid of all the low energy secondary particles
        if (g4particle->get_e() < _reco_e_threshold) continue;

        _mcpart_ID[_nMCPart] = g4particle->get_track_id();
        _mcpart_ID_parent[_nMCPart] = g4particle->get_parent_id();
        _mcpart_PDG[_nMCPart] = g4particle->get_pid();
        _mcpart_E[_nMCPart] = g4particle->get_e();
        _mcpart_px[_nMCPart] = g4particle->get_px();
        _mcpart_py[_nMCPart] = g4particle->get_py();
        _mcpart_pz[_nMCPart] = g4particle->get_pz();
        //BCID added for G4Particle --  HEPMC particle matching
        _mcpart_BCID[_nMCPart] = g4particle->get_barcode();
        // TVector3 projvec(_mcpart_px[0],_mcpart_py[0],_mcpart_pz[0]);
        // float projeta = projvec.Eta();
        _nMCPart++;
      }
      if (Verbosity() > 0)
      {
        cout << "saved\t" << _nMCPart << "\tMC particles" << endl;
      }
    }
    else
    {
      if (Verbosity() > 0)
      {
        cout << PHWHERE << " PHG4TruthInfoContainer node not found on node tree" << endl;
      }
      return;
    }
  }

  //------------------------
  // HEPMC
  //------------------------
  _nHepmcp = 0;
  if (_do_HEPMC)
  {
    PHHepMCGenEventMap* hepmceventmap = findNode::getClass<PHHepMCGenEventMap>(topNode, "PHHepMCGenEventMap");
    if (hepmceventmap)
    {
      if (Verbosity() > 0)
      {
        cout << "saving HepMC output" << endl;
      }
      if (Verbosity() > 0)
      {
        hepmceventmap->Print();
      }

      for (PHHepMCGenEventMap::ConstIter eventIter = hepmceventmap->begin();
           eventIter != hepmceventmap->end();
           ++eventIter)
      {
        PHHepMCGenEvent* hepmcevent = eventIter->second;

        if (hepmcevent)
        {
          HepMC::GenEvent* truthevent = hepmcevent->getEvent();
          if (!truthevent)
          {
            cout << PHWHERE
                 << "no evt pointer under phhepmvgeneventmap found "
                 << endl;
            return;
          }

          HepMC::PdfInfo* pdfinfo = truthevent->pdf_info();

          //     m_partid1 = pdfinfo->id1();
          // m_partid2 = pdfinfo->id2();
          _hepmcp_x1 = pdfinfo->x1();
          _hepmcp_x2 = pdfinfo->x2();

          // m_mpi = truthevent->mpi();

          _hepmcp_procid = truthevent->signal_process_id();

          if (Verbosity() > 2)
          {
            cout << " Iterating over an event" << endl;
          }
          for (HepMC::GenEvent::particle_const_iterator iter = truthevent->particles_begin();
               iter != truthevent->particles_end();
               ++iter)
          {
            _hepmcp_E[_nHepmcp] = (*iter)->momentum().e();
            _hepmcp_PDG[_nHepmcp] = (*iter)->pdg_id();
            _hepmcp_px[_nHepmcp] = (*iter)->momentum().px();
            _hepmcp_py[_nHepmcp] = (*iter)->momentum().py();
            _hepmcp_pz[_nHepmcp] = (*iter)->momentum().pz();
            _hepmcp_status[_nHepmcp] = (*iter)->status();
            _hepmcp_BCID[_nHepmcp] = (*iter)->barcode();
            _hepmcp_m2[_nHepmcp] = 0;
            _hepmcp_m1[_nHepmcp] = 0;
            if ((*iter)->production_vertex())
            {
              for (HepMC::GenVertex::particle_iterator mother = (*iter)->production_vertex()->particles_begin(HepMC::parents);
                   mother != (*iter)->production_vertex()->particles_end(HepMC::parents);
                   ++mother)
              {
                _hepmcp_m2[_nHepmcp] = (*mother)->barcode();
                if (_hepmcp_m1[_nHepmcp] == 0)
                  _hepmcp_m1[_nHepmcp] = (*mother)->barcode();
              }
            }
            if (Verbosity() > 2) cout << "nHepmcp " << _nHepmcp << "\tPDG " << _hepmcp_PDG[_nHepmcp] << "\tEnergy " << _hepmcp_E[_nHepmcp] << "\tbarcode " << _hepmcp_BCID[_nHepmcp] << "\tMother1 " << _hepmcp_m1[_nHepmcp]<< "\tMother2 " << _hepmcp_m2[_nHepmcp] << endl;
            _nHepmcp++;
          }
        }
      }
    }
    else
    {
      if (Verbosity() > 0)
      {
        cout << PHWHERE << " PHHepMCGenEventMap node not found on node tree" << endl;
      }
      return;
    }
  }  //hepmc

  _event_tree->Fill();

  if (Verbosity() > 0){ cout << "Resetting buffer ..." << endl;}
  resetBuffer();
  if (Verbosity() > 0)
  {
    cout << "EventEvaluator buffer reset" << endl;
  }
  return;
}

int EventEvaluator::End(PHCompositeNode* /*topNode*/)
{
  _tfile->cd();

  _event_tree->Write();

  _tfile->Close();

  delete _tfile;

  if(_do_GEOMETRY){
    _tfile_geometry->cd();

    _geometry_tree->Write();

    _tfile_geometry->Close();

    delete _tfile_geometry;
  }
  if (Verbosity() > 0)
  {
    cout << "========================= " << Name() << "::End() ============================" << endl;
    cout << " " << _ievent << " events of output written to: " << _filename << endl;
    cout << "===========================================================================" << endl;
  }

  if (_caloevalstackHCALIN) delete _caloevalstackHCALIN;
  if (_caloevalstackHCALOUT) delete _caloevalstackHCALOUT;
  if (_caloevalstackCEMC) delete _caloevalstackCEMC;

  return Fun4AllReturnCodes::EVENT_OK;
}

int EventEvaluator::GetProjectionIndex(std::string projname)
{
  if (projname.find("HCALIN") != std::string::npos)
    return 1;
  else if (projname.find("HCALOUT") != std::string::npos)
    return 2;
  else if (projname.find("CEMC") != std::string::npos)
    return 3;
  else
    return -1;
  return -1;
}

std::string EventEvaluator::GetProjectionNameFromIndex(int projindex)
{
  switch (projindex)
  {
    case 1:
      return "HCALIN";
    case 2:
      return "HCALOUT";
    case 3:
      return "CEMC";
    default:
      return "NOTHING";
  }
}

void EventEvaluator::resetGeometryArrays()
{
  for (Int_t igeo = 0; igeo < _calo_towers_N; igeo++)
    {
      _calo_towers_iEta[_calo_towers_N] = -10000;
      _calo_towers_iPhi[_calo_towers_N] = -10000;
      _calo_towers_Eta[_calo_towers_N] = -10000;
      _calo_towers_Phi[_calo_towers_N] = -10000;
      _calo_towers_x[_calo_towers_N] = -10000;
      _calo_towers_y[_calo_towers_N] = -10000;
      _calo_towers_z[_calo_towers_N] = -10000;
    }
    _calo_ID = -1;
    _calo_towers_N = 0;
}
void EventEvaluator::resetBuffer()
{
  if (_do_store_event_info)
  {
    _cross_section = 0;
    _event_weight = 0;
    _n_generator_accepted = 0;
    if (Verbosity() > 0){ cout << "\t... event info variables reset" << endl;}
  }
  if (_do_VERTEX)
  {
    _vertex_x = -1000;
    _vertex_y = -1000;
    _vertex_z = -1000;
    _vertex_NCont = 0;
    _vertex_true_x = -1000;
    _vertex_true_y = -1000;
    _vertex_true_z = -1000;
    if (Verbosity() > 0){ cout << "\t... vertex variables reset" << endl;}
  }
  if (_do_HITS)
  {
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
    if (Verbosity() > 0){ cout << "\t... hit variables reset" << endl;}
  }
  if(_do_CEMC){
    _nTowers_CEMC = 0;
    for (Int_t itow = 0; itow < _maxNTowersCentral; itow++)
    {
      _tower_CEMC_E[itow] = 0;
      _tower_CEMC_iEta[itow] = 0;
      _tower_CEMC_iPhi[itow] = 0;
      _tower_CEMC_trueID[itow] = 0;
    }
    if(_do_CLUSTERS){
      _nclusters_CEMC = 0;
      for (Int_t itow = 0; itow < _maxNclustersCentral; itow++)
      {
        _cluster_CEMC_E[itow] = 0;
        _cluster_CEMC_Eta[itow] = 0;
        _cluster_CEMC_Phi[itow] = 0;
        _cluster_CEMC_NTower[itow] = 0;
        _cluster_CEMC_trueID[itow] = 0;
      }
    }
    if (Verbosity() > 0){ cout << "\t... CEMC variables reset" << endl;}
  }
  if(_do_HCALIN){
    _nTowers_HCALIN = 0;
    for (Int_t itow = 0; itow < _maxNTowersCentral; itow++)
    {
      _tower_HCALIN_E[itow] = 0;
      _tower_HCALIN_iEta[itow] = 0;
      _tower_HCALIN_iPhi[itow] = 0;
      _tower_HCALIN_trueID[itow] = 0;
    }
    if(_do_CLUSTERS){
      _nclusters_HCALIN = 0;
      for (Int_t itow = 0; itow < _maxNclustersCentral; itow++)
      {
        _cluster_HCALIN_E[itow] = 0;
        _cluster_HCALIN_Eta[itow] = 0;
        _cluster_HCALIN_Phi[itow] = 0;
        _cluster_HCALIN_NTower[itow] = 0;
        _cluster_HCALIN_trueID[itow] = 0;
      }
    }
    if (Verbosity() > 0){ cout << "\t... HCALIN variables reset" << endl;}
  }
  if(_do_HCALOUT){
    if (Verbosity() > 0){ cout << "\t... resetting HCALOUT variables" << endl;}
    _nTowers_HCALOUT = 0;
    for (Int_t itow = 0; itow < _maxNTowersCentral; itow++)
    {
      _tower_HCALOUT_E[itow] = 0;
      _tower_HCALOUT_iEta[itow] = 0;
      _tower_HCALOUT_iPhi[itow] = 0;
      _tower_HCALOUT_trueID[itow] = 0;
    }
    if(_do_CLUSTERS){
      _nclusters_HCALOUT = 0;
      for (Int_t itow = 0; itow < _maxNclustersCentral; itow++)
      {
        _cluster_HCALOUT_E[itow] = 0;
        _cluster_HCALOUT_Eta[itow] = 0;
        _cluster_HCALOUT_Phi[itow] = 0;
        _cluster_HCALOUT_NTower[itow] = 0;
        _cluster_HCALOUT_trueID[itow] = 0;
      }
    }
    if (Verbosity() > 0){ cout << "\t... HCALOUT variables reset" << endl;}
  }
  if (_do_TRACKS)
  {
    if (Verbosity() > 0){ cout << "\t... resetting Track variables" << endl;}
    _nTracks = 0;
    for (Int_t itrk = 0; itrk < _maxNTracks; itrk++)
    {
      _track_ID[itrk] = 0;
      _track_trueID[itrk] = 0;
      _track_px[itrk] = 0;
      _track_py[itrk] = 0;
      _track_pz[itrk] = 0;
      _track_dca[itrk] = 0;
      _track_dca_2d[itrk] = 0;
      _track_source[itrk] = 0;
    }
    if (_do_PROJECTIONS)
    {
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
    if (Verbosity() > 0){ cout << "\t... track variables reset" << endl;}
  }
  if (_do_MCPARTICLES)
  {
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
      _mcpart_BCID[imcpart] = -10;
    }
  }

  if (_do_HEPMC)
  {
    _nHepmcp = 0;
    for (Int_t iHepmcp = 0; iHepmcp < _maxNHepmcp; iHepmcp++)
    {
      _hepmcp_E[iHepmcp] = 0;
      _hepmcp_PDG[iHepmcp] = 0;
      _hepmcp_px[iHepmcp] = 0;
      _hepmcp_py[iHepmcp] = 0;
      _hepmcp_pz[iHepmcp] = 0;
      _hepmcp_status[iHepmcp] = -10;
      _hepmcp_BCID[iHepmcp] = 0;
      _hepmcp_m2[iHepmcp] = 0;
      _hepmcp_m1[iHepmcp] = 0;
    }
    if (Verbosity() > 0){ cout << "\t... MC variables reset" << endl;}
  }
}
