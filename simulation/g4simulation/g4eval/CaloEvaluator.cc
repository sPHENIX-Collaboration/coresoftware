#include "CaloEvaluator.h"

#include "CaloEvalStack.h"
#include "CaloRawClusterEval.h"
#include "CaloRawTowerEval.h"
#include "CaloTruthEval.h"

#include <g4main/PHG4Particle.h>
#include <g4main/PHG4Shower.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4VtxPoint.h>

#include <g4vertex/GlobalVertex.h>
#include <g4vertex/GlobalVertexMap.h>

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

CaloEvaluator::CaloEvaluator(const string& name, const string& caloname, const string& filename)
  : SubsysReco(name)
  , _caloname(caloname)
  , _ievent(0)
  , _towerID_debug(0)
  , _ieta_debug(0)
  , _iphi_debug(0)
  , _eta_debug(0)
  , _phi_debug(0)
  , _e_debug(0)
  , _x_debug(0)
  , _y_debug(0)
  , _z_debug(0)
  , _truth_trace_embed_flags()
  , _truth_e_threshold(0.0)
  ,  // 0 GeV before reco is traced
  _reco_e_threshold(0.0)
  ,  // 0 GeV before reco is traced
  _caloevalstack(nullptr)
  , _strict(false)
  , _do_gpoint_eval(true)
  , _do_gshower_eval(true)
  , _do_tower_eval(true)
  , _do_cluster_eval(true)
  , _ntp_gpoint(nullptr)
  , _ntp_gshower(nullptr)
  , _ntp_tower(nullptr)
  , _tower_debug(nullptr)
  , _ntp_cluster(nullptr)
  , _filename(filename)
  , _tfile(nullptr)
{
}

int CaloEvaluator::Init(PHCompositeNode* /*topNode*/)
{
  _ievent = 0;

  _tfile = new TFile(_filename.c_str(), "RECREATE");

  if (_do_gpoint_eval) _ntp_gpoint = new TNtuple("ntp_gpoint", "primary vertex => best (first) vertex",
                                                 "event:gvx:gvy:gvz:"
                                                 "vx:vy:vz");

  if (_do_gshower_eval) _ntp_gshower = new TNtuple("ntp_gshower", "truth shower => best cluster",
                                                   "event:gparticleID:gflavor:gnhits:"
                                                   "geta:gphi:ge:gpt:gvx:gvy:gvz:gembed:gedep:"
                                                   "clusterID:ntowers:eta:x:y:z:phi:e:efromtruth");

  //Barak: Added TTree to will allow the TowerID to be set correctly as integer
  if (_do_tower_eval)
  {
    _ntp_tower = new TNtuple("ntp_tower", "tower => max truth primary",
                             "event:towerID:ieta:iphi:eta:phi:e:x:y:z:"
                             "gparticleID:gflavor:gnhits:"
                             "geta:gphi:ge:gpt:gvx:gvy:gvz:"
                             "gembed:gedep:"
                             "efromtruth");

    //Make Tree
    _tower_debug = new TTree("tower_debug", "tower => max truth primary");

    _tower_debug->Branch("event", &_ievent, "event/I");
    _tower_debug->Branch("towerID", &_towerID_debug, "towerID/I");
    _tower_debug->Branch("ieta", &_ieta_debug, "ieta/I");
    _tower_debug->Branch("iphi", &_iphi_debug, "iphi/I");
    _tower_debug->Branch("eta", &_eta_debug, "eta/F");
    _tower_debug->Branch("phi", &_phi_debug, "phi/F");
    _tower_debug->Branch("e", &_e_debug, "e/F");
    _tower_debug->Branch("x", &_x_debug, "x/F");
    _tower_debug->Branch("y", &_y_debug, "y/F");
    _tower_debug->Branch("z", &_z_debug, "z/F");
  }

  if (_do_cluster_eval) _ntp_cluster = new TNtuple("ntp_cluster", "cluster => max truth primary",
                                                   "event:clusterID:ntowers:eta:x:y:z:phi:e:"
                                                   "gparticleID:gflavor:gnhits:"
                                                   "geta:gphi:ge:gpt:gvx:gvy:gvz:"
                                                   "gembed:gedep:"
                                                   "efromtruth");

  return Fun4AllReturnCodes::EVENT_OK;
}

int CaloEvaluator::process_event(PHCompositeNode* topNode)
{
  if (!_caloevalstack)
  {
    _caloevalstack = new CaloEvalStack(topNode, _caloname);
    _caloevalstack->set_strict(_strict);
    _caloevalstack->set_verbosity(Verbosity() + 1);
  }
  else
  {
    _caloevalstack->next_event(topNode);
  }

  //-----------------------------------
  // print what is coming into the code
  //-----------------------------------

  printInputInfo(topNode);

  //---------------------------
  // fill the Evaluator NTuples
  //---------------------------

  fillOutputNtuples(topNode);

  //--------------------------------------------------
  // Print out the ancestry information for this event
  //--------------------------------------------------

  printOutputInfo(topNode);

  ++_ievent;

  return Fun4AllReturnCodes::EVENT_OK;
}

int CaloEvaluator::End(PHCompositeNode* /*topNode*/)
{
  _tfile->cd();

  if (_do_gpoint_eval) _ntp_gpoint->Write();
  if (_do_gshower_eval) _ntp_gshower->Write();
  if (_do_tower_eval)
  {
    _ntp_tower->Write();
    _tower_debug->Write();
  }
  if (_do_cluster_eval) _ntp_cluster->Write();

  _tfile->Close();

  delete _tfile;

  if (Verbosity() > 0)
  {
    cout << "========================= " << Name() << "::End() ============================" << endl;
    cout << " " << _ievent << " events of output written to: " << _filename << endl;
    cout << "===========================================================================" << endl;
  }

  if (_caloevalstack) delete _caloevalstack;

  return Fun4AllReturnCodes::EVENT_OK;
}

void CaloEvaluator::printInputInfo(PHCompositeNode* topNode)
{
  if (Verbosity() > 2) cout << "CaloEvaluator::printInputInfo() entered" << endl;

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
      cerr << PHWHERE << " ERROR: Can't find G4TruthInfo" << endl;
      exit(-1);
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

void CaloEvaluator::printOutputInfo(PHCompositeNode* topNode)
{
  if (Verbosity() > 2) cout << "CaloEvaluator::printOutputInfo() entered" << endl;

  CaloRawClusterEval* clustereval = _caloevalstack->get_rawcluster_eval();
  CaloTruthEval* trutheval = _caloevalstack->get_truth_eval();

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
      cerr << PHWHERE << " ERROR: Can't find G4TruthInfo" << endl;
      exit(-1);
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
      float gvx = vtx->get_x();
      float gvy = vtx->get_y();
      float gvz = vtx->get_z();

      cout << " vtrue = (";
      cout.width(5);
      cout << gvx;
      cout << ",";
      cout.width(5);
      cout << gvy;
      cout << ",";
      cout.width(5);
      cout << gvz;
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

void CaloEvaluator::fillOutputNtuples(PHCompositeNode* topNode)
{
  if (Verbosity() > 2) cout << "CaloEvaluator::fillOutputNtuples() entered" << endl;

  CaloRawClusterEval* clustereval = _caloevalstack->get_rawcluster_eval();
  CaloRawTowerEval* towereval = _caloevalstack->get_rawtower_eval();
  CaloTruthEval* trutheval = _caloevalstack->get_truth_eval();

  //----------------------
  // fill the Event NTuple
  //----------------------

  if (_do_gpoint_eval)
  {
    // need things off of the DST...
    PHG4TruthInfoContainer* truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
    if (!truthinfo)
    {
      cerr << PHWHERE << " ERROR: Can't find G4TruthInfo" << endl;
      exit(-1);
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

    float gpoint_data[7] = {(float) _ievent,
                            gvx,
                            gvy,
                            gvz,
                            vx,
                            vy,
                            vz};

    _ntp_gpoint->Fill(gpoint_data);
  }

  //------------------------
  // fill the Gshower NTuple
  //------------------------

  if (_ntp_gshower)
  {
    if (Verbosity() > 1) cout << Name() << " CaloEvaluator::filling gshower ntuple..." << endl;

    GlobalVertexMap* vertexmap = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");

    PHG4TruthInfoContainer* truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
    if (!truthinfo)
    {
      cerr << PHWHERE << " ERROR: Can't find G4TruthInfo" << endl;
      exit(-1);
    }

    PHG4TruthInfoContainer::ConstRange range = truthinfo->GetPrimaryParticleRange();
    for (PHG4TruthInfoContainer::ConstIterator iter = range.first;
         iter != range.second;
         ++iter)
    {
      PHG4Particle* primary = iter->second;

      if (primary->get_e() < _truth_e_threshold) continue;

      if (!_truth_trace_embed_flags.empty())
      {
        if (_truth_trace_embed_flags.find(trutheval->get_embed(primary)) ==
            _truth_trace_embed_flags.end()) continue;
      }

      float gparticleID = primary->get_track_id();
      float gflavor = primary->get_pid();

      PHG4Shower* shower = trutheval->get_primary_shower(primary);
      float gnhits = NAN;
      if (shower)
        gnhits = shower->get_nhits(trutheval->get_caloid());
      else
        gnhits = 0.0;
      float gpx = primary->get_px();
      float gpy = primary->get_py();
      float gpz = primary->get_pz();
      float ge = primary->get_e();

      float gpt = sqrt(gpx * gpx + gpy * gpy);
      float geta = NAN;
      if (gpt != 0.0) geta = asinh(gpz / gpt);
      float gphi = atan2(gpy, gpx);

      PHG4VtxPoint* vtx = trutheval->get_vertex(primary);
      float gvx = vtx->get_x();
      float gvy = vtx->get_y();
      float gvz = vtx->get_z();

      float gembed = trutheval->get_embed(primary);
      float gedep = trutheval->get_shower_energy_deposit(primary);

      RawCluster* cluster = clustereval->best_cluster_from(primary);

      float clusterID = NAN;
      float ntowers = NAN;
      float eta = NAN;
      float x = NAN;
      float y = NAN;
      float z = NAN;
      float phi = NAN;
      float e = NAN;

      float efromtruth = NAN;

      if (cluster)
      {
        clusterID = cluster->get_id();
        ntowers = cluster->getNTowers();
        x = cluster->get_x();
        y = cluster->get_y();
        z = cluster->get_z();
        phi = cluster->get_phi();
        e = cluster->get_energy();

        efromtruth = clustereval->get_energy_contribution(cluster, primary);

        // require vertex for cluster eta calculation
        if (vertexmap)
        {
          if (!vertexmap->empty())
          {
            GlobalVertex* vertex = (vertexmap->begin()->second);

            eta =
                RawClusterUtility::GetPseudorapidity(
                    *cluster,
                    CLHEP::Hep3Vector(vertex->get_x(), vertex->get_y(), vertex->get_z()));
          }
        }
      }

      float shower_data[] = {(float) _ievent,
                             gparticleID,
                             gflavor,
                             gnhits,
                             geta,
                             gphi,
                             ge,
                             gpt,
                             gvx,
                             gvy,
                             gvz,
                             gembed,
                             gedep,
                             clusterID,
                             ntowers,
                             eta,
                             x,
                             y,
                             z,
                             phi,
                             e,
                             efromtruth};

      _ntp_gshower->Fill(shower_data);
    }
  }

  //----------------------
  // fill the Tower NTuple
  //----------------------

  if (_do_tower_eval)
  {
    if (Verbosity() > 1) cout << "CaloEvaluator::filling tower ntuple..." << endl;

    string towernode = "TOWER_CALIB_" + _caloname;
    RawTowerContainer* towers = findNode::getClass<RawTowerContainer>(topNode, towernode.c_str());
    if (!towers)
    {
      cerr << PHWHERE << " ERROR: Can't find " << towernode << endl;
      exit(-1);
    }

    string towergeomnode = "TOWERGEOM_" + _caloname;
    RawTowerGeomContainer* towergeom = findNode::getClass<RawTowerGeomContainer>(topNode, towergeomnode.c_str());
    if (!towergeom)
    {
      cerr << PHWHERE << " ERROR: Can't find " << towergeomnode << endl;
      exit(-1);
    }

    RawTowerContainer::ConstRange begin_end = towers->getTowers();
    RawTowerContainer::ConstIterator rtiter;
    for (rtiter = begin_end.first; rtiter != begin_end.second; ++rtiter)
    {
      RawTower* tower = rtiter->second;

      if (tower->get_energy() < _reco_e_threshold) continue;

      RawTowerGeom* tower_geom = towergeom->get_tower_geometry(tower->get_id());
      if (!tower_geom)
      {
        cerr << PHWHERE << " ERROR: Can't find tower geometry for this tower hit: ";
        tower->identify();
        exit(-1);
      }

      //cout<<"Tower ID = "<<tower->get_id()<<" for bin(j,k)= "<<tower->get_bineta()<<","<<tower->get_binphi()<<endl; //Added by Barak
      const float towerid = tower->get_id();
      const float ieta = tower->get_bineta();
      const float iphi = tower->get_binphi();
      const float eta = tower_geom->get_eta();
      const float phi = tower_geom->get_phi();
      const float e = tower->get_energy();
      const float x = tower_geom->get_center_x();
      const float y = tower_geom->get_center_y();
      const float z = tower_geom->get_center_z();

      //Added by Barak
      _towerID_debug = tower->get_id();
      _ieta_debug = tower->get_bineta();
      _iphi_debug = tower->get_binphi();
      _eta_debug = tower_geom->get_eta();
      _phi_debug = tower_geom->get_phi();
      _e_debug = tower->get_energy();
      _x_debug = tower_geom->get_center_x();
      _y_debug = tower_geom->get_center_y();
      _z_debug = tower_geom->get_center_z();

      PHG4Particle* primary = towereval->max_truth_primary_particle_by_energy(tower);

      float gparticleID = NAN;
      float gflavor = NAN;
      float gnhits = NAN;
      float gpx = NAN;
      float gpy = NAN;
      float gpz = NAN;
      float ge = NAN;

      float gpt = NAN;
      float geta = NAN;
      float gphi = NAN;

      float gvx = NAN;
      float gvy = NAN;
      float gvz = NAN;

      float gembed = NAN;
      float gedep = NAN;

      float efromtruth = NAN;

      if (primary)
      {
        gparticleID = primary->get_track_id();
        gflavor = primary->get_pid();

        PHG4Shower* shower = trutheval->get_primary_shower(primary);
        if (shower)
          gnhits = shower->get_nhits(trutheval->get_caloid());
        else
          gnhits = 0.0;
        gpx = primary->get_px();
        gpy = primary->get_py();
        gpz = primary->get_pz();
        ge = primary->get_e();

        gpt = sqrt(gpx * gpx + gpy * gpy);
        if (gpt != 0.0) geta = asinh(gpz / gpt);
        gphi = atan2(gpy, gpx);

        PHG4VtxPoint* vtx = trutheval->get_vertex(primary);

        if (vtx)
        {
          gvx = vtx->get_x();
          gvy = vtx->get_y();
          gvz = vtx->get_z();
        }

        gembed = trutheval->get_embed(primary);
        gedep = trutheval->get_shower_energy_deposit(primary);

        efromtruth = towereval->get_energy_contribution(tower, primary);
      }

      float tower_data[] = {(float) _ievent,
                            towerid,
                            ieta,
                            iphi,
                            eta,
                            phi,
                            e,
                            x,
                            y,
                            z,
                            gparticleID,
                            gflavor,
                            gnhits,
                            geta,
                            gphi,
                            ge,
                            gpt,
                            gvx,
                            gvy,
                            gvz,
                            gembed,
                            gedep,
                            efromtruth};

      _ntp_tower->Fill(tower_data);
      _tower_debug->Fill();  //Added by Barak (see above for explanation)
    }
  }

  //------------------------
  // fill the Cluster NTuple
  //------------------------

  if (_do_cluster_eval)
  {
    if (Verbosity() > 1) cout << "CaloEvaluator::filling gcluster ntuple..." << endl;

    GlobalVertexMap* vertexmap = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");

    string clusternode = "CLUSTER_" + _caloname;
    RawClusterContainer* clusters = findNode::getClass<RawClusterContainer>(topNode, clusternode.c_str());
    if (!clusters)
    {
      cerr << PHWHERE << " ERROR: Can't find " << clusternode << endl;
      exit(-1);
    }

    // for every cluster

    for (const auto& iterator : clusters->getClustersMap())
    {
      RawCluster* cluster = iterator.second;

      //    for (unsigned int icluster = 0; icluster < clusters->size(); icluster++)
      //    {
      //      RawCluster* cluster = clusters->getCluster(icluster);

      if (cluster->get_energy() < _reco_e_threshold) continue;

      float clusterID = cluster->get_id();
      float ntowers = cluster->getNTowers();
      float x = cluster->get_x();
      float y = cluster->get_y();
      float z = cluster->get_z();
      float eta = NAN;
      float phi = cluster->get_phi();
      float e = cluster->get_energy();

      // require vertex for cluster eta calculation
      if (vertexmap)
      {
        if (!vertexmap->empty())
        {
          GlobalVertex* vertex = (vertexmap->begin()->second);

          eta =
              RawClusterUtility::GetPseudorapidity(
                  *cluster,
                  CLHEP::Hep3Vector(vertex->get_x(), vertex->get_y(), vertex->get_z()));
        }
      }

      PHG4Particle* primary = clustereval->max_truth_primary_particle_by_energy(cluster);

      float gparticleID = NAN;
      float gflavor = NAN;

      float gnhits = NAN;
      float gpx = NAN;
      float gpy = NAN;
      float gpz = NAN;
      float ge = NAN;

      float gpt = NAN;
      float geta = NAN;
      float gphi = NAN;

      float gvx = NAN;
      float gvy = NAN;
      float gvz = NAN;

      float gembed = NAN;
      float gedep = NAN;

      float efromtruth = NAN;

      if (primary)
      {
        gparticleID = primary->get_track_id();
        gflavor = primary->get_pid();

        PHG4Shower* shower = trutheval->get_primary_shower(primary);
        if (shower)
          gnhits = shower->get_nhits(trutheval->get_caloid());
        else
          gnhits = 0.0;
        gpx = primary->get_px();
        gpy = primary->get_py();
        gpz = primary->get_pz();
        ge = primary->get_e();

        gpt = sqrt(gpx * gpx + gpy * gpy);
        if (gpt != 0.0) geta = asinh(gpz / gpt);
        gphi = atan2(gpy, gpx);

        PHG4VtxPoint* vtx = trutheval->get_vertex(primary);

        if (vtx)
        {
          gvx = vtx->get_x();
          gvy = vtx->get_y();
          gvz = vtx->get_z();
        }

        gembed = trutheval->get_embed(primary);
        gedep = trutheval->get_shower_energy_deposit(primary);

        efromtruth = clustereval->get_energy_contribution(cluster,
                                                          primary);
      }

      float cluster_data[] = {(float) _ievent,
                              clusterID,
                              ntowers,
                              eta,
                              x,
                              y,
                              z,
                              phi,
                              e,
                              gparticleID,
                              gflavor,
                              gnhits,
                              geta,
                              gphi,
                              ge,
                              gpt,
                              gvx,
                              gvy,
                              gvz,
                              gembed,
                              gedep,
                              efromtruth};

      _ntp_cluster->Fill(cluster_data);
    }
  }

  return;
}
