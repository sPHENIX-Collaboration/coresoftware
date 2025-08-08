#include "KFParticle_truthAndDetTools.h"
#include "KFParticle_Tools.h"  // for KFParticle_Tools

#include <g4eval/SvtxEvalStack.h>   // for SvtxEvalStack
#include <g4eval/SvtxTrackEval.h>   // for SvtxTrackEval
#include <g4eval/SvtxTruthEval.h>   // for SvtxTruthEval
#include <g4eval/SvtxVertexEval.h>  // for SvtxVertexEval

#include <globalvertex/SvtxVertex.h>         // for SvtxVertex
#include <globalvertex/SvtxVertexMap.h>      // for SvtxVertexMap, SvtxVer...
#include <trackbase/InttDefs.h>              // for getLadderPhiId, getLad...
#include <trackbase/MvtxDefs.h>              // for getChipId, getStaveId
#include <trackbase/TpcDefs.h>               // for getSectorId, getSide
#include <trackbase/TrkrCluster.h>           // for TrkrCluster
#include <trackbase/TrkrClusterContainer.h>  // for TrkrClusterContainer
#include <trackbase/TrkrDefs.h>              // for getLayer, getTrkrId
#include <trackbase_historic/SvtxPHG4ParticleMap.h>
#include <trackbase_historic/SvtxTrack.h>     // for SvtxTrack, SvtxTrack::...
#include <trackbase_historic/SvtxTrackMap.h>  // for SvtxTrackMap, SvtxTrac...

#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMap.h>

#include <g4main/PHG4Particle.h>            // for PHG4Particle
#include <g4main/PHG4TruthInfoContainer.h>  // for PHG4TruthInfoContainer
#include <g4main/PHG4VtxPoint.h>            // for PHG4VtxPoint

#include <phhepmc/PHHepMCGenEvent.h>     // for PHHepMCGenEvent
#include <phhepmc/PHHepMCGenEventMap.h>  // for PHHepMCGenEventMap

#include <phool/getClass.h>  // for getClass

#include <KFParticle.h>  // for KFParticle

#include <TSystem.h>     // for gSystem->Exit()
#include <TTree.h>       // for TTree

#include <HepMC/GenEvent.h>       // for GenEvent::particle_con...
#include <HepMC/GenParticle.h>    // for GenParticle
#include <HepMC/GenVertex.h>      // for GenVertex::particle_it...
#include <HepMC/IteratorRange.h>  // for parents
#include <HepMC/SimpleVector.h>   // for FourVector

#include <algorithm>  // for max, find
#include <cmath>      // for pow, sqrt
#include <cstdlib>    // for NULL, abs
#include <iostream>   // for operator<<, endl, basi...
#include <iterator>   // for end, begin
#include <map>        // for _Rb_tree_iterator, map
#include <memory>     // for allocator_traits<>::va...
#include <utility>    // for pair

std::map<std::string, int> Use =
    {
        {"MVTX", 1},
        {"INTT", 1},
        {"TPC", 1},
        {"TPOT", 1},
        {"EMCAL", 0},
        {"OHCAL", 0},
        {"IHCAL", 0}};

SvtxTrack *KFParticle_truthAndDetTools::getTrack(unsigned int track_id, SvtxTrackMap *trackmap)
{
  SvtxTrack *matched_track = nullptr;

  for (auto &iter : *trackmap)
  {
    if (iter.first == track_id)
    {
      matched_track = iter.second;
    }
  }

  return matched_track;
}

GlobalVertex *KFParticle_truthAndDetTools::getVertex(unsigned int vertex_id, GlobalVertexMap *vertexmap)
{
  GlobalVertex *matched_vertex = vertexmap->get(vertex_id);

  return matched_vertex;
}

PHG4Particle *KFParticle_truthAndDetTools::getTruthTrack(SvtxTrack *thisTrack, PHCompositeNode *topNode)
{
  /*
   * There are two methods for getting the truth rack from the reco track
   * 1. (recommended) Use the reco -> truth tables (requires SvtxPHG4ParticleMap). Introduced Summer of 2022
   * 2. Get truth track via nClusters. Older method and will work with older DSTs
   */

  PHG4Particle *particle = nullptr;

  SvtxPHG4ParticleMap *dst_reco_truth_map = findNode::getClass<SvtxPHG4ParticleMap>(topNode, "SvtxPHG4ParticleMap");
  if (dst_reco_truth_map)
  {
    m_truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
    if (!m_truthinfo)
    {
      std::cout << "KFParticle truth matching: G4TruthInfo does not exist" << std::endl;
    }

    std::map<float, std::set<int>> truth_set = dst_reco_truth_map->get(thisTrack->get_id());
    if (!truth_set.empty())
    {
      std::pair<float, std::set<int>> best_weight = *truth_set.rbegin();
      int best_truth_id = *best_weight.second.rbegin();
      particle = m_truthinfo->GetParticle(best_truth_id);
    }
  }
  else
  {
    std::cout << __FILE__ << ": SvtxPHG4ParticleMap not found, reverting to max_truth_particle_by_nclusters()" << std::endl;

    if (!m_svtx_evalstack)
    {
      m_svtx_evalstack = new SvtxEvalStack(topNode);
      // clustereval = m_svtx_evalstack->get_cluster_eval();
      // hiteval = m_svtx_evalstack->get_hit_eval();
      trackeval = m_svtx_evalstack->get_track_eval();
      trutheval = m_svtx_evalstack->get_truth_eval();
      vertexeval = m_svtx_evalstack->get_vertex_eval();
    }

    m_svtx_evalstack->next_event(topNode);

    particle = trackeval->max_truth_particle_by_nclusters(thisTrack);
  }
  return particle;
}

void KFParticle_truthAndDetTools::initializeTruthBranches(TTree *m_tree, int daughter_id, const std::string &daughter_number, bool m_constrain_to_vertex_truthMatch)
{
  m_tree->Branch((daughter_number + "_true_ID").c_str(), &m_true_daughter_id[daughter_id], (daughter_number + "_true_ID/I").c_str());
  if (m_constrain_to_vertex_truthMatch)
  {
    m_tree->Branch((daughter_number + "_true_IP").c_str(), &m_true_daughter_ip[daughter_id], (daughter_number + "_true_IP/F").c_str());
    m_tree->Branch((daughter_number + "_true_IP_xy").c_str(), &m_true_daughter_ip_xy[daughter_id], (daughter_number + "_true_IP_xy/F").c_str());
  }
  m_tree->Branch((daughter_number + "_true_px").c_str(), &m_true_daughter_px[daughter_id], (daughter_number + "_true_px/F").c_str());
  m_tree->Branch((daughter_number + "_true_py").c_str(), &m_true_daughter_py[daughter_id], (daughter_number + "_true_py/F").c_str());
  m_tree->Branch((daughter_number + "_true_pz").c_str(), &m_true_daughter_pz[daughter_id], (daughter_number + "_true_pz/F").c_str());
  m_tree->Branch((daughter_number + "_true_p").c_str(), &m_true_daughter_p[daughter_id], (daughter_number + "_true_p/F").c_str());
  m_tree->Branch((daughter_number + "_true_pT").c_str(), &m_true_daughter_pt[daughter_id], (daughter_number + "_true_pT/F").c_str());
  m_tree->Branch((daughter_number + "_true_EV_x").c_str(), &m_true_daughter_vertex_x[daughter_id], (daughter_number + "_true_EV_x/F").c_str());
  m_tree->Branch((daughter_number + "_true_EV_y").c_str(), &m_true_daughter_vertex_y[daughter_id], (daughter_number + "_true_EV_y/F").c_str());
  m_tree->Branch((daughter_number + "_true_EV_z").c_str(), &m_true_daughter_vertex_z[daughter_id], (daughter_number + "_true_EV_z/F").c_str());
  if (m_constrain_to_vertex_truthMatch)
  {
    m_tree->Branch((daughter_number + "_true_PV_x").c_str(), &m_true_daughter_pv_x[daughter_id], (daughter_number + "_true_PV_x/F").c_str());
    m_tree->Branch((daughter_number + "_true_PV_y").c_str(), &m_true_daughter_pv_y[daughter_id], (daughter_number + "_true_PV_y/F").c_str());
    m_tree->Branch((daughter_number + "_true_PV_z").c_str(), &m_true_daughter_pv_z[daughter_id], (daughter_number + "_true_PV_z/F").c_str());
  }
  m_tree->Branch((daughter_number + "_true_track_history_PDG_ID").c_str(), &m_true_daughter_track_history_PDG_ID[daughter_id]);
  m_tree->Branch((daughter_number + "_true_track_history_PDG_mass").c_str(), &m_true_daughter_track_history_PDG_mass[daughter_id]);
  m_tree->Branch((daughter_number + "_true_track_history_px").c_str(), &m_true_daughter_track_history_px[daughter_id]);
  m_tree->Branch((daughter_number + "_true_track_history_py").c_str(), &m_true_daughter_track_history_py[daughter_id]);
  m_tree->Branch((daughter_number + "_true_track_history_pz").c_str(), &m_true_daughter_track_history_pz[daughter_id]);
  m_tree->Branch((daughter_number + "_true_track_history_pE").c_str(), &m_true_daughter_track_history_pE[daughter_id]);
  m_tree->Branch((daughter_number + "_true_track_history_pT").c_str(), &m_true_daughter_track_history_pT[daughter_id]);
}

void KFParticle_truthAndDetTools::fillTruthBranch(PHCompositeNode *topNode, TTree * /*m_tree*/, const KFParticle &daughter, int daughter_id, const KFParticle &kfvertex, bool m_constrain_to_vertex_truthMatch)
{
  float true_px;
  float true_py;
  float true_pz;
  float true_p;
  float true_pt;

  dst_trackmap = findNode::getClass<SvtxTrackMap>(topNode, m_trk_map_node_name_nTuple);
  if (!dst_trackmap)
  {
    std::cout << "KFParticle truth matching: " << m_trk_map_node_name_nTuple << " does not exist" << std::endl;
  }

  if (m_use_mbd_vertex_truth)
  {
    dst_mbdvertexmap = findNode::getClass<MbdVertexMap>(topNode, "MbdVertexMap");
    if (!dst_mbdvertexmap)
    {
      std::cout << "KFParticle truth matching: MbdVertexMap does not exist" << std::endl;
    }
  }
  else
  {
    dst_vertexmap = findNode::getClass<SvtxVertexMap>(topNode, m_vtx_map_node_name_nTuple);
    if (!dst_vertexmap)
    {
      std::cout << "KFParticle truth matching: " << m_vtx_map_node_name_nTuple << " does not exist" << std::endl;
    }
  }

  auto *globalvertexmap = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");
  if (!globalvertexmap)
  {
    std::cout << "KFParticle truth matching: GlobalVertexMap does not exist" << std::endl;
  }

  track = getTrack(daughter.Id(), dst_trackmap);
  g4particle = getTruthTrack(track, topNode);

  bool isParticleValid = g4particle == nullptr ? false : true;

  true_px = isParticleValid ? (Float_t) g4particle->get_px() : 0.;
  true_py = isParticleValid ? (Float_t) g4particle->get_py() : 0.;
  true_pz = isParticleValid ? (Float_t) g4particle->get_pz() : 0.;
  true_p = sqrt(pow(true_px, 2) + pow(true_py, 2) + pow(true_pz, 2));
  true_pt = sqrt(pow(true_px, 2) + pow(true_py, 2));

  m_true_daughter_px[daughter_id] = true_px;
  m_true_daughter_py[daughter_id] = true_py;
  m_true_daughter_pz[daughter_id] = true_pz;
  m_true_daughter_p[daughter_id] = true_p;
  m_true_daughter_pt[daughter_id] = true_pt;
  m_true_daughter_id[daughter_id] = isParticleValid ? g4particle->get_pid() : 0;

  if (!m_svtx_evalstack)
  {
    m_svtx_evalstack = new SvtxEvalStack(topNode);
    // clustereval = m_svtx_evalstack->get_cluster_eval();
    // hiteval = m_svtx_evalstack->get_hit_eval();
    trackeval = m_svtx_evalstack->get_track_eval();
    trutheval = m_svtx_evalstack->get_truth_eval();
    vertexeval = m_svtx_evalstack->get_vertex_eval();
  }

  if (isParticleValid)
  {
    g4vertex_point = trutheval->get_vertex(g4particle);
  }

  m_true_daughter_vertex_x[daughter_id] = isParticleValid ? g4vertex_point->get_x() : 0.;
  m_true_daughter_vertex_y[daughter_id] = isParticleValid ? g4vertex_point->get_y() : 0.;
  m_true_daughter_vertex_z[daughter_id] = isParticleValid ? g4vertex_point->get_z() : 0.;

  if (m_constrain_to_vertex_truthMatch)
  {
    // Calculate true DCA
    GlobalVertex *recoVertex = getVertex(kfvertex.Id(), globalvertexmap);
    GlobalVertex::VTXTYPE whichVtx = m_use_mbd_vertex_truth ? GlobalVertex::MBD : GlobalVertex::SVTX;
    auto svtxviter = recoVertex->find_vertexes(whichVtx);

    // check that it contains a track vertex
    if (svtxviter == recoVertex->end_vertexes())
    {
      std::string vtxType = m_use_mbd_vertex_truth ? "MBD" : "silicon";
      std::cout << "Have a global vertex with no " << vtxType << " vertex... shouldn't happen in KFParticle_truthAndDetTools::fillTruthBranch..." << std::endl;
    }

    auto svtxvertexvector = svtxviter->second;
    MbdVertex *mbdvertex = nullptr;
    SvtxVertex *svtxvertex = nullptr;

    for (auto &vertex_iter : svtxvertexvector)
    {
      if (m_use_mbd_vertex_truth)
      {
        mbdvertex = dst_mbdvertexmap->find(vertex_iter->get_id())->second;
      }
      else
      {
        svtxvertex = dst_vertexmap->find(vertex_iter->get_id())->second;
      }
    }

    PHG4VtxPoint *truePoint = nullptr;
    if (m_use_mbd_vertex_truth)
    {
      std::set<PHG4VtxPoint *> truePointSet = vertexeval->all_truth_points(mbdvertex);
      truePoint = *truePointSet.begin();
    }
    else
    {
      truePoint = vertexeval->max_truth_point_by_ntracks(svtxvertex);
    }

    if (truePoint == nullptr && isParticleValid)
    {
      // PHG4Particle *g4mother = m_truthinfo->GetParticle(g4particle->get_parent_id());
      PHG4Particle *g4mother = m_truthinfo->GetPrimaryParticle(g4particle->get_parent_id());
      if (!g4mother)
      {
        std::cout << "KFParticle truth matching: True mother not found!\n";
        std::cout << "Your truth track DCA will be measured wrt a reconstructed vertex!" << std::endl;
        truePoint = nullptr;
      }
      else
      {
        truePoint = m_truthinfo->GetVtx(g4mother->get_vtx_id());  // Note, this may not be the PV for a decay with tertiaries
      }
    }

    KFParticle trueKFParticleVertex;

    float f_vertexParameters[6] = {0};

    if (truePoint == nullptr)
    {
      std::cout << "KFParticle truth matching: This event has no PHG4VtxPoint information!\n";
      std::cout << "Your truth track DCA will be measured wrt a reconstructed vertex!" << std::endl;

      f_vertexParameters[0] = recoVertex->get_x();
      f_vertexParameters[1] = recoVertex->get_y();
      f_vertexParameters[2] = recoVertex->get_z();
    }
    else
    {
      f_vertexParameters[0] = truePoint->get_x();
      f_vertexParameters[1] = truePoint->get_y();
      f_vertexParameters[2] = truePoint->get_z();
    }

    float f_vertexCovariance[21] = {0};

    trueKFParticleVertex.Create(f_vertexParameters, f_vertexCovariance, 0, -1);

    KFParticle trueKFParticle;

    float f_trackParameters[6] = {m_true_daughter_vertex_x[daughter_id],
                                  m_true_daughter_vertex_y[daughter_id],
                                  m_true_daughter_vertex_z[daughter_id],
                                  true_px,
                                  true_py,
                                  true_pz};

    float f_trackCovariance[21] = {0};

    trueKFParticle.Create(f_trackParameters, f_trackCovariance, 1, -1);

    m_true_daughter_ip[daughter_id] = trueKFParticle.GetDistanceFromVertex(trueKFParticleVertex);
    m_true_daughter_ip_xy[daughter_id] = trueKFParticle.GetDistanceFromVertexXY(trueKFParticleVertex);

    m_true_daughter_pv_x[daughter_id] = truePoint == nullptr ? -99. : truePoint->get_x();
    m_true_daughter_pv_y[daughter_id] = truePoint == nullptr ? -99. : truePoint->get_y();
    m_true_daughter_pv_z[daughter_id] = truePoint == nullptr ? -99. : truePoint->get_z();
  }
}

void KFParticle_truthAndDetTools::fillGeant4Branch(PHG4Particle *particle, int daughter_id)
{
  Float_t pT = sqrt(pow(particle->get_px(), 2) + pow(particle->get_py(), 2));

  m_true_daughter_track_history_PDG_ID[daughter_id].push_back(particle->get_pid());
  m_true_daughter_track_history_PDG_mass[daughter_id].push_back(0);
  m_true_daughter_track_history_px[daughter_id].push_back((Float_t) particle->get_px());
  m_true_daughter_track_history_py[daughter_id].push_back((Float_t) particle->get_py());
  m_true_daughter_track_history_pz[daughter_id].push_back((Float_t) particle->get_pz());
  m_true_daughter_track_history_pE[daughter_id].push_back((Float_t) particle->get_e());
  m_true_daughter_track_history_pT[daughter_id].push_back(pT);
}

void KFParticle_truthAndDetTools::fillHepMCBranch(HepMC::GenParticle *particle, int daughter_id)
{
  const HepMC::FourVector &myFourVector = particle->momentum();

  m_true_daughter_track_history_PDG_ID[daughter_id].push_back(particle->pdg_id());
  m_true_daughter_track_history_PDG_mass[daughter_id].push_back((Float_t) particle->generatedMass());
  m_true_daughter_track_history_px[daughter_id].push_back((Float_t) myFourVector.px());
  m_true_daughter_track_history_py[daughter_id].push_back((Float_t) myFourVector.py());
  m_true_daughter_track_history_pz[daughter_id].push_back((Float_t) myFourVector.pz());
  m_true_daughter_track_history_pE[daughter_id].push_back((Float_t) myFourVector.e());
  m_true_daughter_track_history_pT[daughter_id].push_back((Float_t) myFourVector.perp());
}

int KFParticle_truthAndDetTools::getHepMCInfo(PHCompositeNode *topNode, TTree * /*m_tree*/, const KFParticle &daughter, int daughter_id)
{
  // Make dummy particle for null pointers and missing nodes
  HepMC::GenParticle *dummyParticle = new HepMC::GenParticle();
  HepMC::FourVector dummyFourVector(0, 0, 0, 0);
  dummyParticle->set_momentum(dummyFourVector);
  dummyParticle->set_pdg_id(0);
  dummyParticle->set_generated_mass(0.);

  dst_trackmap = findNode::getClass<SvtxTrackMap>(topNode, m_trk_map_node_name_nTuple);
  if (!dst_trackmap)
  {
    std::cout << "KFParticle truth matching: " << m_trk_map_node_name_nTuple << " does not exist" << std::endl;
    gSystem->Exit(1);
    exit(1);
  }

  track = getTrack(daughter.Id(), dst_trackmap);
  g4particle = getTruthTrack(track, topNode);

  bool isParticleValid = g4particle == nullptr ? false : true;

  if (!isParticleValid)
  {
    std::cout << "KFParticle truth matching: this track is a ghost" << std::endl;
    fillHepMCBranch(dummyParticle, daughter_id);
    return 0;
  }

  m_geneventmap = findNode::getClass<PHHepMCGenEventMap>(topNode, "PHHepMCGenEventMap");
  if (!m_geneventmap)
  {
    std::cout << "KFParticle truth matching: Missing node PHHepMCGenEventMap" << std::endl;
    std::cout << "You will have no mother information" << std::endl;
    fillHepMCBranch(dummyParticle, daughter_id);
    return 0;
  }

  m_genevt = m_geneventmap->get(1);
  if (!m_genevt)
  {
    std::cout << "KFParticle truth matching: Missing node PHHepMCGenEvent" << std::endl;
    std::cout << "You will have no mother information" << std::endl;
    fillHepMCBranch(dummyParticle, daughter_id);
    return 0;
  }

  // Start by looking for our particle in the Geant record
  // Any decay that Geant4 handles will not be in the HepMC record
  // This can happen if you limit the decay volume in the generator
  if (g4particle->get_parent_id() != 0)
  {
    m_truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
    if (!m_truthinfo)
    {
      std::cout << "KFParticle truth matching: G4TruthInfo does not exist" << std::endl;
    }
    while (g4particle->get_parent_id() != 0)
    {
      g4particle = m_truthinfo->GetParticle(g4particle->get_parent_id());
      fillGeant4Branch(g4particle, daughter_id);
    }
  }

  HepMC::GenEvent *theEvent = m_genevt->getEvent();
  HepMC::GenParticle *prevParticle = nullptr;

  // int forbiddenPDGIDs[] = {21, 22};  //Stop tracing history when we reach quarks, gluons and photons
  int forbiddenPDGIDs[] = {0};  // 20230921 - Request made to have gluon information to see about gluon-splitting

  for (HepMC::GenEvent::particle_const_iterator p = theEvent->particles_begin(); p != theEvent->particles_end(); ++p)
  {
    if (((*p)->barcode() == g4particle->get_barcode()))
    {
      prevParticle = *p;
      while (!prevParticle->is_beam())
      {
        bool breakOut = false;
        for (HepMC::GenVertex::particle_iterator mother = prevParticle->production_vertex()->particles_begin(HepMC::parents);
             mother != prevParticle->production_vertex()->particles_end(HepMC::parents); ++mother)
        {
          if (std::find(std::begin(forbiddenPDGIDs), std::end(forbiddenPDGIDs),
                        abs((*mother)->pdg_id())) != std::end(forbiddenPDGIDs))
          {
            breakOut = true;
            break;
          }

          fillHepMCBranch((*mother), daughter_id);
          prevParticle = *mother;
        }
        if (breakOut)
        {
          break;
        }
      }
    }
  }

  return 0;
}  // End of function

void KFParticle_truthAndDetTools::initializeCaloBranches(TTree *m_tree, int daughter_id, const std::string &daughter_number)
{
  m_tree->Branch((daughter_number + "_EMCAL_DeltaPhi").c_str(), &detector_emcal_deltaphi[daughter_id], (daughter_number + "_EMCAL_DeltaPhi/F").c_str());
  m_tree->Branch((daughter_number + "_EMCAL_DeltaEta").c_str(), &detector_emcal_deltaeta[daughter_id], (daughter_number + "_EMCAL_DeltaEta/F").c_str());
  m_tree->Branch((daughter_number + "_EMCAL_DeltaZ").c_str(), &detector_emcal_deltaz[daughter_id], (daughter_number + "_EMCAL_DeltaZ/F").c_str());
  m_tree->Branch((daughter_number + "_EMCAL_energy_3x3").c_str(), &detector_emcal_energy_3x3[daughter_id], (daughter_number + "_EMCAL_energy_3x3/F").c_str());
  m_tree->Branch((daughter_number + "_EMCAL_energy_5x5").c_str(), &detector_emcal_energy_5x5[daughter_id], (daughter_number + "_EMCAL_energy_5x5/F").c_str());
  m_tree->Branch((daughter_number + "_EMCAL_energy_cluster").c_str(), &detector_emcal_cluster_energy[daughter_id], (daughter_number + "_EMCAL_energy_cluster/F").c_str());
  m_tree->Branch((daughter_number + "_IHCAL_DeltaPhi").c_str(), &detector_ihcal_deltaphi[daughter_id], (daughter_number + "_IHCAL_DeltaPhi/F").c_str());
  m_tree->Branch((daughter_number + "_IHCAL_DeltaEta").c_str(), &detector_ihcal_deltaeta[daughter_id], (daughter_number + "_IHCAL_DeltaEta/F").c_str());
  m_tree->Branch((daughter_number + "_IHCAL_energy_3x3").c_str(), &detector_ihcal_energy_3x3[daughter_id], (daughter_number + "_IHCAL_energy_3x3/F").c_str());
  m_tree->Branch((daughter_number + "_IHCAL_energy_5x5").c_str(), &detector_ihcal_energy_5x5[daughter_id], (daughter_number + "_IHCAL_energy_5x5/F").c_str());
  m_tree->Branch((daughter_number + "_IHCAL_energy_cluster").c_str(), &detector_ihcal_cluster_energy[daughter_id], (daughter_number + "_IHCAL_energy_cluster/F").c_str());
  m_tree->Branch((daughter_number + "_OHCAL_DeltaPhi").c_str(), &detector_ohcal_deltaphi[daughter_id], (daughter_number + "_OHCAL_DeltaEta/F").c_str());
  m_tree->Branch((daughter_number + "_OHCAL_DeltaEta").c_str(), &detector_ohcal_deltaeta[daughter_id], (daughter_number + "_OHCAL_DeltaEta/F").c_str());
  m_tree->Branch((daughter_number + "_OHCAL_energy_3x3").c_str(), &detector_ohcal_energy_3x3[daughter_id], (daughter_number + "_OHCAL_energy_3x3/F").c_str());
  m_tree->Branch((daughter_number + "_OHCAL_energy_5x5").c_str(), &detector_ohcal_energy_5x5[daughter_id], (daughter_number + "_OHCAL_energy_5x5/F").c_str());
  m_tree->Branch((daughter_number + "_OHCAL_energy_cluster").c_str(), &detector_ohcal_cluster_energy[daughter_id], (daughter_number + "_OHCAL_energy_cluster/F").c_str());
}

/*
The following function matches tracks to calo clusters. As of 7/1/2025, this only extends to the EMCal. HCal matching is in development.

To run EMCal matching, DST_CALO files must be read into the Fun4All server. 
*/
void KFParticle_truthAndDetTools::fillCaloBranch(PHCompositeNode *topNode,
                                                 TTree * /*m_tree*/, const KFParticle &daughter, int daughter_id, bool &isTrackEMCalmatch)
{
  dst_trackmap = findNode::getClass<SvtxTrackMap>(topNode, m_trk_map_node_name_nTuple);
  if (!dst_trackmap)
  {
    std::cout << "KFParticle truth matching: " << m_trk_map_node_name_nTuple << " does not exist" << std::endl;
    gSystem->Exit(1);
    exit(1);
  }

  track = getTrack(daughter.Id(), dst_trackmap);

  if (!clustersEM)
  {
    clustersEM = findNode::getClass<RawClusterContainer>(topNode, "CLUSTERINFO_CEMC");
    if (!clustersEM)
    {
      std::cout << __FILE__ << "::" << __func__ << " : FATAL ERROR, cannot find cluster container " << "CLUSTERINFO_CEMC" << std::endl;
    }
  }
  // if (!EMCalGeo)
  // {
  //   EMCalGeo = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_CEMC");
  //   if (!EMCalGeo)
  //   {
  //     std::cout << __FILE__ << "::" << __func__ << " : FATAL ERROR, cannot find cluster container " << "TOWERGEOM_CEMC" << std::endl;
  //     // return Fun4AllReturnCodes::ABORTEVENT;
  //   }
  // }
  // if(!_towersEM)
  // {
  //   _towersEM = findNode::getClass<RawTowerContainer>(topNode, "TOWER_CALIB_CEMC");
  //   if(!_towersEM)
  //   {
  //     std::cout << __FILE__ << "::" << __func__ << " : FATAL ERROR, cannot find cluster container " << "TOWER_CALIB_CEMC" << std::endl;
  //     //return Fun4AllReturnCodes::ABORTEVENT;
  //   }
  // }
  if (!clustersIH)
  {
    clustersIH = findNode::getClass<RawClusterContainer>(topNode, "CLUSTER_HCALIN");
    // if (!clustersIH)
    // {
    //   std::cout << __FILE__ << "::" << __func__ << " : FATAL ERROR, cannot find cluster container " << "CLUSTER_HCALIN" << std::endl;
    //   //return Fun4AllReturnCodes::ABORTEVENT;
    // }
  }
  if (!IHCalGeo)
  {
    IHCalGeo = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALIN");
    // if(!IHCalGeo)
    // {
    //   std::cout << __FILE__ << "::" << __func__ << " : FATAL ERROR, cannot find cluster container " << "TOWERGEOM_HCALIN" << std::endl;
    //   //return Fun4AllReturnCodes::ABORTEVENT;
    // }
  }
  // if(!_towersIH)
  // {
  //   _towersIH = findNode::getClass<RawTowerContainer>(topNode, "TOWER_CALIB_HCALIN");
  //   if(!_towersIH)
  //   {
  //     std::cout << __FILE__ << "::" << __func__ << " : FATAL ERROR, cannot find cluster container " << "TOWER_CALIB_HCALIN" << std::endl;
  //     //return Fun4AllReturnCodes::ABORTEVENT;
  //   }
  // }
  if (!clustersOH)
  {
    clustersOH = findNode::getClass<RawClusterContainer>(topNode, "CLUSTER_HCALOUT");
    // if (!clustersOH)
    // {
    //   std::cout << __FILE__ << "::" << __func__ << " : FATAL ERROR, cannot find cluster container " << "CLUSTER_HCALOUT" << std::endl;
    //   //return Fun4AllReturnCodes::ABORTEVENT;
    // }
  }
  if (!OHCalGeo)
  {
    OHCalGeo = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALOUT");
    // if(!OHCalGeo)
    // {
    //   std::cout << __FILE__ << "::" << __func__ << " : FATAL ERROR, cannot find cluster container " << "TOWERGEOM_HCALOUT" << std::endl;
    //   //return Fun4AllReturnCodes::ABORTEVENT;
    // }
  }
  // if(!_towersOH)
  // {
  //   _towersOH = findNode::getClass<RawTowerContainer>(topNode, "TOWER_CALIB_HCALOUT");
  //   if(!_towersOH)
  //   {
  //     std::cout << __FILE__ << "::" << __func__ << " : FATAL ERROR, cannot find cluster container " << "TOWER_CALIB_HCALOUT" << std::endl;
  //     //return Fun4AllReturnCodes::ABORTEVENT;
  //   }
  // }

  // Radii for track projections
  double caloRadiusEMCal;
  caloRadiusEMCal = m_emcal_radius_user;  //Use function set_emcal_radius_user(float set_variable) to set this in your Fun4All macro
  // double caloRadiusIHCal;
  // double caloRadiusOHCal;
  // caloRadiusEMCal = 100.70;
  // caloRadiusEMCal = EMCalGeo->get_radius(); //This requires DST_CALOFITTING 
  // caloRadiusOHCal = OHCalGeo->get_radius();
  // caloRadiusIHCal = IHCalGeo->get_radius();

  // Create variables and containers etc.
  bool is_match;
  int index;
  float radius_scale;
  RawCluster *cluster;

  // EMCAL******************************************************
  //  std::cout << "Starting EMCAL-track matching!" << std::endl;

  SvtxTrackState *thisState = nullptr;
  thisState = track->get_state(caloRadiusEMCal);
  // thisState->identify();
  // std::cout << "size states " << (size_t)track->size_states() << std::endl;
  // for (auto state_iter = track->begin_states();
  //       state_iter != track->end_states();
  //       ++state_iter)
  //   {
  //     SvtxTrackState* tstate = state_iter->second;
  //     tstate->identify();
  //     std::cout<< "radius " << sqrt((tstate->get_x())*(tstate->get_x()) + (tstate->get_y())*(tstate->get_y())) << std::endl;
  //   }

  // assert(thisState);
  float _track_phi_emc = std::numeric_limits<float>::quiet_NaN();
  float _track_eta_emc = std::numeric_limits<float>::quiet_NaN();
  float _track_x_emc = std::numeric_limits<float>::quiet_NaN();
  float _track_y_emc = std::numeric_limits<float>::quiet_NaN();
  float _track_z_emc = std::numeric_limits<float>::quiet_NaN();

  // EMCal variables and vectors
  float _emcal_phi = std::numeric_limits<float>::quiet_NaN();
  float _emcal_eta = std::numeric_limits<float>::quiet_NaN();
  float _emcal_x = std::numeric_limits<float>::quiet_NaN();
  float _emcal_y = std::numeric_limits<float>::quiet_NaN();
  float _emcal_z = std::numeric_limits<float>::quiet_NaN();
  radius_scale = std::numeric_limits<float>::quiet_NaN();
  float _emcal_3x3 = std::numeric_limits<float>::quiet_NaN();
  float _emcal_5x5 = std::numeric_limits<float>::quiet_NaN();
  float _emcal_clusE = std::numeric_limits<float>::quiet_NaN();
  std::vector<float> v_emcal_phi;
  std::vector<float> v_emcal_eta;
  std::vector<float> v_emcal_x;
  std::vector<float> v_emcal_y;
  std::vector<float> v_emcal_z;
  std::vector<float> v_emcal_dphi;
  std::vector<float> v_emcal_deta;
  std::vector<float> v_emcal_dr;
  std::vector<float> v_emcal_dz;
  std::vector<float> v_emcal_3x3;
  std::vector<float> v_emcal_5x5;
  std::vector<float> v_emcal_clusE;

  // Set variables for matching
  is_match = false;
  index = -1;
  // int ijk = 0; // nothing is being done with this variable in the end

  //clustersEM->identify();

  if (thisState != nullptr)
  {
    _track_x_emc = thisState->get_x();
    _track_y_emc = thisState->get_y();
    _track_z_emc = thisState->get_z();
    _track_phi_emc = std::atan2(_track_y_emc, _track_x_emc);
    _track_eta_emc = std::asinh(_track_z_emc / std::sqrt((_track_x_emc * _track_x_emc) + (_track_y_emc * _track_y_emc)));

    // Create objects, containers, iterators for clusters
    cluster = nullptr;
    RawClusterContainer::Range begin_end_EMC = clustersEM->getClusters();
    RawClusterContainer::Iterator clusIter_EMC;

    // Loop over the EMCal clusters
    for (clusIter_EMC = begin_end_EMC.first; clusIter_EMC != begin_end_EMC.second; ++clusIter_EMC)
    {
      // Minimum energy cut
      cluster = clusIter_EMC->second;
      float cluster_energy = cluster->get_energy();

      if (cluster_energy < m_emcal_e_low_cut)
      {
        // ijk++;
        continue;
      }

      // Get cluster information
      _emcal_x = cluster->get_x();
      _emcal_y = cluster->get_y();
      _emcal_z = cluster->get_z();
      radius_scale = m_emcal_radius_user / std::sqrt((_emcal_x * _emcal_x) + (_emcal_y * _emcal_y));
      _emcal_x *= radius_scale;
      _emcal_y *= radius_scale;
      _emcal_z *= radius_scale;
      _emcal_phi = std::atan2(_emcal_y, _emcal_y);
      _emcal_eta = std::asinh(_emcal_z / std::sqrt((_emcal_x * _emcal_x) + (_emcal_y * _emcal_y)));
      // _emcal_3x3 = get_e3x3(cluster, _towersEM, 0); //0 for emcal
      // _emcal_5x5 = get_e5x5(cluster, _towersEM, 0); //0 for emcal
      _emcal_3x3 = std::numeric_limits<float>::quiet_NaN();
      _emcal_5x5 = std::numeric_limits<float>::quiet_NaN();
      _emcal_clusE = cluster_energy;

      // Variables to determine potential matches
      float dphi = PiRange(_track_phi_emc - _emcal_phi);
      float dz = _track_z_emc - _emcal_z;
      float deta = _track_eta_emc - _emcal_eta;
      float tmparg = caloRadiusEMCal * dphi;
      float dr = std::sqrt((tmparg * tmparg) + (dz * dz));  // sqrt((R*dphi)^2 + (dz)^2
      // float dr = sqrt((dphi*dphi + deta*deta)); //previous version

      // Requirements for a possible match
      if (dz > m_dz_cut_high || dz < m_dz_cut_low)
      {
        continue;
      }
      if (dphi > m_dphi_cut_high || dphi < m_dphi_cut_low)
      {
        continue;
      }

      // std::cout << "**********DELTA INFORMATION************" << std::endl;
      // std::cout << "dphi = " << dphi << std::endl;
      // std::cout << "deta = " << deta << std::endl;
      // std::cout << "dz =   " << dz <<   std::endl;
      // std::cout << "dr =   " << dr <<   std::endl;
      // std::cout << "****************************************" << std::endl;
      // Add potential match's information to vectors
      v_emcal_phi.push_back(_emcal_phi);
      v_emcal_eta.push_back(_emcal_eta);
      v_emcal_x.push_back(_emcal_x);
      v_emcal_y.push_back(_emcal_y);
      v_emcal_z.push_back(_emcal_z);
      v_emcal_dphi.push_back(dphi);
      v_emcal_dz.push_back(dz);
      v_emcal_deta.push_back(deta);
      v_emcal_dr.push_back(dr);
      v_emcal_3x3.push_back(_emcal_3x3);
      v_emcal_5x5.push_back(_emcal_5x5);
      v_emcal_clusE.push_back(_emcal_clusE);
      is_match = true;
      // ijk++;
    }

    // Find the closest match from all potential matches
    if (is_match == true)
    {
      float tmp = 99999;
      for (long unsigned int i = 0; i < v_emcal_dr.size(); i++)
      {
        if (v_emcal_dr[i] < tmp)
        {
          index = i;
          tmp = v_emcal_dr[i];
        }
      }
    }
  }

  /*
  // Print out statements
  if (index != -1)
  {
    std::cout << "matched tracks!!!" << std::endl;
    std::cout << "emcal x = " << v_emcal_x[index] << " , y = " << v_emcal_y[index] << " , z = " << v_emcal_z[index] << " , phi = " << v_emcal_phi[index] << " , eta = " << v_emcal_eta[index] << std::endl;
    std::cout << "track projected x = " << _track_x_emc << " , y = " << _track_y_emc << " , z = " << _track_z_emc << " , phi = " << _track_phi_emc << " , eta = " << _track_eta_emc << std::endl;
    std::cout << "track px = " << track->get_px() << " , py = " << track->get_py() << " , pz = " << track->get_pz() << " , pt = " << track->get_pt() << " , p = " << track->get_p() << " , charge = " << track->get_charge() << std::endl;
  }
  */

  // Save values to the branches!
  if (index == -1)
  {
    detector_emcal_deltaphi[daughter_id] = std::numeric_limits<float>::quiet_NaN();
    detector_emcal_deltaeta[daughter_id] = std::numeric_limits<float>::quiet_NaN();
    detector_emcal_deltaeta[daughter_id] = std::numeric_limits<float>::quiet_NaN();
    detector_emcal_energy_3x3[daughter_id] = std::numeric_limits<float>::quiet_NaN();
    detector_emcal_energy_5x5[daughter_id] = std::numeric_limits<float>::quiet_NaN();
    detector_emcal_cluster_energy[daughter_id] = std::numeric_limits<float>::quiet_NaN();
    isTrackEMCalmatch = false;
  }
  else
  {
    detector_emcal_deltaphi[daughter_id] = v_emcal_dphi[index];
    detector_emcal_deltaeta[daughter_id] = v_emcal_deta[index];
    detector_emcal_deltaz[daughter_id] = v_emcal_dz[index];
    detector_emcal_energy_3x3[daughter_id] = std::numeric_limits<float>::quiet_NaN();
    detector_emcal_energy_5x5[daughter_id] = std::numeric_limits<float>::quiet_NaN();
    detector_emcal_cluster_energy[daughter_id] = v_emcal_clusE[index];
    isTrackEMCalmatch = true;
  }

  // HCAL*******************************************************

  // INNER
  /*
  std::cout << "Starting IHCAL-track matching!" << std::endl;

  // Track projection
  thisState = nullptr;
  thisState = track->get_state(caloRadiusIHCal);
  float _track_phi_ihc = std::numeric_limits<float>::quiet_NaN();
  float _track_eta_ihc = std::numeric_limits<float>::quiet_NaN();
  float _track_x_ihc = std::numeric_limits<float>::quiet_NaN();
  float _track_y_ihc = std::numeric_limits<float>::quiet_NaN();
  float _track_z_ihc = std::numeric_limits<float>::quiet_NaN();

  // IHCal variables and vectors
  float _ihcal_phi = std::numeric_limits<float>::quiet_NaN();
  float _ihcal_eta = std::numeric_limits<float>::quiet_NaN();
  float _ihcal_x = std::numeric_limits<float>::quiet_NaN();
  float _ihcal_y = std::numeric_limits<float>::quiet_NaN();
  float _ihcal_z = std::numeric_limits<float>::quiet_NaN();
  float _ihcal_3x3 = std::numeric_limits<float>::quiet_NaN();
  //float _ihcal_5x5 = std::numeric_limits<float>::quiet_NaN();
  float _ihcal_clusE = std::numeric_limits<float>::quiet_NaN();
  radius_scale = std::numeric_limits<float>::quiet_NaN();
  std::vector<float> v_ihcal_phi;
  std::vector<float> v_ihcal_eta;
  std::vector<float> v_ihcal_x;
  std::vector<float> v_ihcal_y;
  std::vector<float> v_ihcal_z;
  std::vector<float> v_ihcal_dphi;
  std::vector<float> v_ihcal_deta;
  std::vector<float> v_ihcal_dr;
  std::vector<float> v_ihcal_3x3;
  //std::vector<float> v_ihcal_5x5;
  std::vector<float> v_ihcal_clusE;

  if(thisState != nullptr){

  // Reset variables for matching
  is_match = false;
  index = -1;

  _track_phi_ihc = atan2(thisState->get_y(), thisState->get_x());
  _track_eta_ihc = asinh(thisState->get_z()/sqrt(thisState->get_x()*thisState->get_x() + thisState->get_y()*thisState->get_y()));
  _track_x_ihc = thisState->get_x();
  _track_y_ihc = thisState->get_y();
  _track_z_ihc = thisState->get_z();


  // Create objects, containers, iterators for clusters
  cluster = nullptr;
  RawClusterContainer::Range begin_end_IHC = clustersIH->getClusters();
  RawClusterContainer::Iterator clusIter_IHC;

  // Loop over the IHCal clusters
  for (clusIter_IHC = begin_end_IHC.first; clusIter_IHC != begin_end_IHC.second; ++clusIter_IHC)
  {

  // Minimum energy cut
  cluster = clusIter_IHC->second;
  if(cluster->get_energy() < m_ihcal_e_low_cut)
  {
  continue;
  }

  std::cout << "Found a cluster about threshhold energy!" << std::endl;

  // Get cluster information
  _ihcal_phi = atan2(cluster->get_y(), cluster->get_x());
  _ihcal_eta = asinh(cluster->get_z()/sqrt(cluster->get_x()*cluster->get_x() + cluster->get_y()*cluster->get_y()));
  _ihcal_x = cluster->get_x();
  _ihcal_y = cluster->get_y();
  radius_scale = m_ihcal_radius_user / sqrt(_ihcal_x*_ihcal_x+_ihcal_y*_ihcal_y);
  _ihcal_z = radius_scale*cluster->get_z();
  // _ihcal_3x3 = get_e3x3(cluster, _towersIH, 0); //1 for ihcal
  // _ihcal_5x5 = get_e5x5(cluster, _towersIH, 0); //1 for ihcal
  _ihcal_3x3 = std::numeric_limits<float>::quiet_NaN();
  _ihcal_5x5 = std::numeric_limits<float>::quiet_NaN();
  _ihcal_clusE = cluster->get_energy();


  // Variables to determine potential matches
  float dphi = abs(PiRange(_track_phi_ihc - _ihcal_phi));
  float dz = abs(_track_z_ihc - _ihcal_z);
  float deta = abs(_ihcal_eta - _track_eta_ihc);
  float dr = sqrt((dphi*dphi + deta*deta));

  // Requirements for a possible match
  if(dphi<m_dphi_cut && dz<m_dz_cut)
  {
  //Add potential match's information to vectors
  v_ihcal_phi.push_back(_ihcal_phi);
  v_ihcal_eta.push_back(_ihcal_eta);
  v_ihcal_x.push_back(_ihcal_x);
  v_ihcal_y.push_back(_ihcal_y);
  v_ihcal_z.push_back(_ihcal_z);
  v_ihcal_dphi.push_back(dphi);
  v_ihcal_deta.push_back(deta);
  v_ihcal_dr.push_back(dr);
  v_ihcal_3x3.push_back(_ihcal_3x3);
  v_ihcal_clusE.push_back(_ihcal_clusE);

  is_match = true;
  }
  }

  // Find the closest match from all potential matches
  if (is_match == true)
  {
  float tmp = 99999;
  for(long unsigned int i = 0; i < v_ihcal_dr.size(); i++){
  if(v_ihcal_dr[i] < tmp){
  index = i;
  tmp = v_ihcal_dr[i];
  }
  }
  }
  }

  // Print out statihents
  if(index != -1){
  std::cout<<"matched tracks!!!"<<std::endl;
  std::cout<<"ihcal x = "<<v_ihcal_x[index]<<" , y = "<<v_ihcal_y[index]<<" , z = "<<v_ihcal_z[index]<<" , phi = "<<v_ihcal_phi[index]<<" , eta = "<<v_ihcal_eta[index]<<std::endl;
  std::cout<<"track projected x = "<<_track_x_ihc<<" , y = "<<_track_y_ihc<<" , z = "<<_track_z_ihc<<" , phi = "<<_track_phi_ihc<<" , eta = "<<_track_eta_ihc<<std::endl;
  std::cout<<"track px = "<<track->get_px()<<" , py = "<<track->get_py()<<" , pz = "<<track->get_pz()<<" , pt = "<<track->get_pt()<<" , p = "<<track->get_p()<<" , charge = "<<track->get_charge()<<std::endl;
  }
  */

  // Save values to the branches!
  // if(index == -1){
  detector_ihcal_deltaphi[daughter_id] = std::numeric_limits<float>::quiet_NaN();
  detector_ihcal_deltaeta[daughter_id] = std::numeric_limits<float>::quiet_NaN();
  detector_ihcal_energy_3x3[daughter_id] = std::numeric_limits<float>::quiet_NaN();
  detector_ihcal_energy_5x5[daughter_id] = std::numeric_limits<float>::quiet_NaN();
  detector_ihcal_cluster_energy[daughter_id] = std::numeric_limits<float>::quiet_NaN();
  // }
  // else{
  // detector_ihcal_deltaphi[daughter_id] = v_ihcal_dphi[index];
  // detector_ihcal_deltaeta[daughter_id] = v_ihcal_deta[index];
  // detector_ihcal_energy_3x3[daughter_id] = std::numeric_limits<float>::quiet_NaN();
  // detector_ihcal_energy_5x5[daughter_id] = std::numeric_limits<float>::quiet_NaN();
  // detector_ihcal_cluster_energy[daughter_id] = v_ihcal_clusE[index];
  // }

  // std::cout << "IHCAL CLLUSTERS MATCHED TO TRACK" << std::endl;

  // OUTER

  // std::cout << "Starting OHCAL-track matching!" << std::endl;

  /*

  // Track projection
  thisState = nullptr;
  thisState = track->get_state(caloRadiusOHCal);
  float _track_phi_ohc = std::numeric_limits<float>::quiet_NaN();
  float _track_eta_ohc = std::numeric_limits<float>::quiet_NaN();
  float _track_x_ohc = std::numeric_limits<float>::quiet_NaN();
  float _track_y_ohc = std::numeric_limits<float>::quiet_NaN();
  float _track_z_ohc = std::numeric_limits<float>::quiet_NaN();

  // OHCal variables and vectors
  float _ohcal_phi = std::numeric_limits<float>::quiet_NaN();
  float _ohcal_eta = std::numeric_limits<float>::quiet_NaN();
  float _ohcal_x = std::numeric_limits<float>::quiet_NaN();
  float _ohcal_y = std::numeric_limits<float>::quiet_NaN();
  float _ohcal_z = std::numeric_limits<float>::quiet_NaN();
  float _ohcal_3x3 = std::numeric_limits<float>::quiet_NaN();
  //float _ohcal_5x5 = std::numeric_limits<float>::quiet_NaN();
  float _ohcal_clusE = std::numeric_limits<float>::quiet_NaN();
  radius_scale = std::numeric_limits<float>::quiet_NaN();
  std::vector<float> v_ohcal_phi;
  std::vector<float> v_ohcal_eta;
  std::vector<float> v_ohcal_x;
  std::vector<float> v_ohcal_y;
  std::vector<float> v_ohcal_z;
  std::vector<float> v_ohcal_dphi;
  std::vector<float> v_ohcal_deta;
  std::vector<float> v_ohcal_dr;
  std::vector<float> v_ohcal_3x3;
  //std::vector<float> v_ohcal_5x5;
  std::vector<float> v_ohcal_clusE;

  // Reset variables for matching
  is_match = false;
  index = -1;

  if(thisState != nullptr)
  {
  _track_phi_ohc = atan2(thisState->get_y(), thisState->get_x());
  _track_eta_ohc = asinh(thisState->get_z()/sqrt(thisState->get_x()*thisState->get_x() + thisState->get_y()*thisState->get_y()));
  _track_x_ohc = thisState->get_x();
  _track_y_ohc = thisState->get_y();
  _track_z_ohc = thisState->get_z();


  // Create objects, containers, iterators for clusters
  cluster = nullptr;
  RawClusterContainer::Range begin_end_OHC = clustersOH->getClusters();
  RawClusterContainer::Iterator clusIter_OHC;

  // Loop over the OHCal clusters
  for (clusIter_OHC = begin_end_OHC.first; clusIter_OHC != begin_end_OHC.second; ++clusIter_OHC)
  {

  // Minimum energy cut
  cluster = clusIter_OHC->second;
  if(cluster->get_energy() < m_ohcal_e_low_cut)
  {
  continue;
  }

  // Get cluster information
  _ohcal_phi = atan2(cluster->get_y(), cluster->get_x());
  _ohcal_eta = asinh(cluster->get_z()/sqrt(cluster->get_x()*cluster->get_x() + cluster->get_y()*cluster->get_y()));
  _ohcal_x = cluster->get_x();
  _ohcal_y = cluster->get_y();
  radius_scale = m_ohcal_radius_user / sqrt(_ohcal_x*_ohcal_x+_ohcal_y*_ohcal_y);
  _ohcal_z = radius_scale*cluster->get_z();
  // _ohcal_3x3 = get_e3x3(cluster, _towersOH, 0); //2 for ohcal
  // _ohcal_5x5 = get_e5x5(cluster, _towersOH, 0); //2 for ohcal
  _ohcal_3x3 = std::numeric_limits<float>::quiet_NaN();
  _ohcal_5x5 = std::numeric_limits<float>::quiet_NaN();
  _ohcal_clusE = cluster->get_energy();

  // Variables to determine potential matches
  float dphi = abs(PiRange(_track_phi_ohc - _ohcal_phi));
  float dz = abs(_track_z_ohc - _ohcal_z);
  float deta = abs(_ohcal_eta - _track_eta_ohc);
  float dr = sqrt((dphi*dphi + deta*deta));

  // Requirohents for a possible match
  if(dphi<m_dphi_cut && dz<m_dz_cut)
  {
  //Add potential match's information to vectors
  v_ohcal_phi.push_back(_ohcal_phi);
  v_ohcal_eta.push_back(_ohcal_eta);
  v_ohcal_x.push_back(_ohcal_x);
  v_ohcal_y.push_back(_ohcal_y);
  v_ohcal_z.push_back(_ohcal_z);
  v_ohcal_dphi.push_back(dphi);
  v_ohcal_deta.push_back(deta);
  v_ohcal_dr.push_back(dr);
  v_ohcal_3x3.push_back(_ohcal_3x3);
  v_ohcal_clusE.push_back(_ohcal_clusE);

  is_match = true;
  }
  }


  // Find the closest match from all potential matches
  if (is_match == true)
  {
  float tmp = 99999;
  for(long unsigned int i = 0; i < v_ohcal_dr.size(); i++){
  if(v_ohcal_dr[i] < tmp){
  index = i;
  tmp = v_ohcal_dr[i];
  }
  }
  }
  }

  std::cout << "match identified" << std::endl;

  // Print out statements
  if(index != -1){
  std::cout<<"matched tracks!!!"<<std::endl;
  std::cout<<"ohcal x = "<<v_ohcal_x[index]<<" , y = "<<v_ohcal_y[index]<<" , z = "<<v_ohcal_z[index]<<" , phi = "<<v_ohcal_phi[index]<<" , eta = "<<v_ohcal_eta[index]<<std::endl;
  std::cout<<"track projected x = "<<_track_x_ohc<<" , y = "<<_track_y_ohc<<" , z = "<<_track_z_ohc<<" , phi = "<<_track_phi_ohc<<" , eta = "<<_track_eta_ohc<<std::endl;
  std::cout<<"track px = "<<track->get_px()<<" , py = "<<track->get_py()<<" , pz = "<<track->get_pz()<<" , pt = "<<track->get_pt()<<" , p = "<<track->get_p()<<" , charge = "<<track->get_charge()<<std::endl;
  }
  */

  // Save values to the branches!
  // if(index == -1){
  detector_ohcal_deltaphi[daughter_id] = std::numeric_limits<float>::quiet_NaN();
  detector_ohcal_deltaeta[daughter_id] = std::numeric_limits<float>::quiet_NaN();
  detector_ohcal_energy_3x3[daughter_id] = std::numeric_limits<float>::quiet_NaN();
  detector_ohcal_energy_5x5[daughter_id] = std::numeric_limits<float>::quiet_NaN();
  detector_ohcal_cluster_energy[daughter_id] = std::numeric_limits<float>::quiet_NaN();
  // }
  // else{
  // detector_ohcal_deltaphi[daughter_id] = v_ohcal_dphi[index];
  // detector_ohcal_deltaeta[daughter_id] = v_ohcal_deta[index];
  // detector_ohcal_energy_3x3[daughter_id] = std::numeric_limits<float>::quiet_NaN();
  // detector_ohcal_energy_5x5[daughter_id] = std::numeric_limits<float>::quiet_NaN();
  // detector_ohcal_cluster_energy[daughter_id] = v_ohcal_clusE[index];
  // }
  // std::cout << "OHCAL CLLUSTERS MATCHED TO TRACK" << std::endl;
}

void KFParticle_truthAndDetTools::initializeDetectorBranches(TTree *m_tree, int daughter_id, const std::string &daughter_number)
{
  m_tree->Branch((daughter_number + "_residual_x").c_str(), &residual_x[daughter_id]);
  m_tree->Branch((daughter_number + "_residual_y").c_str(), &residual_y[daughter_id]);
  m_tree->Branch((daughter_number + "_residual_z").c_str(), &residual_z[daughter_id]);
  m_tree->Branch((daughter_number + "_layer").c_str(), &detector_layer[daughter_id]);

  for (auto const &subdetector : Use)
  {
    if (subdetector.second)
    {
      initializeSubDetectorBranches(m_tree, subdetector.first, daughter_id, daughter_number);
    }
  }
}

void KFParticle_truthAndDetTools::initializeSubDetectorBranches(TTree *m_tree, const std::string &detectorName, int daughter_id, const std::string &daughter_number)
{
  if (detectorName == "MVTX")
  {
    m_tree->Branch((daughter_number + "_" + detectorName + "_staveID").c_str(), &mvtx_staveID[daughter_id]);
    m_tree->Branch((daughter_number + "_" + detectorName + "_chipID").c_str(), &mvtx_chipID[daughter_id]);
    m_tree->Branch((daughter_number + "_" + detectorName + "_nHits").c_str(), &detector_nHits_MVTX[daughter_id]);
    m_tree->Branch((daughter_number + "_" + detectorName + "_nStates").c_str(), &detector_nStates_MVTX[daughter_id]);
  }
  if (detectorName == "INTT")
  {
    m_tree->Branch((daughter_number + "_" + detectorName + "_ladderZID").c_str(), &intt_ladderZID[daughter_id]);
    m_tree->Branch((daughter_number + "_" + detectorName + "_ladderPhiID").c_str(), &intt_ladderPhiID[daughter_id]);
    m_tree->Branch((daughter_number + "_" + detectorName + "_nHits").c_str(), &detector_nHits_INTT[daughter_id]);
    m_tree->Branch((daughter_number + "_" + detectorName + "_nStates").c_str(), &detector_nStates_INTT[daughter_id]);
  }
  if (detectorName == "TPC")
  {
    m_tree->Branch((daughter_number + "_" + detectorName + "_sectorID").c_str(), &tpc_sectorID[daughter_id]);
    m_tree->Branch((daughter_number + "_" + detectorName + "_side").c_str(), &tpc_side[daughter_id]);
    m_tree->Branch((daughter_number + "_" + detectorName + "_nHits").c_str(), &detector_nHits_TPC[daughter_id]);
    m_tree->Branch((daughter_number + "_" + detectorName + "_nStates").c_str(), &detector_nStates_TPC[daughter_id]);
  }
  if (detectorName == "TPOT")
  {
    m_tree->Branch((daughter_number + "_" + detectorName + "_nHits").c_str(), &detector_nHits_TPOT[daughter_id]);
    m_tree->Branch((daughter_number + "_" + detectorName + "_nStates").c_str(), &detector_nStates_TPOT[daughter_id]);
  }
}

void KFParticle_truthAndDetTools::fillDetectorBranch(PHCompositeNode *topNode,
                                                     TTree * /*m_tree*/, const KFParticle &daughter, int daughter_id)
{
  dst_trackmap = findNode::getClass<SvtxTrackMap>(topNode, m_trk_map_node_name_nTuple);
  if (!dst_trackmap)
  {
    std::cout << "KFParticle truth matching: " << m_trk_map_node_name_nTuple << " does not exist" << std::endl;
  }

  std::string geoName = "ActsGeometry";
  geometry = findNode::getClass<ActsGeometry>(topNode, geoName);
  if (!geometry)
  {
    std::cout << "KFParticle detector info: " << geoName << " does not exist" << std::endl;
  }

  dst_clustermap = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  if (!dst_clustermap)
  {
    std::cout << "KFParticle detector info: TRKR_CLUSTER does not exist" << std::endl;
  }

  track = getTrack(daughter.Id(), dst_trackmap);
  detector_nHits_MVTX[daughter_id] = 0;
  detector_nHits_INTT[daughter_id] = 0;
  detector_nHits_TPC[daughter_id] = 0;
  detector_nHits_TPOT[daughter_id] = 0;

  TrackSeed *silseed = track->get_silicon_seed();
  TrackSeed *tpcseed = track->get_tpc_seed();

  if (silseed)
  {
    for (auto cluster_iter = silseed->begin_cluster_keys(); cluster_iter != silseed->end_cluster_keys(); ++cluster_iter)
    {
      const auto &cluster_key = *cluster_iter;
      const auto trackerID = TrkrDefs::getTrkrId(cluster_key);

      detector_layer[daughter_id].push_back(TrkrDefs::getLayer(cluster_key));

      unsigned int staveId;
      unsigned int chipId;
      unsigned int ladderZId;
      unsigned int ladderPhiId;
      unsigned int sectorId;
      unsigned int side;
      staveId = chipId = ladderZId = ladderPhiId = sectorId = side = std::numeric_limits<unsigned int>::quiet_NaN();

      if (Use["MVTX"] && trackerID == TrkrDefs::mvtxId)
      {
        staveId = MvtxDefs::getStaveId(cluster_key);
        chipId = MvtxDefs::getChipId(cluster_key);
        ++detector_nHits_MVTX[daughter_id];
      }
      else if (Use["INTT"] && trackerID == TrkrDefs::inttId)
      {
        ladderZId = InttDefs::getLadderZId(cluster_key);
        ladderPhiId = InttDefs::getLadderPhiId(cluster_key);
        ++detector_nHits_INTT[daughter_id];
      }

      mvtx_staveID[daughter_id].push_back(staveId);
      mvtx_chipID[daughter_id].push_back(chipId);
      intt_ladderZID[daughter_id].push_back(ladderZId);
      intt_ladderPhiID[daughter_id].push_back(ladderPhiId);
      tpc_sectorID[daughter_id].push_back(sectorId);
      tpc_side[daughter_id].push_back(side);
    }
  }

  if (tpcseed)
  {
    for (auto cluster_iter = tpcseed->begin_cluster_keys(); cluster_iter != tpcseed->end_cluster_keys(); ++cluster_iter)
    {
      const auto &cluster_key = *cluster_iter;
      const auto trackerID = TrkrDefs::getTrkrId(cluster_key);

      detector_layer[daughter_id].push_back(TrkrDefs::getLayer(cluster_key));

      unsigned int staveId;
      unsigned int chipId;
      unsigned int ladderZId;
      unsigned int ladderPhiId;
      unsigned int sectorId;
      unsigned int side;
      staveId = chipId = ladderZId = ladderPhiId = sectorId = side = std::numeric_limits<unsigned int>::quiet_NaN();

      if (Use["TPC"] && trackerID == TrkrDefs::tpcId)
      {
        sectorId = TpcDefs::getSectorId(cluster_key);
        side = TpcDefs::getSide(cluster_key);
        ++detector_nHits_TPC[daughter_id];
      }
      else if (Use["TPOT"] && trackerID == TrkrDefs::micromegasId)
      {
        ++detector_nHits_TPOT[daughter_id];
      }

      mvtx_staveID[daughter_id].push_back(staveId);
      mvtx_chipID[daughter_id].push_back(chipId);
      intt_ladderZID[daughter_id].push_back(ladderZId);
      intt_ladderPhiID[daughter_id].push_back(ladderPhiId);
      tpc_sectorID[daughter_id].push_back(sectorId);
      tpc_side[daughter_id].push_back(side);
    }
  }

  for (auto state_iter = track->begin_states();
       state_iter != track->end_states();
       ++state_iter)
  {
    SvtxTrackState *tstate = state_iter->second;
    if (tstate->get_pathlength() != 0)  // The first track state is an extrapolation so has no cluster
    {
      auto stateckey = tstate->get_cluskey();
      TrkrCluster *cluster = dst_clustermap->findCluster(stateckey);
      if (!cluster)
      {
	// do not have associated cluster, could be track states projected to calo system
        continue;
      }
      auto global = geometry->getGlobalPosition(stateckey, cluster);

      residual_x[daughter_id].push_back(global.x() - tstate->get_x());
      residual_y[daughter_id].push_back(global.y() - tstate->get_y());
      residual_z[daughter_id].push_back(global.z() - tstate->get_z());

      uint8_t id = TrkrDefs::getTrkrId(stateckey);

      switch (id)
      {
      case TrkrDefs::mvtxId:
        ++detector_nStates_MVTX[daughter_id];
        break;
      case TrkrDefs::inttId:
        ++detector_nStates_INTT[daughter_id];
        break;
      case TrkrDefs::tpcId:
        ++detector_nStates_TPC[daughter_id];
        break;
      case TrkrDefs::micromegasId:
        ++detector_nStates_TPOT[daughter_id];
        break;
      default:
        //std::cout << "Cluster key doesnt match a tracking system, could be related with projected track state to calorimeter system" << std::endl;
        break;
      }
    }
  }
}

int KFParticle_truthAndDetTools::getPVID(PHCompositeNode *topNode, const KFParticle &kfpvertex)
{
  if (m_dont_use_global_vertex_truth)
  {
    if (m_use_mbd_vertex_truth)
    {
      dst_mbdvertexmap = findNode::getClass<MbdVertexMap>(topNode, "MbdVertexMap");
      if (dst_mbdvertexmap)
      {
        MbdVertex *m_dst_vertex = dst_mbdvertexmap->get(kfpvertex.Id());
        return m_dst_vertex->get_beam_crossing();
      }

      std::cout << "KFParticle vertex matching: " << m_vtx_map_node_name_nTuple << " does not exist" << std::endl;
    }
    else
    {
      dst_vertexmap = findNode::getClass<SvtxVertexMap>(topNode, m_vtx_map_node_name_nTuple);
      if (dst_vertexmap)
      {
        SvtxVertex *m_dst_vertex = dst_vertexmap->get(kfpvertex.Id());
        return m_dst_vertex->get_beam_crossing();
      }

      std::cout << "KFParticle vertex matching: " << m_vtx_map_node_name_nTuple << " does not exist" << std::endl;
    }
  }
  else
  {
    auto globalvertexmap = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");
    if (!globalvertexmap)
    {
      return -100;
    }

    GlobalVertex *gvertex = globalvertexmap->get(kfpvertex.Id());
    return gvertex->get_beam_crossing();
  }

  return -100;
}

void KFParticle_truthAndDetTools::allPVInfo(PHCompositeNode *topNode,
                                            TTree * /*m_tree*/,
                                            const KFParticle &motherParticle,
                                            std::vector<KFParticle> daughters,
                                            std::vector<KFParticle> intermediates)
{
  KFParticle_Tools kfpTupleTools;
  kfpTupleTools.set_dont_use_global_vertex(m_dont_use_global_vertex_truth);
  std::vector<KFParticle> primaryVertices = kfpTupleTools.makeAllPrimaryVertices(topNode, m_vtx_map_node_name_nTuple);

  for (auto &primaryVertice : primaryVertices)
  {
    allPV_x.push_back(primaryVertice.GetX());
    allPV_y.push_back(primaryVertice.GetY());
    allPV_z.push_back(primaryVertice.GetZ());

    allPV_mother_IP.push_back(motherParticle.GetDistanceFromVertex(primaryVertice));
    allPV_mother_IPchi2.push_back(motherParticle.GetDeviationFromVertex(primaryVertice));

    for (unsigned int j = 0; j < daughters.size(); ++j)
    {
      allPV_daughter_IP[j].push_back(daughters[j].GetDistanceFromVertex(primaryVertice));
      allPV_daughter_IPchi2[j].push_back(daughters[j].GetDeviationFromVertex(primaryVertice));
    }

    for (unsigned int j = 0; j < intermediates.size(); ++j)
    {
      allPV_intermediates_IP[j].push_back(intermediates[j].GetDistanceFromVertex(primaryVertice));
      allPV_intermediates_IPchi2[j].push_back(intermediates[j].GetDeviationFromVertex(primaryVertice));
    }
  }
}

void KFParticle_truthAndDetTools::clearVectors()
{
  for (int i = 0; i < m_num_tracks_nTuple; ++i)
  {
    // Truth vectors
    m_true_daughter_track_history_PDG_ID[i].clear();
    m_true_daughter_track_history_PDG_mass[i].clear();
    m_true_daughter_track_history_px[i].clear();
    m_true_daughter_track_history_py[i].clear();
    m_true_daughter_track_history_pz[i].clear();
    m_true_daughter_track_history_pE[i].clear();
    m_true_daughter_track_history_pT[i].clear();

    // Detector vectors
    residual_x[i].clear();
    residual_y[i].clear();
    residual_z[i].clear();
    detector_layer[i].clear();
    mvtx_staveID[i].clear();
    mvtx_chipID[i].clear();
    intt_ladderZID[i].clear();
    intt_ladderPhiID[i].clear();
    tpc_sectorID[i].clear();
    tpc_side[i].clear();

    detector_nStates_MVTX[i] = 0;
    detector_nStates_INTT[i] = 0;
    detector_nStates_TPC[i] = 0;
    detector_nStates_TPOT[i] = 0;

    // PV vectors
    allPV_daughter_IP[i].clear();
    allPV_daughter_IPchi2[i].clear();
  }

  allPV_x.clear();
  allPV_y.clear();
  allPV_z.clear();
  allPV_z.clear();

  allPV_mother_IP.clear();
  allPV_mother_IPchi2.clear();

  for (int i = 0; i < m_num_intermediate_states_nTuple; ++i)
  {
    allPV_intermediates_IP[i].clear();
    allPV_intermediates_IPchi2[i].clear();
  }
}
