#include "KFParticle_truthAndDetTools.h"

#include "KFParticle_Tools.h"  // for KFParticle_Tools

#include <g4eval/SvtxEvalStack.h>   // for SvtxEvalStack
#include <g4eval/SvtxTrackEval.h>   // for SvtxTrackEval
#include <g4eval/SvtxTruthEval.h>   // for SvtxTruthEval
#include <g4eval/SvtxVertexEval.h>  // for SvtxVertexEval

#include <trackbase/InttDefs.h>                // for getLadderPhiId, getLad...
#include <trackbase/MvtxDefs.h>                // for getChipId, getStaveId
#include <trackbase/TpcDefs.h>                 // for getSectorId, getSide
#include <trackbase/TrkrCluster.h>             // for TrkrCluster
#include <trackbase/TrkrClusterContainer.h>    // for TrkrClusterContainer
#include <trackbase/TrkrDefs.h>                // for getLayer, getTrkrId
#include <trackbase_historic/SvtxPHG4ParticleMap_v1.h>
#include <trackbase_historic/SvtxTrack.h>      // for SvtxTrack, SvtxTrack::...
#include <trackbase_historic/SvtxTrackMap.h>   // for SvtxTrackMap, SvtxTrac...
#include <trackbase_historic/SvtxVertex.h>     // for SvtxVertex
#include <trackbase_historic/SvtxVertexMap.h>  // for SvtxVertexMap, SvtxVer...

#include <g4main/PHG4Particle.h>            // for PHG4Particle
#include <g4main/PHG4TruthInfoContainer.h>  // for PHG4TruthInfoContainer
#include <g4main/PHG4VtxPoint.h>            // for PHG4VtxPoint

#include <phhepmc/PHHepMCGenEvent.h>     // for PHHepMCGenEvent
#include <phhepmc/PHHepMCGenEventMap.h>  // for PHHepMCGenEventMap

#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/getClass.h>        // for getClass

#include <KFParticle.h>  // for KFParticle
#include <TString.h>     // for TString, operator+
#include <TTree.h>       // for TTree

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <HepMC/GenEvent.h>       // for GenEvent::particle_con...
#include <HepMC/GenVertex.h>      // for GenVertex::particle_it...
#pragma GCC diagnostic pop

#include <HepMC/GenParticle.h>    // for GenParticle
#include <HepMC/IteratorRange.h>  // for parents
#include <HepMC/SimpleVector.h>   // for FourVector

#include <math.h>     // for pow, sqrt
#include <stdlib.h>   // for NULL, abs
#include <algorithm>  // for max, find
#include <iostream>   // for operator<<, endl, basi...
#include <iterator>   // for end, begin
#include <map>        // for _Rb_tree_iterator, map
#include <memory>     // for allocator_traits<>::va...
#include <utility>    // for pair

class PHNode;

std::map<std::string, int> Use =
    {
        {"MVTX", 1},
        {"INTT", 1},
        {"TPC", 1},
        {"EMCAL", 0},
        {"OHCAL", 0},
        {"IHCAL", 0}};

KFParticle_truthAndDetTools::KFParticle_truthAndDetTools()
  : m_svtx_evalstack(nullptr)
{
}  //Constructor

KFParticle_truthAndDetTools::~KFParticle_truthAndDetTools() {}  //Destructor

SvtxTrack *KFParticle_truthAndDetTools::getTrack(unsigned int track_id, SvtxTrackMap *trackmap)
{
  SvtxTrack *matched_track = NULL;

  for (SvtxTrackMap::Iter iter = trackmap->begin();
       iter != trackmap->end();
       ++iter)
  {
    if (iter->first == track_id) matched_track = iter->second;
  }

  return matched_track;
}

SvtxVertex *KFParticle_truthAndDetTools::getVertex(unsigned int vertex_id, SvtxVertexMap *vertexmap)
{
  SvtxVertex *matched_vertex = vertexmap->get(vertex_id);

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

  PHNodeIterator nodeIter(topNode);
  PHNode *findNode = dynamic_cast<PHNode *>(nodeIter.findFirst("SvtxPHG4ParticleMap"));
  if (findNode)
  {
    findNode = dynamic_cast<PHNode *>(nodeIter.findFirst("G4TruthInfo"));
    if (findNode)
    {
      m_truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
    }
    else
    {
      std::cout << "KFParticle truth matching: G4TruthInfo does not exist" << std::endl;
    }

    SvtxPHG4ParticleMap_v1 *dst_reco_truth_map = findNode::getClass<SvtxPHG4ParticleMap_v1>(topNode, "SvtxPHG4ParticleMap");

    std::map<float, std::set<int>> truth_set = dst_reco_truth_map->get(thisTrack->get_id()); 
    const auto& best_weight = truth_set.rbegin();
    int best_truth_id =  *best_weight->second.rbegin();
    particle = m_truthinfo->GetParticle(best_truth_id);
  }
  else
  {
    std::cout << __FILE__ << ": SvtxPHG4ParticleMap not found, reverting to max_truth_particle_by_nclusters()" << std::endl;
  
    if (!m_svtx_evalstack)
    {
      m_svtx_evalstack = new SvtxEvalStack(topNode);
      //clustereval = m_svtx_evalstack->get_cluster_eval();
      //hiteval = m_svtx_evalstack->get_hit_eval();
      trackeval = m_svtx_evalstack->get_track_eval();
      trutheval = m_svtx_evalstack->get_truth_eval();
      vertexeval = m_svtx_evalstack->get_vertex_eval();
    }
  
    m_svtx_evalstack->next_event(topNode);
  
    particle = trackeval->max_truth_particle_by_nclusters(thisTrack);
  }
  return particle;
}

void KFParticle_truthAndDetTools::initializeTruthBranches(TTree *m_tree, int daughter_id, std::string daughter_number, bool m_constrain_to_vertex_truthMatch)
{
  m_tree->Branch(TString(daughter_number) + "_true_ID", &m_true_daughter_id[daughter_id], TString(daughter_number) + "_true_ID/I");
  if (m_constrain_to_vertex_truthMatch)
  {
    m_tree->Branch(TString(daughter_number) + "_true_IP", &m_true_daughter_ip[daughter_id], TString(daughter_number) + "_true_IP/F");
    m_tree->Branch(TString(daughter_number) + "_true_IP_xy", &m_true_daughter_ip_xy[daughter_id], TString(daughter_number) + "_true_IP_xy/F");
  }
  m_tree->Branch(TString(daughter_number) + "_true_px", &m_true_daughter_px[daughter_id], TString(daughter_number) + "_true_px/F");
  m_tree->Branch(TString(daughter_number) + "_true_py", &m_true_daughter_py[daughter_id], TString(daughter_number) + "_true_py/F");
  m_tree->Branch(TString(daughter_number) + "_true_pz", &m_true_daughter_pz[daughter_id], TString(daughter_number) + "_true_pz/F");
  m_tree->Branch(TString(daughter_number) + "_true_p", &m_true_daughter_p[daughter_id], TString(daughter_number) + "_true_p/F");
  m_tree->Branch(TString(daughter_number) + "_true_pT", &m_true_daughter_pt[daughter_id], TString(daughter_number) + "_true_pT/F");
  m_tree->Branch(TString(daughter_number) + "_true_EV_x", &m_true_daughter_vertex_x[daughter_id], TString(daughter_number) + "_true_EV_x/F");
  m_tree->Branch(TString(daughter_number) + "_true_EV_y", &m_true_daughter_vertex_y[daughter_id], TString(daughter_number) + "_true_EV_y/F");
  m_tree->Branch(TString(daughter_number) + "_true_EV_z", &m_true_daughter_vertex_z[daughter_id], TString(daughter_number) + "_true_EV_z/F");
  if (m_constrain_to_vertex_truthMatch)
  {
    m_tree->Branch(TString(daughter_number) + "_true_PV_x", &m_true_daughter_pv_x[daughter_id], TString(daughter_number) + "_true_PV_x/F");
    m_tree->Branch(TString(daughter_number) + "_true_PV_y", &m_true_daughter_pv_y[daughter_id], TString(daughter_number) + "_true_PV_y/F");
    m_tree->Branch(TString(daughter_number) + "_true_PV_z", &m_true_daughter_pv_z[daughter_id], TString(daughter_number) + "_true_PV_z/F");
  }
  m_tree->Branch(TString(daughter_number) + "_true_track_history_PDG_ID", &m_true_daughter_track_history_PDG_ID[daughter_id]);
  m_tree->Branch(TString(daughter_number) + "_true_track_history_PDG_mass", &m_true_daughter_track_history_PDG_mass[daughter_id]);
  m_tree->Branch(TString(daughter_number) + "_true_track_history_px", &m_true_daughter_track_history_px[daughter_id]);
  m_tree->Branch(TString(daughter_number) + "_true_track_history_py", &m_true_daughter_track_history_py[daughter_id]);
  m_tree->Branch(TString(daughter_number) + "_true_track_history_pz", &m_true_daughter_track_history_pz[daughter_id]);
  m_tree->Branch(TString(daughter_number) + "_true_track_history_pE", &m_true_daughter_track_history_pE[daughter_id]);
  m_tree->Branch(TString(daughter_number) + "_true_track_history_pT", &m_true_daughter_track_history_pT[daughter_id]);
}

void KFParticle_truthAndDetTools::fillTruthBranch(PHCompositeNode *topNode, TTree * /*m_tree*/, KFParticle daughter, int daughter_id, KFParticle vertex, bool m_constrain_to_vertex_truthMatch)
{
  float true_px, true_py, true_pz, true_p, true_pt;

  PHNodeIterator nodeIter(topNode);
  PHNode *findNode = dynamic_cast<PHNode *>(nodeIter.findFirst(m_trk_map_node_name_nTuple));
  if (findNode)
  {
    dst_trackmap = findNode::getClass<SvtxTrackMap>(topNode, m_trk_map_node_name_nTuple);
  }
  else
  {
    std::cout << "KFParticle truth matching: " << m_trk_map_node_name_nTuple << " does not exist" << std::endl;
  }
  findNode = dynamic_cast<PHNode *>(nodeIter.findFirst(m_vtx_map_node_name_nTuple));
  if (findNode)
  {
    dst_vertexmap = findNode::getClass<SvtxVertexMap>(topNode, m_vtx_map_node_name_nTuple);
  }
  else
  {
    std::cout << "KFParticle truth matching: " << m_vtx_map_node_name_nTuple << " does not exist" << std::endl;
  }

  track = getTrack(daughter.Id(), dst_trackmap);
  g4particle = getTruthTrack(track, topNode);

  bool isParticleValid = g4particle == nullptr ? 0 : 1;

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
    //clustereval = m_svtx_evalstack->get_cluster_eval();
    //hiteval = m_svtx_evalstack->get_hit_eval();
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
    //Calculate true DCA
    SvtxVertex *recoVertex = getVertex(vertex.Id(), dst_vertexmap);
    PHG4VtxPoint *truePoint = vertexeval->max_truth_point_by_ntracks(recoVertex);

    KFParticle trueKFParticleVertex;

    float f_vertexParameters[6] = {0};

    if (truePoint == NULL)
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

    m_true_daughter_pv_x[daughter_id] = truePoint == NULL ? -99. : truePoint->get_x();
    m_true_daughter_pv_y[daughter_id] = truePoint == NULL ? -99. : truePoint->get_y();
    m_true_daughter_pv_z[daughter_id] = truePoint == NULL ? -99. : truePoint->get_z();
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
  m_true_daughter_track_history_pT[daughter_id].push_back((Float_t) pT);
}

void KFParticle_truthAndDetTools::fillHepMCBranch(HepMC::GenParticle *particle, int daughter_id)
{
  HepMC::FourVector myFourVector = particle->momentum();

  m_true_daughter_track_history_PDG_ID[daughter_id].push_back(particle->pdg_id());
  m_true_daughter_track_history_PDG_mass[daughter_id].push_back((Float_t) particle->generatedMass());
  m_true_daughter_track_history_px[daughter_id].push_back((Float_t) myFourVector.px());
  m_true_daughter_track_history_py[daughter_id].push_back((Float_t) myFourVector.py());
  m_true_daughter_track_history_pz[daughter_id].push_back((Float_t) myFourVector.pz());
  m_true_daughter_track_history_pE[daughter_id].push_back((Float_t) myFourVector.e());
  m_true_daughter_track_history_pT[daughter_id].push_back((Float_t) myFourVector.perp());
}

int KFParticle_truthAndDetTools::getHepMCInfo(PHCompositeNode *topNode, TTree * /*m_tree*/, KFParticle daughter, int daughter_id)
{
  //Make dummy particle for null pointers and missing nodes
  HepMC::GenParticle *dummyParticle = new HepMC::GenParticle();
  HepMC::FourVector dummyFourVector(0, 0, 0, 0);
  dummyParticle->set_momentum(dummyFourVector);
  dummyParticle->set_pdg_id(0);
  dummyParticle->set_generated_mass(0.);

  PHNodeIterator nodeIter(topNode);
  PHNode *findNode = dynamic_cast<PHNode *>(nodeIter.findFirst(m_trk_map_node_name_nTuple));
  if (findNode)
  {
    dst_trackmap = findNode::getClass<SvtxTrackMap>(topNode, m_trk_map_node_name_nTuple);
  }
  else
  {
    std::cout << "KFParticle truth matching: " << m_trk_map_node_name_nTuple << " does not exist" << std::endl;
  }

  track = getTrack(daughter.Id(), dst_trackmap);
  g4particle = getTruthTrack(track, topNode);

  bool isParticleValid = g4particle == nullptr ? 0 : 1;

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

  //Start by looking for our particle in the Geant record
  //Any decay that Geant4 handles will not be in the HepMC record
  //This can happen if you limit the decay volume in the generator
  if (g4particle->get_parent_id() != 0)
  {
    PHNode *findNode = dynamic_cast<PHNode *>(nodeIter.findFirst("G4TruthInfo"));
    if (findNode)
    {
      m_truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
    }
    else
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

  int forbiddenPDGIDs[] = {21, 22};  //Stop tracing history when we reach quarks, gluons and photons

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
        if (breakOut) break;
      }
    }
  }

  return 0;
}  //End of function

void KFParticle_truthAndDetTools::initializeCaloBranches(TTree *m_tree, int daughter_id, std::string daughter_number)
{
  m_tree->Branch(TString(daughter_number) + "_EMCAL_DeltaPhi", &detector_emcal_deltaphi[daughter_id], TString(daughter_number) + "_EMCAL_DeltaPhi/F");
  m_tree->Branch(TString(daughter_number) + "_EMCAL_DeltaEta", &detector_emcal_deltaeta[daughter_id], TString(daughter_number) + "_EMCAL_DeltaEta/F");
  m_tree->Branch(TString(daughter_number) + "_EMCAL_energy_3x3", &detector_emcal_energy_3x3[daughter_id], TString(daughter_number) + "_EMCAL_energy_3x3/F");
  m_tree->Branch(TString(daughter_number) + "_EMCAL_energy_5x5", &detector_emcal_energy_5x5[daughter_id], TString(daughter_number) + "_EMCAL_energy_5x5/F");
  m_tree->Branch(TString(daughter_number) + "_EMCAL_energy_cluster", &detector_emcal_cluster_energy[daughter_id], TString(daughter_number) + "_EMCAL_energy_cluster/F");
  m_tree->Branch(TString(daughter_number) + "_IHCAL_DeltaPhi", &detector_ihcal_deltaphi[daughter_id], TString(daughter_number) + "_IHCAL_DeltaPhi/F");
  m_tree->Branch(TString(daughter_number) + "_IHCAL_DeltaEta", &detector_ihcal_deltaeta[daughter_id], TString(daughter_number) + "_IHCAL_DeltaEta/F");
  m_tree->Branch(TString(daughter_number) + "_IHCAL_energy_3x3", &detector_ihcal_energy_3x3[daughter_id], TString(daughter_number) + "_IHCAL_energy_3x3/F");
  m_tree->Branch(TString(daughter_number) + "_IHCAL_energy_5x5", &detector_ihcal_energy_5x5[daughter_id], TString(daughter_number) + "_IHCAL_energy_5x5/F");
  m_tree->Branch(TString(daughter_number) + "_IHCAL_energy_cluster", &detector_ihcal_cluster_energy[daughter_id], TString(daughter_number) + "_IHCAL_energy_cluster/F");
  m_tree->Branch(TString(daughter_number) + "_OHCAL_DeltaPhi", &detector_ohcal_deltaphi[daughter_id], TString(daughter_number) + "_OHCAL_DeltaEta/F");
  m_tree->Branch(TString(daughter_number) + "_OHCAL_DeltaEta", &detector_ohcal_deltaeta[daughter_id], TString(daughter_number) + "_OHCAL_DeltaEta/F");
  m_tree->Branch(TString(daughter_number) + "_OHCAL_energy_3x3", &detector_ohcal_energy_3x3[daughter_id], TString(daughter_number) + "_OHCAL_energy_3x3/F");
  m_tree->Branch(TString(daughter_number) + "_OHCAL_energy_5x5", &detector_ohcal_energy_5x5[daughter_id], TString(daughter_number) + "_OHCAL_energy_5x5/F");
  m_tree->Branch(TString(daughter_number) + "_OHCAL_energy_cluster", &detector_ohcal_cluster_energy[daughter_id], TString(daughter_number) + "_OHCAL_energy_cluster/F");
}

void KFParticle_truthAndDetTools::fillCaloBranch(PHCompositeNode *topNode,
                                                 TTree * /*m_tree*/, KFParticle daughter, int daughter_id)
{
  PHNodeIterator nodeIter(topNode);
  PHNode *findNode = dynamic_cast<PHNode *>(nodeIter.findFirst(m_trk_map_node_name_nTuple));
  if (findNode)
  {
    dst_trackmap = findNode::getClass<SvtxTrackMap>(topNode, m_trk_map_node_name_nTuple);
  }
  else
  {
    std::cout << "KFParticle truth matching: " << m_trk_map_node_name_nTuple << " does not exist" << std::endl;
  }

  track = getTrack(daughter.Id(), dst_trackmap);

  detector_emcal_deltaphi[daughter_id] = track->get_cal_dphi(SvtxTrack::CAL_LAYER(1));
  detector_emcal_deltaeta[daughter_id] = track->get_cal_deta(SvtxTrack::CAL_LAYER(1));
  detector_emcal_energy_3x3[daughter_id] = track->get_cal_energy_3x3(SvtxTrack::CAL_LAYER(1));
  detector_emcal_energy_5x5[daughter_id] = track->get_cal_energy_5x5(SvtxTrack::CAL_LAYER(1));
  detector_emcal_cluster_energy[daughter_id] = track->get_cal_cluster_e(SvtxTrack::CAL_LAYER(1));

  detector_ihcal_deltaphi[daughter_id] = track->get_cal_dphi(SvtxTrack::CAL_LAYER(2));
  detector_ihcal_deltaeta[daughter_id] = track->get_cal_deta(SvtxTrack::CAL_LAYER(2));
  detector_ihcal_energy_3x3[daughter_id] = track->get_cal_energy_3x3(SvtxTrack::CAL_LAYER(2));
  detector_ihcal_energy_5x5[daughter_id] = track->get_cal_energy_5x5(SvtxTrack::CAL_LAYER(2));
  detector_ihcal_cluster_energy[daughter_id] = track->get_cal_cluster_e(SvtxTrack::CAL_LAYER(2));

  detector_ohcal_deltaphi[daughter_id] = track->get_cal_dphi(SvtxTrack::CAL_LAYER(3));
  detector_ohcal_deltaeta[daughter_id] = track->get_cal_deta(SvtxTrack::CAL_LAYER(3));
  detector_ohcal_energy_3x3[daughter_id] = track->get_cal_energy_3x3(SvtxTrack::CAL_LAYER(3));
  detector_ohcal_energy_5x5[daughter_id] = track->get_cal_energy_5x5(SvtxTrack::CAL_LAYER(3));
  detector_ohcal_cluster_energy[daughter_id] = track->get_cal_cluster_e(SvtxTrack::CAL_LAYER(3));
}

void KFParticle_truthAndDetTools::initializeDetectorBranches(TTree *m_tree, int daughter_id, std::string daughter_number)
{
  m_tree->Branch(TString(daughter_number) + "_local_x", &detector_local_x[daughter_id]);
  m_tree->Branch(TString(daughter_number) + "_local_y", &detector_local_y[daughter_id]);
  m_tree->Branch(TString(daughter_number) + "_local_z", &detector_local_z[daughter_id]);
  m_tree->Branch(TString(daughter_number) + "_layer", &detector_layer[daughter_id]);

  for (auto const &subdetector : Use)
  {
    if (subdetector.second) initializeSubDetectorBranches(m_tree, subdetector.first, daughter_id, daughter_number);
  }
}

void KFParticle_truthAndDetTools::initializeSubDetectorBranches(TTree *m_tree, std::string detectorName, int daughter_id, std::string daughter_number)
{
  if (detectorName == "MVTX")
  {
    m_tree->Branch(TString(daughter_number) + "_" + TString(detectorName) + "_staveID", &mvtx_staveID[daughter_id]);
    m_tree->Branch(TString(daughter_number) + "_" + TString(detectorName) + "_chipID", &mvtx_chipID[daughter_id]);
  }
  if (detectorName == "INTT")
  {
    m_tree->Branch(TString(daughter_number) + "_" + TString(detectorName) + "_ladderZID", &intt_ladderZID[daughter_id]);
    m_tree->Branch(TString(daughter_number) + "_" + TString(detectorName) + "_ladderPhiID", &intt_ladderPhiID[daughter_id]);
  }
  if (detectorName == "TPC")
  {
    m_tree->Branch(TString(daughter_number) + "_" + TString(detectorName) + "_sectorID", &tpc_sectorID[daughter_id]);
    m_tree->Branch(TString(daughter_number) + "_" + TString(detectorName) + "_side", &tpc_side[daughter_id]);
  }
}

void KFParticle_truthAndDetTools::fillDetectorBranch(PHCompositeNode *topNode,
                                                     TTree * /*m_tree*/, KFParticle daughter, int daughter_id)
{
  PHNodeIterator nodeIter(topNode);

  PHNode *findNode = dynamic_cast<PHNode *>(nodeIter.findFirst(m_trk_map_node_name_nTuple));
  if (findNode)
  {
    dst_trackmap = findNode::getClass<SvtxTrackMap>(topNode, m_trk_map_node_name_nTuple);
  }
  else
  {
    std::cout << "KFParticle truth matching: " << m_trk_map_node_name_nTuple << " does not exist" << std::endl;
  }

  findNode = dynamic_cast<PHNode *>(nodeIter.findFirst("TRKR_CLUSTER"));
  if (findNode)
  {
    dst_clustermap = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  }
  else
  {
    std::cout << "KFParticle detector info: TRKR_CLUSTER does not exist" << std::endl;
  }

  track = getTrack(daughter.Id(), dst_trackmap);

  for (SvtxTrack::ConstClusterKeyIter iter = track->begin_cluster_keys();
       iter != track->end_cluster_keys();
       ++iter)
  {
    TrkrDefs::cluskey clusKey = *iter;
    TrkrCluster *cluster = dst_clustermap->findCluster(clusKey);
    const unsigned int trkrId = TrkrDefs::getTrkrId(clusKey);

    detector_local_x[daughter_id].push_back(cluster->getX());
    detector_local_y[daughter_id].push_back(cluster->getY());
    detector_local_z[daughter_id].push_back(cluster->getZ());
    detector_layer[daughter_id].push_back(TrkrDefs::getLayer(clusKey));
    unsigned int staveId, chipId, ladderZId, ladderPhiId, sectorId, side;
    staveId = chipId = ladderZId = ladderPhiId = sectorId = side = -99;

    if (Use["MVTX"] && trkrId == TrkrDefs::mvtxId)
    {
      staveId = MvtxDefs::getStaveId(clusKey);
      chipId = MvtxDefs::getChipId(clusKey);
    }
    else if (Use["INTT"] && trkrId == TrkrDefs::inttId)
    {
      ladderZId = InttDefs::getLadderZId(clusKey);
      ladderPhiId = InttDefs::getLadderPhiId(clusKey);
    }
    else if (Use["TPC"] && trkrId == TrkrDefs::tpcId)
    {
      sectorId = TpcDefs::getSectorId(clusKey);
      side = TpcDefs::getSide(clusKey);
    }

    mvtx_staveID[daughter_id].push_back(staveId);
    mvtx_chipID[daughter_id].push_back(chipId);
    intt_ladderZID[daughter_id].push_back(ladderZId);
    intt_ladderPhiID[daughter_id].push_back(ladderPhiId);
    tpc_sectorID[daughter_id].push_back(sectorId);
    tpc_side[daughter_id].push_back(side);
  }
}

void KFParticle_truthAndDetTools::allPVInfo(PHCompositeNode *topNode,
                                            TTree * /*m_tree*/,
                                            KFParticle motherParticle,
                                            std::vector<KFParticle> daughters,
                                            std::vector<KFParticle> intermediates)
{
  KFParticle_Tools kfpTupleTools;
  std::vector<KFParticle> primaryVertices = kfpTupleTools.makeAllPrimaryVertices(topNode, m_vtx_map_node_name_nTuple);

  for (unsigned int i = 0; i < primaryVertices.size(); ++i)
  {
    allPV_x.push_back(primaryVertices[i].GetX());
    allPV_y.push_back(primaryVertices[i].GetY());
    allPV_z.push_back(primaryVertices[i].GetZ());

    allPV_mother_IP.push_back(motherParticle.GetDistanceFromVertex(primaryVertices[i]));
    allPV_mother_IPchi2.push_back(motherParticle.GetDeviationFromVertex(primaryVertices[i]));

    for (unsigned int j = 0; j < daughters.size(); ++j)
    {
      allPV_daughter_IP[j].push_back(daughters[j].GetDistanceFromVertex(primaryVertices[i]));
      allPV_daughter_IPchi2[j].push_back(daughters[j].GetDeviationFromVertex(primaryVertices[i]));
    }

    for (unsigned int j = 0; j < intermediates.size(); ++j)
    {
      allPV_intermediates_IP[j].push_back(intermediates[j].GetDistanceFromVertex(primaryVertices[i]));
      allPV_intermediates_IPchi2[j].push_back(intermediates[j].GetDeviationFromVertex(primaryVertices[i]));
    }
  }
}

void KFParticle_truthAndDetTools::clearVectors()
{
  for (int i = 0; i < m_num_tracks_nTuple; ++i)
  {
    //Truth vectors
    m_true_daughter_track_history_PDG_ID[i].clear();
    m_true_daughter_track_history_PDG_mass[i].clear();
    m_true_daughter_track_history_px[i].clear();
    m_true_daughter_track_history_py[i].clear();
    m_true_daughter_track_history_pz[i].clear();
    m_true_daughter_track_history_pE[i].clear();
    m_true_daughter_track_history_pT[i].clear();

    //Detector vectors
    detector_local_x[i].clear();
    detector_local_y[i].clear();
    detector_local_z[i].clear();
    detector_layer[i].clear();
    mvtx_staveID[i].clear();
    mvtx_chipID[i].clear();
    intt_ladderZID[i].clear();
    intt_ladderPhiID[i].clear();
    tpc_sectorID[i].clear();
    tpc_side[i].clear();

    //PV vectors
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
