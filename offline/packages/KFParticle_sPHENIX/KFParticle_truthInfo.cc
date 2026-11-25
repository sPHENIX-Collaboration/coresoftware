#include "KFParticle_truthInfo.h"
#include "KFParticle_util.h"

#include <g4eval/SvtxEvalStack.h>   // for SvtxEvalStack
#include <g4eval/SvtxTrackEval.h>   // for SvtxTrackEval
#include <g4eval/SvtxTruthEval.h>   // for SvtxTruthEval
#include <g4eval/SvtxVertexEval.h>  // for SvtxVertexEval

#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMap.h>
#include <globalvertex/MbdVertex.h>
#include <globalvertex/MbdVertexMap.h>
#include <globalvertex/SvtxVertex.h>     // for SvtxVertex
#include <globalvertex/SvtxVertexMap.h>  // for SvtxVertexMap, SvtxVer...

#include <trackbase_historic/SvtxPHG4ParticleMap.h>
#include <trackbase_historic/SvtxTrack.h>     // for SvtxTrack
#include <trackbase_historic/SvtxTrackMap.h>  // for SvtxTrackMap

#include <g4main/PHG4Particle.h>            // for PHG4Particle
#include <g4main/PHG4TruthInfoContainer.h>  // for PHG4TruthInfoContainer
#include <g4main/PHG4VtxPoint.h>            // for PHG4VtxPoint

#include <phhepmc/PHHepMCGenEvent.h>     // for PHHepMCGenEvent
#include <phhepmc/PHHepMCGenEventMap.h>  // for PHHepMCGenEventMap

#include <phool/getClass.h>  // for getClass

#include <KFParticle.h>  // for KFParticle

#include <TSystem.h>  // for gSystem->Exit()
#include <TTree.h>    // for TTree

#include <HepMC/GenEvent.h>       // for GenEvent::particle_con...
#include <HepMC/GenParticle.h>    // for GenParticle
#include <HepMC/GenVertex.h>      // for GenVertex::particle_it...
#include <HepMC/IteratorRange.h>  // for parents
#include <HepMC/SimpleVector.h>   // for FourVector

#include <algorithm>  // for find
#include <cmath>      // for pow, sqrt
#include <cstdlib>    // for abs, exit
#include <iostream>   // for operator<<, endl
#include <iterator>   // for begin, end
#include <set>
#include <string>

GlobalVertex *KFParticle_truthInfo::getVertex(unsigned int vertex_id, GlobalVertexMap *vertexmap)
{
  return vertexmap->get(vertex_id);
}

PHG4Particle *KFParticle_truthInfo::getTruthTrack(SvtxTrack *thisTrack, PHCompositeNode *topNode)
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

    const auto &truth_set = dst_reco_truth_map->get(thisTrack->get_id());
    if (!truth_set.empty())
    {
      const auto &best_weight = *truth_set.rbegin();
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

void KFParticle_truthInfo::initializeBranches(TTree *m_tree, int daughter_id, const std::string &daughter_number, bool m_constrain_to_vertex_truthMatch)
{
  m_tree->Branch((daughter_number + "_true_ID").c_str(), &m_true_daughter_id[daughter_id], (daughter_number + "_true_ID/I").c_str());
  if (m_constrain_to_vertex_truthMatch)
  {
    m_tree->Branch((daughter_number + "_true_PV_DCA").c_str(), &m_true_daughter_ip[daughter_id], (daughter_number + "_true_PV_DCA/F").c_str());
    m_tree->Branch((daughter_number + "_true_PV_DCA_xy").c_str(), &m_true_daughter_ip_xy[daughter_id], (daughter_number + "_true_PV_DCA_xy/F").c_str());
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

void KFParticle_truthInfo::fill(PHCompositeNode *topNode, const KFParticle &daughter, int daughter_id, const KFParticle &kfvertex, bool m_constrain_to_vertex_truthMatch)
{
  float true_px;
  float true_py;
  float true_pz;
  float true_p;
  float true_pt;

  auto *dst_trackmap = findNode::getClass<SvtxTrackMap>(topNode, m_config.trackMapNodeName);
  if (!dst_trackmap)
  {
    std::cout << "KFParticle truth matching: " << m_config.trackMapNodeName << " does not exist" << std::endl;
  }

  if (m_config.useMbdVertex)
  {
    dst_mbdvertexmap = findNode::getClass<MbdVertexMap>(topNode, "MbdVertexMap");
    if (!dst_mbdvertexmap)
    {
      std::cout << "KFParticle truth matching: MbdVertexMap does not exist" << std::endl;
    }
  }
  else
  {
    dst_vertexmap = findNode::getClass<SvtxVertexMap>(topNode, m_config.vertexMapNodeName);
    if (!dst_vertexmap)
    {
      std::cout << "KFParticle truth matching: " << m_config.vertexMapNodeName << " does not exist" << std::endl;
    }
  }

  auto *globalvertexmap = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");
  if (!globalvertexmap)
  {
    std::cout << "KFParticle truth matching: GlobalVertexMap does not exist" << std::endl;
  }

  SvtxTrack *track = KFParticle_util::getTrack(daughter.Id(), dst_trackmap);
  g4particle = getTruthTrack(track, topNode);

  const bool isParticleValid = g4particle != nullptr;

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

  m_svtx_evalstack->next_event(topNode);
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
    GlobalVertex::VTXTYPE whichVtx = m_config.useMbdVertex ? GlobalVertex::MBD : GlobalVertex::SVTX;
    auto svtxviter = recoVertex->find_vertexes(whichVtx);

    // check that it contains a track vertex
    if (svtxviter == recoVertex->end_vertexes())
    {
      std::string vtxType = m_config.useMbdVertex ? "MBD" : "silicon";
      std::cout << "Have a global vertex with no " << vtxType << " vertex... shouldn't happen in KFParticle_truthInfo::fill..." << std::endl;
    }

    const auto &svtxvertexvector = svtxviter->second;
    const Vertex *dstVertex = nullptr;

    for (const auto &vertex_iter : svtxvertexvector)
    {
      if (m_config.useMbdVertex)
      {
        dstVertex = dst_mbdvertexmap->find(vertex_iter->get_id())->second;
      }
      else
      {
        dstVertex = dst_vertexmap->find(vertex_iter->get_id())->second;
      }
    }

    PHG4VtxPoint *truePoint = nullptr;
    if (m_config.useMbdVertex)
    {
      std::set<PHG4VtxPoint *> truePointSet = vertexeval->all_truth_points(dstVertex);
      if (!truePointSet.empty())
      {
        truePoint = *truePointSet.begin();
      }
    }
    else
    {
      truePoint = vertexeval->max_truth_point_by_ntracks(dstVertex);
    }

    if (truePoint == nullptr && isParticleValid)
    {
      PHG4Particle *g4mother = trutheval->get_parent_particle(g4particle);
      if (!g4mother)
      {
        std::cout << "KFParticle truth matching: True mother not found!\n";
        std::cout << "Your truth track DCA will be measured wrt a reconstructed vertex!" << std::endl;
        truePoint = nullptr;
      }
      else
      {
        truePoint = trutheval->get_vertex(g4mother);
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

void KFParticle_truthInfo::fillGeant4Branch(PHG4Particle *particle, int daughter_id)
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

void KFParticle_truthInfo::fillHepMCBranch(HepMC::GenParticle *particle, int daughter_id)
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

int KFParticle_truthInfo::getHepMCInfo(PHCompositeNode *topNode, const KFParticle &daughter, int daughter_id)
{
  // Make dummy particle for null pointers and missing nodes
  HepMC::GenParticle dummyParticle;
  HepMC::FourVector dummyFourVector(0, 0, 0, 0);
  dummyParticle.set_momentum(dummyFourVector);
  dummyParticle.set_pdg_id(0);
  dummyParticle.set_generated_mass(0.);

  auto *dst_trackmap = findNode::getClass<SvtxTrackMap>(topNode, m_config.trackMapNodeName);
  if (!dst_trackmap)
  {
    std::cout << "KFParticle truth matching: " << m_config.trackMapNodeName << " does not exist" << std::endl;
    gSystem->Exit(1);
    exit(1);
  }

  SvtxTrack *track = KFParticle_util::getTrack(daughter.Id(), dst_trackmap);
  g4particle = getTruthTrack(track, topNode);

  const bool isParticleValid = g4particle != nullptr;

  if (!isParticleValid)
  {
    std::cout << "KFParticle truth matching: this track is a ghost" << std::endl;
    fillHepMCBranch(&dummyParticle, daughter_id);
    return 0;
  }

  m_geneventmap = findNode::getClass<PHHepMCGenEventMap>(topNode, "PHHepMCGenEventMap");
  if (!m_geneventmap)
  {
    std::cout << "KFParticle truth matching: Missing node PHHepMCGenEventMap" << std::endl;
    std::cout << "You will have no mother information" << std::endl;
    fillHepMCBranch(&dummyParticle, daughter_id);
    return 0;
  }

  m_genevt = m_geneventmap->get(1);
  if (!m_genevt)
  {
    std::cout << "KFParticle truth matching: Missing node PHHepMCGenEvent" << std::endl;
    std::cout << "You will have no mother information" << std::endl;
    fillHepMCBranch(&dummyParticle, daughter_id);
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
                        std::abs((*mother)->pdg_id())) != std::end(forbiddenPDGIDs))
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

void KFParticle_truthInfo::clear(int nTracks)
{
  for (int i = 0; i < nTracks; ++i)
  {
    // Truth vectors
    m_true_daughter_track_history_PDG_ID[i].clear();
    m_true_daughter_track_history_PDG_mass[i].clear();
    m_true_daughter_track_history_px[i].clear();
    m_true_daughter_track_history_py[i].clear();
    m_true_daughter_track_history_pz[i].clear();
    m_true_daughter_track_history_pE[i].clear();
    m_true_daughter_track_history_pT[i].clear();
  }
}
