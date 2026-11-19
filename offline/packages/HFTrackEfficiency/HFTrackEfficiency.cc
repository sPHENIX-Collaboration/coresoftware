#include "HFTrackEfficiency.h"

#include <decayfinder/DecayFinderContainerBase.h>

#include <trackbase/TrkrDefs.h>

#include <trackbase_historic/PHG4ParticleSvtxMap.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrackMap_v2.h>
#include <trackbase_historic/TrackSeed.h>
#include <trackbase_historic/TrackSeedContainer.h>

#include <g4eval/SvtxClusterEval.h>
#include <g4eval/SvtxTruthEval.h>

#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4VtxPoint.h>

#include <phhepmc/PHHepMCGenEvent.h>
#include <phhepmc/PHHepMCGenEventMap.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <TDatabasePDG.h>
#include <TFile.h>
#include <TTree.h>

#include <CLHEP/Vector/LorentzRotation.h>
#include <CLHEP/Vector/LorentzVector.h>
#include <CLHEP/Vector/ThreeVector.h>

#include <HepMC/GenEvent.h>
#include <HepMC/GenParticle.h>
#include <HepMC/GenVertex.h>  // for GenVertex::particle_iterator
#include <HepMC/SimpleVector.h>
#include <HepMC/Units.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <set>

namespace
{
  // Speed of light in cm/ps, to turn a flight distance into a proper decay time
  constexpr double c_cm_per_ps = 0.0299792458;

  // Follow the first incoming particle of every vertex back to the start of the HepMC record
  HepMC::GenVertex *findRootVertex(const HepMC::GenParticle *particle)
  {
    HepMC::GenVertex *vtx = particle->production_vertex();
    constexpr int maxDepth = 200;  // Guard against a malformed record
    for (int depth = 0; vtx && depth < maxDepth; ++depth)
    {
      if (vtx->particles_in_size() == 0)
      {
        break;
      }
      const HepMC::GenParticle *parent = *(vtx->particles_in_const_begin());
      if (!parent || !parent->production_vertex())
      {
        break;
      }
      vtx = parent->production_vertex();
    }
    return vtx;
  }

  // Position of a HepMC vertex in cm in the lab frame, the same transformation that HepMCNodeReader applies
  // (generator frame -> lab frame, units to cm, plus the collision vertex) but without any extra reader offset
  CLHEP::Hep3Vector getLabFrameVertex(const HepMC::GenVertex &vtx, const HepMC::GenEvent &evt, const CLHEP::HepLorentzRotation &evtGen2Lab, const HepMC::FourVector &collisionVertex)
  {
    const double lengthFactor = HepMC::Units::conversion_factor(evt.length_unit(), HepMC::Units::CM);
    const CLHEP::HepLorentzVector genVertex(vtx.position().x(), vtx.position().y(), vtx.position().z(), vtx.position().t());
    const CLHEP::HepLorentzVector labVertex = evtGen2Lab(genVertex);
    return CLHEP::Hep3Vector(labVertex.x() * lengthFactor + collisionVertex.x(),
                             labVertex.y() * lengthFactor + collisionVertex.y(),
                             labVertex.z() * lengthFactor + collisionVertex.z());
  }
}  // namespace

//____________________________________________________________________________..
HFTrackEfficiency::HFTrackEfficiency(const std::string &name)
  : SubsysReco(name)
  , m_input_track_map_node_name("SvtxTrackMap")
  , m_output_track_map_node_name("HFSelected")
  , m_outfile_name("outputHFTrackEff.root")
  , m_truthRecoMatchPercent(5.)
  , m_nDaughters(2)
{
}

//____________________________________________________________________________..
int HFTrackEfficiency::Init(PHCompositeNode *topNode)
{
  m_truthRecoMatchPercent /= 100.;

  if (m_write_track_map)
  {
    PHNodeIterator nodeIter(topNode);

    PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(nodeIter.findFirst("PHCompositeNode", "DST"));
    if (!dstNode)
    {
      dstNode = new PHCompositeNode("DST");
      topNode->addNode(dstNode);
      std::cout << "DST node added" << std::endl;
    }

    outputNodeName = m_output_track_map_node_name + "_SvtxTrackMap";
    m_output_trackMap = new SvtxTrackMap_v2();
    PHIODataNode<PHObject> *outputTrackNode = new PHIODataNode<PHObject>(m_output_trackMap, outputNodeName, "PHObject");
    dstNode->addNode(outputTrackNode);
    std::cout << outputNodeName << " node added" << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int HFTrackEfficiency::process_event(PHCompositeNode *topNode)
{
  PHNodeIterator nodeIter(topNode);
  PHCompositeNode *dstNode = static_cast<PHCompositeNode *>(nodeIter.findFirst("PHCompositeNode", "DST"));
  assert(dstNode);

  std::string df_node_name = m_df_module_name + "_DecayMap";
  m_decayMap = findNode::getClass<DecayFinderContainerBase>(topNode, df_node_name);
  if (!m_decayMap)
  {
    std::cout << __FILE__ << ": Missing node " << df_node_name << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_input_trackMap = findNode::getClass<SvtxTrackMap>(topNode, m_input_track_map_node_name);
  if (!m_input_trackMap)
  {
    std::cout << __FILE__ << ": Missing node " << m_input_track_map_node_name << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_geneventmap = findNode::getClass<PHHepMCGenEventMap>(topNode, "PHHepMCGenEventMap");
  if (!m_geneventmap)
  {
    std::cout << __FILE__ << ": Missing node PHHepMCGenEventMap" << std::endl;
  }

  m_truthInfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  if (!m_truthInfo)
  {
    std::cout << __FILE__ << ": Missing node G4TruthInfo" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // Not required, only used to find the clusters of tracks that store seed IDs rather than seed pointers
  m_siliconSeeds = findNode::getClass<TrackSeedContainer>(topNode, "SiliconTrackSeedContainer");
  m_tpcSeeds = findNode::getClass<TrackSeedContainer>(topNode, "TpcTrackSeedContainer");

  if (!m_svtx_evalstack)
  {
    m_svtx_evalstack = new SvtxEvalStack(topNode);
    trackeval = m_svtx_evalstack->get_track_eval();
    trackeval->set_track_nodename(m_input_track_map_node_name);
  }
  m_svtx_evalstack->next_event(topNode);

  if (m_decay_descriptor.empty() && !m_decayMap->empty())
  {
    getDecayDescriptor();
    getNDaughters();
    if (m_nDaughters > m_maxTracks)
    {
      std::cout << __FILE__ << ": The decay has " << m_nDaughters << " trackable daughters but at most " << m_maxTracks << " are supported" << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
    if (m_nDaughters != 2 && m_nDaughters != 3)
    {
      std::cout << __FILE__ << ": Warning, " << m_nDaughters << " trackable daughters. Decays found in the Geant4 record only support 2 or 3" << std::endl;
    }
  }

  bool reconstructableDecay = false;
  for (auto &iter : *m_decayMap)
  {
    Decay decay = iter.second;
    ++m_counter_allDecays;

    m_all_tracks_reconstructed = findTracks(topNode, decay);
    if (m_all_tracks_reconstructed)
    {
      ++m_counter_acceptedDecays;
      reconstructableDecay = true;
    }

    if (m_write_nTuple)
    {
      if (m_counter_allDecays == 1)
      {
        initializeBranches();
      }
      m_tree->Fill();
    }
  }

  if (m_triggerOnDecay && !reconstructableDecay)
  {
    if (Verbosity() >= VERBOSITY_MORE)
    {
      std::cout << "No decays were reconstructable in this event, skipping" << std::endl;
    }
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int HFTrackEfficiency::End(PHCompositeNode * /*topNode*/)
{
  PrintEff();

  if (m_write_nTuple && m_counter_allDecays != 0)
  {
    m_outfile->Write();
    m_outfile->Close();
    delete m_outfile;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

PHG4Particle *HFTrackEfficiency::findGeant4Primary(const int embeddingID, const int barcode, const int pdgID, const CLHEP::Hep3Vector &labMomentum) const
{
  if (!m_truthInfo)
  {
    return nullptr;
  }

  // HepMCNodeReader copies the HepMC barcode to the Geant4 primary and registers it with the embedding ID of the HepMC event
  // so that is the unambiguous match. The barcode is only unique within one HepMC event, hence the embedding ID.
  // If that fails fall back on the closest momentum (in the lab frame, as the Geant4 particles are)
  PHG4Particle *closest = nullptr;
  double closestDelta = std::numeric_limits<double>::max();
  const double tolerance = std::max(5e-3, 5e-3 * labMomentum.mag());

  PHG4TruthInfoContainer::ConstRange range = m_truthInfo->GetPrimaryParticleRange();
  for (PHG4TruthInfoContainer::ConstIterator iter = range.first; iter != range.second; ++iter)
  {
    PHG4Particle *candidate = iter->second;
    if (candidate->get_pid() != pdgID)
    {
      continue;
    }

    if (candidate->get_barcode() == barcode && m_truthInfo->isEmbeded(candidate->get_track_id()) == embeddingID)
    {
      return candidate;
    }

    const double delta = std::max({std::abs(candidate->get_px() - labMomentum.x()),
                                   std::abs(candidate->get_py() - labMomentum.y()),
                                   std::abs(candidate->get_pz() - labMomentum.z())});
    if (delta <= tolerance && delta < closestDelta)
    {
      closest = candidate;
      closestDelta = delta;
    }
  }

  return closest;
}

bool HFTrackEfficiency::getGeant4Vertices(const PHG4Particle *motherG4, const PHG4Particle *daughterG4, CLHEP::Hep3Vector &primary, CLHEP::Hep3Vector &production, CLHEP::Hep3Vector &decay) const
{
  // Follow the parents of the mother back to the particle that started the chain
  const PHG4Particle *root = motherG4;
  for (int depth = 0; depth < 200 && root->get_parent_id() != 0; ++depth)
  {
    const PHG4Particle *parent = m_truthInfo->GetParticle(root->get_parent_id());
    if (!parent)
    {
      break;
    }
    root = parent;
  }

  const PHG4VtxPoint *primaryVtx = m_truthInfo->GetVtx(root->get_vtx_id());
  const PHG4VtxPoint *productionVtx = m_truthInfo->GetVtx(motherG4->get_vtx_id());
  const PHG4VtxPoint *decayVtx = m_truthInfo->GetVtx(daughterG4->get_vtx_id());
  if (!primaryVtx || !productionVtx || !decayVtx)
  {
    return false;
  }

  primary.set(primaryVtx->get_x(), primaryVtx->get_y(), primaryVtx->get_z());
  production.set(productionVtx->get_x(), productionVtx->get_y(), productionVtx->get_z());
  decay.set(decayVtx->get_x(), decayVtx->get_y(), decayVtx->get_z());

  return true;
}

void HFTrackEfficiency::countClusters(SvtxTrack *track, PHG4Particle *truthParticle, const unsigned int index)
{
  // Same cluster keys that SvtxTrackEval uses, but also valid for tracks that only store the seed IDs
  std::set<TrkrDefs::cluskey> clusterKeys;
  auto collectKeys = [&clusterKeys](const TrackSeed *seed)
  {
    if (seed)
    {
      for (auto iter = seed->begin_cluster_keys(); iter != seed->end_cluster_keys(); ++iter)
      {
        clusterKeys.insert(*iter);
      }
    }
  };

  const TrackSeed *siliconSeed = track->get_silicon_seed();
  if (!siliconSeed && m_siliconSeeds && track->get_silicon_seed_id() < m_siliconSeeds->size())
  {
    siliconSeed = m_siliconSeeds->get(track->get_silicon_seed_id());
  }
  collectKeys(siliconSeed);

  const TrackSeed *tpcSeed = track->get_tpc_seed();
  if (!tpcSeed && m_tpcSeeds && track->get_tpc_seed_id() < m_tpcSeeds->size())
  {
    tpcSeed = m_tpcSeeds->get(track->get_tpc_seed_id());
  }
  collectKeys(tpcSeed);

  SvtxClusterEval *clusterEval = trackeval->get_cluster_eval();
  SvtxTruthEval *truthEval = trackeval->get_truth_eval();

  for (int iDet = 0; iDet < m_nDetectors; ++iDet)
  {
    m_reco_track_nclusters[index][iDet] = 0;
    m_reco_track_ntrue_clusters[index][iDet] = truthParticle ? 0 : -1;
  }

  for (const auto &key : clusterKeys)
  {
    int iDet = -1;
    switch (TrkrDefs::getTrkrId(key))
    {
    case TrkrDefs::mvtxId:
      iDet = 0;
      break;
    case TrkrDefs::inttId:
      iDet = 1;
      break;
    case TrkrDefs::tpcId:
      iDet = 2;
      break;
    case TrkrDefs::micromegasId:
      iDet = 3;
      break;
    default:
      break;
    }
    if (iDet < 0)
    {
      continue;
    }

    ++m_reco_track_nclusters[index][iDet];

    if (truthParticle)
    {
      // A cluster is true if the truth daughter left energy in it, counted once however many particles share the cluster
      for (PHG4Particle *contributor : clusterEval->all_truth_particles(key))
      {
        if (truthEval->are_same_particle(contributor, truthParticle))
        {
          ++m_reco_track_ntrue_clusters[index][iDet];
          break;
        }
      }
    }
  }
}

bool HFTrackEfficiency::findTracks(PHCompositeNode *topNode, const Decay &decay)
{
  resetBranches();

  int trackableParticles[] = {11, 13, 211, 321, 2212};
  bool recoTrackFound = false;
  std::set<SvtxTrack *> selectedTracks;

  CLHEP::HepLorentzVector motherRecoLV;
  CLHEP::HepLorentzVector daughterSumTrueLV;
  CLHEP::HepLorentzVector motherTrueLV;
  CLHEP::HepLorentzVector daughterTrueLV;

  // True vertices, all in cm and in the lab frame that Geant4 (and so the reconstruction) uses
  bool haveVertices = false;
  CLHEP::Hep3Vector primaryVtx;
  CLHEP::Hep3Vector productionVtx;
  CLHEP::Hep3Vector decayVtx;

  // The HepMC record is in the event generator frame. HepMCNodeReader rotates it to the lab frame (crossing angle) and shifts it
  // by the collision vertex, possibly plus a vertex offset/smearing of its own that is not visible from here.
  // The latter is recovered by comparing one HepMC daughter vertex to the Geant4 vertex of the same particle.
  CLHEP::HepLorentzRotation evtGen2Lab;
  CLHEP::Hep3Vector vertexOffset;
  bool vertexOffsetKnown = false;

  HepMC::GenEvent *theEvent = nullptr;
  if (m_geneventmap)
  {
    m_genevt = m_geneventmap->get(decay[0].first.first);
    assert(m_genevt);

    theEvent = m_genevt->getEvent();
    HepMC::GenParticle *mother = theEvent->barcode_to_particle(decay[0].first.second);
    assert(mother);
    if (Verbosity() >= VERBOSITY_MORE)
    {
      mother->print();
    }

    // These are the generator frame values, the crossing angle rotation is not applied. Compare with the Geant4 truth at the 1e-3 level
    m_true_mother_pT = mother->momentum().perp();
    m_true_mother_p = std::sqrt(std::pow(mother->momentum().px(), 2) + std::pow(mother->momentum().py(), 2) + std::pow(mother->momentum().pz(), 2));  // Must have an old HepMC build, no mag function
    m_true_mother_eta = mother->momentum().eta();
    m_true_mother_phi = mother->momentum().phi();
    if (mother->momentum().e() > std::fabs(mother->momentum().pz()))
    {
      m_true_mother_rapidity = 0.5 * log((mother->momentum().e() + mother->momentum().pz()) / (mother->momentum().e() - mother->momentum().pz()));
    }
    else
    {
      m_true_mother_rapidity = -999.;
    }

    evtGen2Lab = m_genevt->get_LorentzRotation_EvtGen2Lab();
    HepMC::GenVertex *motherProductionVtx = mother->production_vertex();
    HepMC::GenVertex *motherDecayVtx = mother->end_vertex();
    HepMC::GenVertex *rootVtx = findRootVertex(mother);
    if (motherProductionVtx && motherDecayVtx && rootVtx)
    {
      const HepMC::FourVector &collisionVtx = m_genevt->get_collision_vertex();
      primaryVtx = getLabFrameVertex(*rootVtx, *theEvent, evtGen2Lab, collisionVtx);
      productionVtx = getLabFrameVertex(*motherProductionVtx, *theEvent, evtGen2Lab, collisionVtx);
      decayVtx = getLabFrameVertex(*motherDecayVtx, *theEvent, evtGen2Lab, collisionVtx);
      haveVertices = true;

      // A prompt mother is produced where the record starts. Compared in the HepMC length unit, no reliance on the record starting at the origin
      constexpr double epsilon = 1e-6;
      m_is_primary = std::abs(motherProductionVtx->position().x() - rootVtx->position().x()) < epsilon &&
                     std::abs(motherProductionVtx->position().y() - rootVtx->position().y()) < epsilon &&
                     std::abs(motherProductionVtx->position().z() - rootVtx->position().z()) < epsilon;
    }
    else if (Verbosity() > 0)
    {
      std::cout << __FILE__ << ": The mother in the HepMC record has no production or decay vertex, no true vertices will be filled" << std::endl;
    }
  }

  // Fill the mother and daughter truth from a Geant4 mother/daughter pair, used when the decay was found in the Geant4 record
  auto fillFromGeant4 = [&](PHG4Particle *motherG4, PHG4Particle *daughterG4, const int daughterPDG, const int daughterIndex)
  {
    if (Verbosity() >= VERBOSITY_MORE)
    {
      daughterG4->identify();
    }

    m_is_primary = m_truthInfo->is_sPHENIX_primary(motherG4);

    CLHEP::Hep3Vector mother3Vector(motherG4->get_px(), motherG4->get_py(), motherG4->get_pz());
    motherTrueLV.setVectM(mother3Vector, getParticleMass(decay[0].second));
    m_true_mother_pT = motherTrueLV.perp();
    m_true_mother_p = mother3Vector.mag();
    m_true_mother_eta = motherTrueLV.pseudoRapidity();
    m_true_mother_phi = motherTrueLV.phi();
    m_true_mother_rapidity = motherTrueLV.rapidity();

    daughterTrueLV.setVectM(CLHEP::Hep3Vector(daughterG4->get_px(), daughterG4->get_py(), daughterG4->get_pz()), getParticleMass(daughterPDG));
    daughterSumTrueLV += daughterTrueLV;

    // A HepMC mother already provided the vertices
    if (!haveVertices)
    {
      haveVertices = getGeant4Vertices(motherG4, daughterG4, primaryVtx, productionVtx, decayVtx);
    }

    m_true_track_PID[daughterIndex] = daughterG4->get_pid();
  };

  int index = -1;

  for (unsigned int i = 1; i < decay.size(); ++i)
  {
    m_dst_track = nullptr;

    if (std::find(std::begin(trackableParticles), std::end(trackableParticles),
                  std::abs(decay[i].second)) != std::end(trackableParticles))
    {
      ++index;
      if (index >= m_maxTracks)
      {
        std::cout << __FILE__ << ": Decay has more than " << m_maxTracks << " trackable daughters, increase m_maxTracks" << std::endl;
        return false;
      }

      PHG4Particle *daughterG4 = nullptr;  // Must not carry over from the previous daughter
      bool haveTruthKinematics = false;

      HepMC::GenParticle *daughterHepMC = nullptr;
      if (theEvent && decay[i].first.second > -1)
      {
        daughterHepMC = theEvent->barcode_to_particle(decay[i].first.second);
      }

      if (daughterHepMC)
      {
        if (Verbosity() >= VERBOSITY_MORE)
        {
          daughterHepMC->print();
        }

        daughterTrueLV.setVectM(CLHEP::Hep3Vector(daughterHepMC->momentum().px(), daughterHepMC->momentum().py(), daughterHepMC->momentum().pz()), getParticleMass(decay[i].second));
        daughterSumTrueLV += daughterTrueLV;
        haveTruthKinematics = true;

        m_true_track_PID[index] = daughterHepMC->pdg_id();

        // We need the G4 particle, not the HepMC one, to use the truth/reco matching.
        // The Geant4 momentum is in the lab frame so rotate the HepMC one before comparing
        const CLHEP::HepLorentzVector daughterLabLV = evtGen2Lab(CLHEP::HepLorentzVector(daughterHepMC->momentum().px(), daughterHepMC->momentum().py(), daughterHepMC->momentum().pz(), daughterHepMC->momentum().e()));
        daughterG4 = findGeant4Primary(decay[i].first.first, decay[i].first.second, decay[i].second, daughterLabLV.vect());

        if (daughterG4 && haveVertices && !vertexOffsetKnown)
        {
          const HepMC::GenVertex *daughterProductionVtx = daughterHepMC->production_vertex();
          const PHG4VtxPoint *daughterG4Vtx = m_truthInfo->GetVtx(daughterG4->get_vtx_id());
          if (daughterProductionVtx && daughterG4Vtx)
          {
            const CLHEP::Hep3Vector expectedVtx = getLabFrameVertex(*daughterProductionVtx, *theEvent, evtGen2Lab, m_genevt->get_collision_vertex());
            vertexOffset.set(daughterG4Vtx->get_x() - expectedVtx.x(), daughterG4Vtx->get_y() - expectedVtx.y(), daughterG4Vtx->get_z() - expectedVtx.z());
            vertexOffsetKnown = true;
            if (Verbosity() > 0 && vertexOffset.mag() > 1e-3)
            {
              std::cout << __FILE__ << ": The Geant4 vertices are shifted by " << vertexOffset << " cm relative to the HepMC record in the lab frame, correcting the true vertices" << std::endl;
            }
          }
        }
      }
      else if (m_truthInfo)
      {
        PHG4TruthInfoContainer::ConstRange range = m_truthInfo->GetParticleRange();

        for (PHG4TruthInfoContainer::ConstIterator iter = range.first; iter != range.second; ++iter)
        {
          PHG4Particle *candidate = iter->second;

          PHG4Particle *motherG4 = nullptr;
          if (candidate->get_parent_id() != 0)
          {
            motherG4 = m_truthInfo->GetParticle(candidate->get_parent_id());
          }
          if (!motherG4)
          {
            continue;
          }

          if (motherG4->get_pid() == decay[0].second && motherG4->get_barcode() == decay[0].first.second && candidate->get_pid() == decay[i].second && candidate->get_barcode() == decay[i].first.second && m_nDaughters == 2)
          {
            daughterG4 = candidate;
            fillFromGeant4(motherG4, daughterG4, decay[i].second, index);
            haveTruthKinematics = true;
            break;
          }

          if (m_nDaughters == 3)
          {
            if (i != 4 && motherG4->get_pid() == decay[3].second && motherG4->get_barcode() == decay[3].first.second && candidate->get_pid() == decay[i].second && candidate->get_barcode() == decay[i].first.second)
            {
              PHG4Particle *motherG4_temp = nullptr;
              if (motherG4->get_parent_id() != 0)
              {
                motherG4_temp = m_truthInfo->GetParticle(motherG4->get_parent_id());
              }

              if (motherG4_temp && motherG4_temp->get_pid() == decay[0].second && motherG4_temp->get_barcode() == decay[0].first.second)
              {
                motherG4 = motherG4_temp;
              }
              else
              {
                continue;
              }
            }
            else if (i != 4)
            {
              continue;
            }

            if (i == 4 && !(motherG4->get_pid() == decay[0].second && motherG4->get_barcode() == decay[0].first.second && candidate->get_pid() == decay[i].second && candidate->get_barcode() == decay[i].first.second))
            {
              continue;
            }

            daughterG4 = candidate;
            fillFromGeant4(motherG4, daughterG4, decay[i].second, index);
            haveTruthKinematics = true;
            break;
          }
        }
      }

      if (!daughterG4)
      {
        ++m_counter_unmatchedTruth;
        if (Verbosity() > 0)
        {
          std::cout << __FILE__ << ": No Geant4 particle found for daughter " << index + 1 << " (PDG " << decay[i].second << ")" << std::endl;
        }
      }

      if (!haveTruthKinematics)
      {
        continue;  // Nothing to compare a reco track to
      }

      m_true_track_pT[index] = (float) daughterTrueLV.perp();
      m_true_track_eta[index] = (float) daughterTrueLV.pseudoRapidity();
      m_true_track_rapidity[index] = (float) daughterTrueLV.rapidity();
      m_true_track_phi[index] = (float) daughterTrueLV.phi();
      m_min_true_track_pT = std::min(m_true_track_pT[index], m_min_true_track_pT);
      m_max_true_track_pT = std::max(m_true_track_pT[index], m_max_true_track_pT);

      if (trackeval && daughterG4)
      {
        m_dst_track = trackeval->best_track_from(daughterG4);
        if (m_dst_track)
        {
          m_used_truth_reco_map[index] = true;
          recoTrackFound = true;
        }
      }
      else
      {
        for (auto &iter : *m_input_trackMap)
        {
          m_dst_track = iter.second;
          float delta_px = (m_dst_track->get_px() - daughterTrueLV.px()) / daughterTrueLV.px();
          float delta_py = (m_dst_track->get_py() - daughterTrueLV.py()) / daughterTrueLV.py();
          float delta_pz = (m_dst_track->get_pz() - daughterTrueLV.pz()) / daughterTrueLV.pz();

          if (std::abs(delta_px) <= m_truthRecoMatchPercent && std::abs(delta_py) <= m_truthRecoMatchPercent && std::abs(delta_pz) <= m_truthRecoMatchPercent)
          {
            recoTrackFound = true;
            break;
          }
        }
      }

      if (recoTrackFound)
      {
        selectedTracks.insert(m_dst_track);
        if (Verbosity() >= VERBOSITY_MORE)
        {
          m_dst_track->identify();
        }
        m_reco_track_exists[index] = true;
        m_reco_track_pT[index] = m_dst_track->get_pt();
        m_reco_track_eta[index] = m_dst_track->get_eta();
        m_reco_track_phi[index] = m_dst_track->get_phi();
        m_reco_track_chi2nDoF[index] = m_dst_track->get_chisq() / m_dst_track->get_ndf();

        countClusters(m_dst_track, daughterG4, index);

        m_min_reco_track_pT = std::min(m_reco_track_pT[index], m_min_reco_track_pT);
        m_max_reco_track_pT = std::max(m_reco_track_pT[index], m_max_reco_track_pT);

        CLHEP::HepLorentzVector daughterRecoLV;
        daughterRecoLV.setVectM(CLHEP::Hep3Vector(m_dst_track->get_px(), m_dst_track->get_py(), m_dst_track->get_pz()), getParticleMass(m_true_track_PID[index]));

        motherRecoLV += daughterRecoLV;
      }

      recoTrackFound = false;
    }
  }

  m_true_mother_mass = daughterSumTrueLV.m();

  if (haveVertices)
  {
    if (vertexOffsetKnown)
    {
      primaryVtx += vertexOffset;
      productionVtx += vertexOffset;
      decayVtx += vertexOffset;
    }

    m_primary_vtx_x = primaryVtx.x();
    m_primary_vtx_y = primaryVtx.y();
    m_primary_vtx_z = primaryVtx.z();
    m_production_vtx_x = productionVtx.x();
    m_production_vtx_y = productionVtx.y();
    m_production_vtx_z = productionVtx.z();
    m_decay_vtx_x = decayVtx.x();
    m_decay_vtx_y = decayVtx.y();
    m_decay_vtx_z = decayVtx.z();

    // Flight of the mother from where it was produced to where it decays. Not measured from the primary vertex,
    // that would include the flight of any parent (e.g. the B for a non-prompt D0)
    const CLHEP::Hep3Vector flight = decayVtx - productionVtx;
    m_true_decay_length = flight.mag();

    // Proper time of the mother, t = L m / (p c)
    if (m_true_mother_p > 0)
    {
      m_true_decay_time = m_true_decay_length * getParticleMass(decay[0].second) / (m_true_mother_p * c_cm_per_ps);
    }
  }

  bool foundDecay = true;
  if (selectedTracks.size() == m_nDaughters)
  {
    m_reco_mother_mass = motherRecoLV.m();
    m_reco_mother_pT = motherRecoLV.perp();
    if (m_write_track_map)
    {
      m_output_trackMap = findNode::getClass<SvtxTrackMap>(topNode, outputNodeName);
      for (const auto &track : selectedTracks)
      {
        m_output_trackMap->insertWithKey(track, track->get_id());
      }
    }
  }
  else
  {
    foundDecay = false;
  }

  selectedTracks.clear();

  return foundDecay;
}

void HFTrackEfficiency::initializeBranches()
{
  m_outfile = new TFile(m_outfile_name.c_str(), "RECREATE");
  m_tree = new TTree("HFTrackEfficiency", "HFTrackEfficiency");
  m_tree->OptimizeBaskets();
  m_tree->SetAutoSave(-5e6);  // Save the output file every 5MB

  m_tree->Branch("all_tracks_reconstructed", &m_all_tracks_reconstructed, "all_tracks_reconstructed/O");
  m_tree->Branch("is_primary", &m_is_primary, "is_primary/O");
  m_tree->Branch("true_mother_mass", &m_true_mother_mass, "true_mother_mass/F");
  m_tree->Branch("reco_mother_mass", &m_reco_mother_mass, "reco_mother_mass/F");
  m_tree->Branch("true_mother_pT", &m_true_mother_pT, "true_mother_pT/F");
  m_tree->Branch("reco_mother_pT", &m_reco_mother_pT, "reco_mother_pT/F");
  m_tree->Branch("true_mother_p", &m_true_mother_p, "true_mother_p/F");
  m_tree->Branch("true_mother_eta", &m_true_mother_eta, "true_mother_eta/F");
  m_tree->Branch("true_mother_rapidity", &m_true_mother_rapidity, "true_mother_rapidity/F");
  m_tree->Branch("true_mother_phi", &m_true_mother_phi, "true_mother_phi/F");
  m_tree->Branch("min_true_track_pT", &m_min_true_track_pT, "min_true_track_pT/F");
  m_tree->Branch("min_reco_track_pT", &m_min_reco_track_pT, "min_reco_track_pT/F");
  m_tree->Branch("max_true_track_pT", &m_max_true_track_pT, "max_true_track_pT/F");
  m_tree->Branch("max_reco_track_pT", &m_max_reco_track_pT, "max_reco_track_pT/F");

  const char *detectorNames[m_nDetectors] = {"mvtx", "intt", "tpc", "tpot"};

  for (unsigned int iTrack = 0; iTrack < m_nDaughters; ++iTrack)
  {
    std::string daughter_number = "track_" + std::to_string(iTrack + 1);
    m_tree->Branch("reco_" + TString(daughter_number) + "_exists", &m_reco_track_exists[iTrack], "reco_" + TString(daughter_number) + "_exists/O");
    m_tree->Branch("reco_" + TString(daughter_number) + "_used_truth_reco_map", &m_used_truth_reco_map[iTrack], "reco_" + TString(daughter_number) + "_used_truth_reco_map/O");
    m_tree->Branch("true_" + TString(daughter_number) + "_pT", &m_true_track_pT[iTrack], "true_" + TString(daughter_number) + "_pT/F");
    m_tree->Branch("reco_" + TString(daughter_number) + "_pT", &m_reco_track_pT[iTrack], "reco_" + TString(daughter_number) + "_pT/F");
    m_tree->Branch("true_" + TString(daughter_number) + "_eta", &m_true_track_eta[iTrack], "true_" + TString(daughter_number) + "_eta/F");
    m_tree->Branch("reco_" + TString(daughter_number) + "_eta", &m_reco_track_eta[iTrack], "reco_" + TString(daughter_number) + "_eta/F");
    m_tree->Branch("true_" + TString(daughter_number) + "_rapidity", &m_true_track_rapidity[iTrack], "true_" + TString(daughter_number) + "_rapidity/F");
    m_tree->Branch("true_" + TString(daughter_number) + "_phi", &m_true_track_phi[iTrack], "true_" + TString(daughter_number) + "_phi/F");
    m_tree->Branch("reco_" + TString(daughter_number) + "_phi", &m_reco_track_phi[iTrack], "reco_" + TString(daughter_number) + "_phi/F");
    m_tree->Branch("true_" + TString(daughter_number) + "_PID", &m_true_track_PID[iTrack], "true_" + TString(daughter_number) + "_PID/F");
    m_tree->Branch("reco_" + TString(daughter_number) + "_chi2nDoF", &m_reco_track_chi2nDoF[iTrack], "reco_" + TString(daughter_number) + "_chi2nDoF/F");
    for (int iDet = 0; iDet < m_nDetectors; ++iDet)
    {
      const TString nClusters = "reco_" + TString(daughter_number) + "_nclusters_" + detectorNames[iDet];
      const TString nTrueClusters = "reco_" + TString(daughter_number) + "_ntrue_clusters_" + detectorNames[iDet];
      m_tree->Branch(nClusters, &m_reco_track_nclusters[iTrack][iDet], nClusters + "/I");
      m_tree->Branch(nTrueClusters, &m_reco_track_ntrue_clusters[iTrack][iDet], nTrueClusters + "/I");
    }
  }

  m_tree->Branch("true_primary_vertex_x", &m_primary_vtx_x, "true_primary_vertex_x/F");
  m_tree->Branch("true_primary_vertex_y", &m_primary_vtx_y, "true_primary_vertex_y/F");
  m_tree->Branch("true_primary_vertex_z", &m_primary_vtx_z, "true_primary_vertex_z/F");
  m_tree->Branch("true_production_vertex_x", &m_production_vtx_x, "true_production_vertex_x/F");
  m_tree->Branch("true_production_vertex_y", &m_production_vtx_y, "true_production_vertex_y/F");
  m_tree->Branch("true_production_vertex_z", &m_production_vtx_z, "true_production_vertex_z/F");
  m_tree->Branch("true_decay_vertex_x", &m_decay_vtx_x, "true_decay_vertex_x/F");
  m_tree->Branch("true_decay_vertex_y", &m_decay_vtx_y, "true_decay_vertex_y/F");
  m_tree->Branch("true_decay_vertex_z", &m_decay_vtx_z, "true_decay_vertex_z/F");
  m_tree->Branch("true_decay_length", &m_true_decay_length, "true_decay_length/F");
  m_tree->Branch("true_decay_time", &m_true_decay_time, "true_decay_time/F");
}

void HFTrackEfficiency::resetBranches()
{
  m_all_tracks_reconstructed = false;
  m_is_primary = false;
  m_true_mother_mass = std::numeric_limits<float>::quiet_NaN();
  m_reco_mother_mass = std::numeric_limits<float>::quiet_NaN();
  m_true_mother_pT = std::numeric_limits<float>::quiet_NaN();
  m_reco_mother_pT = std::numeric_limits<float>::quiet_NaN();
  m_true_mother_p = std::numeric_limits<float>::quiet_NaN();
  m_true_mother_eta = std::numeric_limits<float>::quiet_NaN();
  m_true_mother_rapidity = std::numeric_limits<float>::quiet_NaN();
  m_true_mother_phi = std::numeric_limits<float>::quiet_NaN();
  m_min_true_track_pT = std::numeric_limits<float>::max();
  m_min_reco_track_pT = std::numeric_limits<float>::max();
  m_max_true_track_pT = std::numeric_limits<float>::lowest();
  m_max_reco_track_pT = std::numeric_limits<float>::lowest();
  for (int iTrack = 0; iTrack < m_maxTracks; ++iTrack)
  {
    m_reco_track_exists[iTrack] = false;
    m_used_truth_reco_map[iTrack] = false;
    m_true_track_pT[iTrack] = std::numeric_limits<float>::quiet_NaN();
    m_reco_track_pT[iTrack] = std::numeric_limits<float>::quiet_NaN();
    m_true_track_eta[iTrack] = std::numeric_limits<float>::quiet_NaN();
    m_reco_track_eta[iTrack] = std::numeric_limits<float>::quiet_NaN();
    m_true_track_rapidity[iTrack] = std::numeric_limits<float>::quiet_NaN();
    m_true_track_phi[iTrack] = std::numeric_limits<float>::quiet_NaN();
    m_reco_track_phi[iTrack] = std::numeric_limits<float>::quiet_NaN();
    m_true_track_PID[iTrack] = std::numeric_limits<float>::quiet_NaN();
    m_reco_track_chi2nDoF[iTrack] = std::numeric_limits<float>::quiet_NaN();
    for (int iDet = 0; iDet < m_nDetectors; ++iDet)
    {
      m_reco_track_nclusters[iTrack][iDet] = -1;
      m_reco_track_ntrue_clusters[iTrack][iDet] = -1;
    }
  }

  m_primary_vtx_x = std::numeric_limits<float>::quiet_NaN();
  m_primary_vtx_y = std::numeric_limits<float>::quiet_NaN();
  m_primary_vtx_z = std::numeric_limits<float>::quiet_NaN();
  m_production_vtx_x = std::numeric_limits<float>::quiet_NaN();
  m_production_vtx_y = std::numeric_limits<float>::quiet_NaN();
  m_production_vtx_z = std::numeric_limits<float>::quiet_NaN();
  m_decay_vtx_x = std::numeric_limits<float>::quiet_NaN();
  m_decay_vtx_y = std::numeric_limits<float>::quiet_NaN();
  m_decay_vtx_z = std::numeric_limits<float>::quiet_NaN();
  m_true_decay_length = std::numeric_limits<float>::quiet_NaN();
  m_true_decay_time = std::numeric_limits<float>::quiet_NaN();
}

void HFTrackEfficiency::getDecayDescriptor()
{
  Decay decay = m_decayMap->begin()->second;
  m_decay_descriptor = getParticleName(decay[0].second);
  m_decay_descriptor += " ->";
  for (unsigned int i = 1; i < decay.size(); ++i)
  {
    m_decay_descriptor += " ";
    m_decay_descriptor += getParticleName(decay[i].second);
  }
}

void HFTrackEfficiency::getNDaughters()
{
  int trackableParticles[] = {11, 13, 211, 321, 2212};
  m_nDaughters = 0;
  Decay decay = m_decayMap->begin()->second;
  for (unsigned int i = 1; i < decay.size(); ++i)
  {
    if (std::find(std::begin(trackableParticles), std::end(trackableParticles),
                  std::abs(decay[i].second)) != std::end(trackableParticles))
    {
      ++m_nDaughters;
    }
  }
}

std::string HFTrackEfficiency::getParticleName(const int PDGID)
{
  return TDatabasePDG::Instance()->GetParticle(PDGID)->GetName();
}

float HFTrackEfficiency::getParticleMass(const int PDGID)
{
  return TDatabasePDG::Instance()->GetParticle(PDGID)->Mass();
}

//____________________________________________________________________________..
void HFTrackEfficiency::PrintEff()
{
  std::cout << "\n--------------- Heavy Flavor Tracking Efficiency ---------------" << std::endl;
  std::cout << "Tracking Efficiency for " << m_decay_descriptor << std::endl;
  std::cout << "Number of decays fully in acceptance: " << m_counter_allDecays << std::endl;
  std::cout << "Number of decays with all tracks reconstructed: " << m_counter_acceptedDecays << std::endl;
  if (m_counter_allDecays != 0)
  {
    std::cout << "Efficiency: " << std::setw(4) << (float) m_counter_acceptedDecays * 100 / float(m_counter_allDecays) << "%" << std::endl;
  }
  if (m_counter_unmatchedTruth != 0)
  {
    std::cout << "Warning: " << m_counter_unmatchedTruth << " daughters had no Geant4 particle, these use the momentum matching to find a reco track" << std::endl;
  }
  std::cout << "----------------------------------------------------------------\n"
            << std::endl;
}
