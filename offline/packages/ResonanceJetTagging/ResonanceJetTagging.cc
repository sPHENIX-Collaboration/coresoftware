#include "ResonanceJetTagging.h"

/// Cluster/Calorimeter includes
#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawClusterUtility.h>

/// Jet includes
#include <g4jets/Jet.h>
#include <g4jets/JetMapv1.h>
#include <g4jets/Jetv1.h>

/// Tracking includes
#include <g4vertex/GlobalVertex.h>
#include <g4vertex/GlobalVertexMap.h>

#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>

/// HEPMC truth includes
#include <phhepmc/PHHepMCGenEvent.h>
#include <phhepmc/PHHepMCGenEventMap.h>

#include <g4main/PHG4Particle.h>            // for PHG4Particle
#include <g4main/PHG4Particlev2.h>            // for PHG4Particle

// Particle Flow
#include <particleflowreco/ParticleFlowElement.h>
#include <particleflowreco/ParticleFlowElementContainer.h>

/// Fun4All includes
#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>
#include <phool/phool.h>

#include <kfparticle_sphenix/KFParticle_Container.h>

#include <KFParticle.h>

#include <fastjet/ClusterSequence.hh>
#include <fastjet/FunctionOfPseudoJet.hh>
#include <fastjet/JetDefinition.hh>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <HepMC/GenEvent.h>
#include <HepMC/GenVertex.h>
#pragma GCC diagnostic pop

/// ROOT includes
#include <TDatabasePDG.h>

/// C++ includes
#include <algorithm>  // for max
#include <cassert>
#include <cmath>
#include <cstddef>   // for size_t
#include <iostream>  // for operator<<
#include <iterator>  // for reverse_i...
#include <memory>    // for allocator...
#include <set>       // for set
#include <string>

#include <g4main/PHG4TruthInfoContainer.h>
class PHG4TruthInfoContainer;

/**
 * ResonanceJetTagging is a class developed to reconstruct jets containing a D-meson
 * The class can be adapted to tag jets using any kind of particle
 * Author: Antonio Silva (antonio.sphenix@gmail.com)
 * Contributor: Jakub Kvapil (jakub.kvapil@cern.ch)
 */

/**
 * Constructor of module
 */
ResonanceJetTagging::ResonanceJetTagging(const std::string &name, const TAG tag, const std::string &KFparticle_Container_name)
  : SubsysReco(name)
  , m_particleflow_mineta(-1.1)
  , m_particleflow_maxeta(1.1)
  , m_track_minpt(0.)
  , m_track_maxpt(9999.)
  , m_track_mineta(-1.1)
  , m_track_maxeta(1.1)
  , m_EMCal_cluster_minpt(0.)
  , m_EMCal_cluster_maxpt(9999.)
  , m_EMCal_cluster_mineta(-1.1)
  , m_EMCal_cluster_maxeta(1.1)
  , m_HCal_cluster_minpt(0.)
  , m_HCal_cluster_maxpt(9999.)
  , m_HCal_cluster_mineta(-1.1)
  , m_HCal_cluster_maxeta(1.1)
  , m_add_particleflow(true)
  , m_add_tracks(false)
  , m_add_EMCal_clusters(false)
  , m_add_HCal_clusters(false)
  , m_jetr(0.4)
  , m_jetalgo(fastjet::antikt_algorithm)
  , m_recomb_scheme(fastjet::pt_scheme)
  , m_dorec(true)
  , m_dotruth(false)
  , m_nDaughters(0)
  , m_tag_particle(tag)
  , m_KFparticle_name(KFparticle_Container_name)
{
  switch (m_tag_particle) {
    case ResonanceJetTagging::TAG::D0:
      m_tag_pdg = 421;
      m_nDaughters = 2;
      break;
    case ResonanceJetTagging::TAG::D0TOK3PI:
      m_tag_pdg = 421;
      m_nDaughters = 4;
      break;
    case ResonanceJetTagging::TAG::DPLUS:
      m_tag_pdg = 411;
      m_nDaughters = 3;
      break;
    case ResonanceJetTagging::TAG::DSTAR:
      m_tag_pdg = 413;
      m_nDaughters = 0;
      break;
    case ResonanceJetTagging::TAG::JPSY:
      m_tag_pdg = 433;
      m_nDaughters = 0;
      break;
    case ResonanceJetTagging::TAG::K0:
      m_tag_pdg = 311;
      m_nDaughters = 0;
      break;
    case ResonanceJetTagging::TAG::GAMMA:
      m_tag_pdg = 22;
      m_nDaughters = 0;
      break;
    case ResonanceJetTagging::TAG::ELECTRON:
      m_tag_pdg = 11;
      m_nDaughters = 0;
      break;
    case ResonanceJetTagging::TAG::LAMBDAC:
      m_tag_pdg = 4122;
      m_nDaughters = 3;
      break;
  }

}

/**
 * Destructor of module
 */
ResonanceJetTagging::~ResonanceJetTagging()
{

}

/**
 * Initialize the module and prepare looping over events
 */
int ResonanceJetTagging::Init(PHCompositeNode *topNode)
{
  if (Verbosity() > 5)
  {
    std::cout << "Beginning Init in ResonanceJetTagging" << std::endl;
  }

  createJetNode(topNode);

  return 0;
}

/**
 * Main workhorse function where each event is looped over and
 * data from each event is collected from the node tree for analysis
 */
int ResonanceJetTagging::process_event(PHCompositeNode *topNode)
{
  if(m_nDaughters == 0)
  {
    std::cout<<"ERROR: Number of Decay Daughters Not Set, ABORTING!";
    return Fun4AllReturnCodes::ABORTRUN;
  }

  switch (m_tag_particle) {
    case ResonanceJetTagging::TAG::D0:
      [[fallthrough]];
    case ResonanceJetTagging::TAG::D0TOK3PI:
      [[fallthrough]];
    case ResonanceJetTagging::TAG::DPLUS:
      [[fallthrough]];
    case ResonanceJetTagging::TAG::LAMBDAC:
      return tagHFHadronic(topNode);
      break;
    default:
      std::cout<<"ERROR: Tag Function Not Set, ABORTING!";
      return Fun4AllReturnCodes::ABORTRUN;
      break;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

/**
 * End the module and finish any data collection. Clean up any remaining
 * loose ends
 */
int ResonanceJetTagging::End(PHCompositeNode * /*topNode*/)
{
  if (Verbosity() > 1)
  {
    std::cout << "Finished ResonanceJetTagging analysis package" << std::endl;
  }

  return 0;
}

int ResonanceJetTagging::tagHFHadronic(PHCompositeNode *topNode)
{
  if(m_dorec)
  {
    KFParticle_Container *kfContainer = findNode::getClass<KFParticle_Container>(topNode, m_KFparticle_name);
    if(!kfContainer) return Fun4AllReturnCodes::ABORTEVENT;

    KFParticle *TagCand = nullptr;
    PHG4Particlev2 *Cand = new PHG4Particlev2();
    //const int nDaughters = 2;  // TagCand->NDaughters() is returning 0, bug?
    std::vector<KFParticle*> TagDaughters(m_nDaughters);
    std::vector<PHG4Particlev2*> Daughters(m_nDaughters);
    m_jet_id = 0;

    for (unsigned int i = 0; i < kfContainer->size(); i++)
    {
      TagCand = kfContainer->get(i);
      if (std::abs(TagCand->GetPDG()) == m_tag_pdg)
      {
        for (int idau = 0; idau < m_nDaughters; idau++)
        {
          TagDaughters.at(idau) = kfContainer->get(i + idau + 1);
          Daughters.at(idau) = new PHG4Particlev2();
          Daughters.at(idau)->set_px(TagDaughters.at(idau)->Px());
          Daughters.at(idau)->set_py(TagDaughters.at(idau)->Py());
          Daughters.at(idau)->set_pz(TagDaughters.at(idau)->Pz());
          //For daughters keep ID as the track ID, so they can be removed from the sample
          //given to fastjet
          Daughters.at(idau)->set_barcode(TagDaughters.at(idau)->Id());
        }

        Cand->set_px(TagCand->Px());
        Cand->set_py(TagCand->Py());
        Cand->set_pz(TagCand->Pz());
        Cand->set_e(TagCand->E());
        //For tag particle keep ID as position in the KF Container
        //so that later the tag particle can be recovered from jet constituent
        Cand->set_barcode(i);

        findTaggedJets(topNode, Cand, Daughters);

        for (long unsigned int idau = 0; idau < Daughters.size(); idau++) delete Daughters.at(idau);

        m_jet_id++;
        i += m_nDaughters;  // Go to the next D meson
      }
    }
  delete Cand;
  }

  if (m_dotruth)
  {
    findMCTaggedJets(topNode);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

void ResonanceJetTagging::findTaggedJets(PHCompositeNode *topNode, PHG4Particlev2 *Tag, const std::vector<PHG4Particlev2*> &TagDecays)
{
  std::unique_ptr<fastjet::JetDefinition> jetdef(new fastjet::JetDefinition(m_jetalgo, m_jetr, m_recomb_scheme, fastjet::Best));
  std::vector<fastjet::PseudoJet> particles;
  Jet *taggedJet;
  std::map<int, std::pair<Jet::SRC, int>> fjMap;

  fastjet::PseudoJet fjTag(Tag->get_px(), Tag->get_py(), Tag->get_pz(), Tag->get_e());
  fjTag.set_user_index(0);  // index 0 is the tag particle
  particles.push_back(fjTag);
  fjMap.insert(std::pair<int, std::pair<Jet::SRC, int>>(0, std::make_pair(Jet::SRC::VOID, Tag->get_barcode())));  // Maybe we could have a Jet::SRC::HF for HF-tagging?

  if (m_add_particleflow)
  {
    addParticleFlow(topNode, particles, TagDecays, fjMap);
  }

  if (m_add_tracks)
  {
    addTracks(topNode, particles, TagDecays, fjMap);
  }

  if (m_add_EMCal_clusters || m_add_HCal_clusters)
  {
    addClusters(topNode, particles, fjMap);
  }

  fastjet::ClusterSequence jetFinder(particles, *jetdef);
  std::vector<fastjet::PseudoJet> fastjets = jetFinder.inclusive_jets();

  taggedJet = new Jetv1();

  for (auto &fastjet : fastjets)
  {
    bool isTaggedJet = false;
    std::vector<fastjet::PseudoJet> comps = fastjet.constituents();
    for (auto &comp : comps)
    {
      taggedJet->insert_comp(fjMap[comp.user_index()].first, fjMap[comp.user_index()].second);
      if (comp.user_index() == 0)
      {
        taggedJet->set_px(fastjet.px());
        taggedJet->set_py(fastjet.py());
        taggedJet->set_pz(fastjet.pz());
        taggedJet->set_e(fastjet.e());
        taggedJet->set_id(m_jet_id);
        isTaggedJet = true;
      }
    }
    if (isTaggedJet)
    {
      break;
    }
    else
    {
      taggedJet->clear_comp();
    }
  }
  m_taggedJetMap->insert(taggedJet);

  return;
}

void ResonanceJetTagging::addParticleFlow(PHCompositeNode *topNode, std::vector<fastjet::PseudoJet> &particles, const std::vector<PHG4Particlev2*> &TagDecays, std::map<int, std::pair<Jet::SRC, int>> &fjMap)
{
  ParticleFlowElementContainer *pflowContainer = findNode::getClass<ParticleFlowElementContainer>(topNode, "ParticleFlowElements");

  if(!pflowContainer)
  {
    std::cout << PHWHERE
         << "ParticleFlowElements node is missing, can't collect particles"
         << std::endl;
    return;
  }

  SvtxTrack *track = nullptr;

  std::vector<RawCluster *> pfEMCalClusters;
  std::vector<RawCluster *> pfHCalClusters;

  int idpart = particles.size();

  ParticleFlowElementContainer::ConstRange begin_end = pflowContainer->getParticleFlowElements();
  ParticleFlowElementContainer::ConstIterator rtiter;
  for (rtiter = begin_end.first; rtiter != begin_end.second; ++rtiter)
  {
    ParticleFlowElement *pflow = rtiter->second;

    if (!pflow)
    {
      continue;
    }

    if (!isAcceptableParticleFlow(pflow))
    {
      continue;
    }

    track = pflow->get_track();

    // Remove D0 decay daughter
    if (track && isDecay(track, TagDecays))
    {
      continue;
    }

    fastjet::PseudoJet fjPartFlow(pflow->get_px(), pflow->get_py(), pflow->get_pz(), pflow->get_e());

    fjPartFlow.set_user_index(idpart);
    particles.push_back(fjPartFlow);
    fjMap.insert(std::pair<int, std::pair<Jet::SRC, int>>(idpart, std::make_pair(Jet::SRC::PARTICLE, pflow->get_id())));
    idpart++;
  }
}

bool ResonanceJetTagging::isAcceptableParticleFlow(ParticleFlowElement *pfPart)
{
  // Only eta cut at this moment
  if ((pfPart->get_eta() < m_particleflow_mineta) || (pfPart->get_eta() > m_particleflow_maxeta))
  {
    return false;
  }

  return true;
}
void ResonanceJetTagging::addTracks(PHCompositeNode *topNode, std::vector<fastjet::PseudoJet> &particles, const std::vector<PHG4Particlev2*> &TagDecays, std::map<int, std::pair<Jet::SRC, int>> &fjMap)
{
  SvtxTrackMap *trackmap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");

  if (!trackmap)
  {
    std::cout << PHWHERE
              << "SvtxTrackMap node is missing, can't collect tracks"
              << std::endl;
    return;
  }

  int idpart = particles.size();

  SvtxTrack *track = nullptr;

  for (auto &iter : *trackmap)
  {
    track = iter.second;

    if (!isAcceptableTrack(track))
    {
      continue;
    }

    if (isDecay(track, TagDecays))
    {
      continue;
    }

    fastjet::PseudoJet fjTrack(track->get_px(), track->get_py(), track->get_pz(), 0.);

    fjTrack.set_user_index(idpart);
    particles.push_back(fjTrack);
    fjMap.insert(std::pair<int, std::pair<Jet::SRC, int>>(idpart, std::make_pair(Jet::SRC::TRACK, track->get_id())));
    idpart++;
  }
}

bool ResonanceJetTagging::isAcceptableTrack(SvtxTrack *track)
{
  if ((track->get_pt() < m_track_minpt) || (track->get_pt() > m_track_maxpt))
  {
    return false;
  }
  if ((track->get_eta() < m_track_mineta) || (track->get_eta() > m_track_maxeta))
  {
    return false;
  }

  return true;
}

void ResonanceJetTagging::addClusters(PHCompositeNode *topNode, std::vector<fastjet::PseudoJet> &particles, std::map<int, std::pair<Jet::SRC, int>> &fjMap)
{
  GlobalVertexMap *vertexmap = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");
  if (!vertexmap)
  {
    std::cout << "ResonanceJetTagging::getEmcalClusters - Fatal Error - GlobalVertexMap node is missing. Please turn on the do_global flag in the main macro in order to reconstruct the global vertex." << std::endl;
    assert(vertexmap);  // force quit

    return;
  }

  if (vertexmap->empty())
  {
    std::cout << "ResonanceJetTagging::getEmcalClusters - Fatal Error - GlobalVertexMap node is empty. Please turn on the do_global flag in the main macro in order to reconstruct the global vertex." << std::endl;
    return;
  }

  GlobalVertex *vtx = vertexmap->begin()->second;
  if (vtx == nullptr)
  {
    return;
  }

  // Get the current position in the particles vector
  int idpart = particles.size();

  // EMCAL Clusters
  if (m_add_EMCal_clusters)
  {
    RawClusterContainer *clustersEMC = findNode::getClass<RawClusterContainer>(topNode, "CLUSTER_CEMC");

    if (!clustersEMC)
    {
      std::cout << PHWHERE
                << "EMCal cluster node is missing, can't collect EMCal clusters"
                << std::endl;
      return;
    }

    RawClusterContainer::ConstRange begin_end_EMC = clustersEMC->getClusters();
    RawClusterContainer::ConstIterator clusIter_EMC;

    /// Loop over the EMCal clusters
    for (clusIter_EMC = begin_end_EMC.first;
         clusIter_EMC != begin_end_EMC.second;
         ++clusIter_EMC)
    {
      const RawCluster *cluster = clusIter_EMC->second;

      CLHEP::Hep3Vector vertex(vtx->get_x(), vtx->get_y(), vtx->get_z());
      CLHEP::Hep3Vector E_vec_cluster = RawClusterUtility::GetECoreVec(*cluster, vertex);

      if (!isAcceptableEMCalCluster(E_vec_cluster))
      {
        continue;
      }

      double cluster_e = E_vec_cluster.mag();

      double cluster_pt = E_vec_cluster.perp();

      double cluster_phi = E_vec_cluster.getPhi();

      double cluster_px = cluster_pt * cos(cluster_phi);
      double cluster_py = cluster_pt * sin(cluster_phi);
      double cluster_pz = sqrt(cluster_e * cluster_e - cluster_px * cluster_px - cluster_py * cluster_py);

      fastjet::PseudoJet fjCluster(cluster_px, cluster_py, cluster_pz, cluster_e);

      fjCluster.set_user_index(idpart);
      particles.push_back(fjCluster);
      fjMap.insert(std::pair<int, std::pair<Jet::SRC, int>>(idpart, std::make_pair(Jet::SRC::CEMC_CLUSTER, cluster->get_id())));
      idpart++;
    }
  }

  if (m_add_HCal_clusters)
  {
    // HCAL Clusters
    RawClusterContainer *clustersHCALIN = findNode::getClass<RawClusterContainer>(topNode, "CLUSTER_HCALIN");

    if (!clustersHCALIN)
    {
      std::cout << PHWHERE
                << "EMCal cluster node is missing, can't collect EMCal clusters"
                << std::endl;
      return;
    }

    RawClusterContainer::ConstRange begin_end_HCALIN = clustersHCALIN->getClusters();
    RawClusterContainer::ConstIterator clusIter_HCALIN;

    /// Loop over the EMCal clusters
    for (clusIter_HCALIN = begin_end_HCALIN.first;
         clusIter_HCALIN != begin_end_HCALIN.second;
         ++clusIter_HCALIN)
    {
      /// Get this cluster
      const RawCluster *cluster = clusIter_HCALIN->second;

      CLHEP::Hep3Vector vertex(vtx->get_x(), vtx->get_y(), vtx->get_z());
      CLHEP::Hep3Vector E_vec_cluster = RawClusterUtility::GetECoreVec(*cluster, vertex);

      if (!isAcceptableHCalCluster(E_vec_cluster))
      {
        continue;
      }

      double cluster_e = E_vec_cluster.mag();
      double cluster_pt = E_vec_cluster.perp();
      double cluster_phi = E_vec_cluster.getPhi();

      double cluster_px = cluster_pt * cos(cluster_phi);
      double cluster_py = cluster_pt * sin(cluster_phi);
      double cluster_pz = sqrt(cluster_e * cluster_e - cluster_px * cluster_px - cluster_py * cluster_py);

      fastjet::PseudoJet fjCluster(cluster_px, cluster_py, cluster_pz, cluster_e);

      fjCluster.set_user_index(idpart);
      particles.push_back(fjCluster);
      idpart++;
    }

    RawClusterContainer *clustersHCALOUT = findNode::getClass<RawClusterContainer>(topNode, "CLUSTER_HCALOUT");

    if (!clustersHCALOUT)
    {
      std::cout << PHWHERE
                << "EMCal cluster node is missing, can't collect EMCal clusters"
                << std::endl;
      return;
    }

    RawClusterContainer::ConstRange begin_end_HCALOUT = clustersHCALOUT->getClusters();
    RawClusterContainer::ConstIterator clusIter_HCALOUT;

    /// Loop over the EMCal clusters
    for (clusIter_HCALOUT = begin_end_HCALOUT.first;
         clusIter_HCALOUT != begin_end_HCALOUT.second;
         ++clusIter_HCALOUT)
    {
      /// Get this cluster
      const RawCluster *cluster = clusIter_HCALOUT->second;

      /// Get cluster characteristics
      /// This helper class determines the photon characteristics
      /// depending on the vertex position
      /// This is important for e.g. eta determination and E_T determination
      CLHEP::Hep3Vector vertex(vtx->get_x(), vtx->get_y(), vtx->get_z());
      CLHEP::Hep3Vector E_vec_cluster = RawClusterUtility::GetECoreVec(*cluster, vertex);

      if (!isAcceptableHCalCluster(E_vec_cluster))
      {
        continue;
      }

      double cluster_e = E_vec_cluster.mag();
      double cluster_pt = E_vec_cluster.perp();
      double cluster_phi = E_vec_cluster.getPhi();

      double cluster_px = cluster_pt * cos(cluster_phi);
      double cluster_py = cluster_pt * sin(cluster_phi);
      double cluster_pz = sqrt(cluster_e * cluster_e - cluster_px * cluster_px - cluster_py * cluster_py);

      fastjet::PseudoJet fjCluster(cluster_px, cluster_py, cluster_pz, cluster_e);

      fjCluster.set_user_index(idpart);
      particles.push_back(fjCluster);
      idpart++;
    }
  }

}

bool ResonanceJetTagging::isAcceptableEMCalCluster(CLHEP::Hep3Vector &E_vec_cluster)
{
  if ((E_vec_cluster.perp() < m_EMCal_cluster_minpt) || (E_vec_cluster.perp() > m_EMCal_cluster_maxpt))
  {
    return false;
  }
  if ((E_vec_cluster.pseudoRapidity() < m_EMCal_cluster_mineta) || (E_vec_cluster.pseudoRapidity() > m_EMCal_cluster_maxeta))
  {
    return false;
  }

  return true;
}

bool ResonanceJetTagging::isAcceptableHCalCluster(CLHEP::Hep3Vector &E_vec_cluster)
{
  if ((E_vec_cluster.perp() < m_HCal_cluster_minpt) || (E_vec_cluster.perp() > m_HCal_cluster_maxpt))
  {
    return false;
  }
  if ((E_vec_cluster.pseudoRapidity() < m_HCal_cluster_mineta) || (E_vec_cluster.pseudoRapidity() > m_HCal_cluster_maxeta))
  {
    return false;
  }

  return true;
}

bool ResonanceJetTagging::isDecay(SvtxTrack *track,  const std::vector<PHG4Particlev2*> &decays)
{
  for (long unsigned int idecay = 0; idecay < decays.size(); idecay++)
  {
    if(int(track->get_id()) == decays.at(idecay)->get_barcode())
    {
      return true;
    }
  }
  return false;
}

bool ResonanceJetTagging::isDecay(HepMC::GenParticle *particle, const std::vector<PHG4Particlev2*> &decays)
{
  for (long unsigned int idecay = 0; idecay < decays.size(); idecay++)
  {
    if (particle->barcode() == decays.at(idecay)->get_barcode())
    {
      return true;
    }
  }
  return false;
}

void ResonanceJetTagging::findMCTaggedJets(PHCompositeNode *topNode)
{
  PHHepMCGenEventMap *hepmceventmap = findNode::getClass<PHHepMCGenEventMap>(topNode, "PHHepMCGenEventMap");
  /// If the node was not properly put on the tree, return
  if (!hepmceventmap)
  {
    std::cout << PHWHERE
              << "HEPMC event map node is missing, can't collected HEPMC truth particles"
              << std::endl;
    return;
  }

  PHHepMCGenEvent *hepmcevent = hepmceventmap->get(1);

  if (!hepmcevent)
  {
    return;
  }

  HepMC::GenEvent *hepMCGenEvent = hepmcevent->getEvent();

  if (!hepMCGenEvent)
  {
    return;
  }

  TDatabasePDG *database = TDatabasePDG::Instance();
  TParticlePDG *partPDG = nullptr;

  PHG4TruthInfoContainer *m_truthinfo = nullptr;

  m_truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");

  m_truth_jet_id = 0;

  Jet *mcTaggedJet;

  for (HepMC::GenEvent::particle_const_iterator tag = hepMCGenEvent->particles_begin(); tag != hepMCGenEvent->particles_end(); ++tag)
  {
    if (((*tag)->momentum().eta() < m_track_mineta) || ((*tag)->momentum().eta() > m_track_maxeta))
    {
      continue;
    }

    if (std::abs((*tag)->pdg_id()) == m_tag_pdg)
    {
      std::vector<int> decayIDs;

      HepMC::GenVertex *TagVertex = (*tag)->end_vertex();
      //if HEPMC vertex exists
      if(TagVertex){
        for (HepMC::GenVertex::particle_iterator it = TagVertex->particles_begin(HepMC::descendants); it != TagVertex->particles_end(HepMC::descendants); ++it)
        {
          decayIDs.push_back((*it)->barcode());
        }
      //if not, look into GEANT
      } 
      else
      {
        PHG4TruthInfoContainer::ConstRange range = m_truthinfo->GetParticleRange();
        for(PHG4TruthInfoContainer::ConstIterator iter = range.first; iter != range.second; ++iter)
        {  
          PHG4Particle* g4particle = iter->second;
          PHG4Particle* mother = nullptr;
          if (g4particle->get_parent_id() != 0) mother = m_truthinfo->GetParticle(g4particle->get_parent_id());
          else continue;
          if (mother->get_barcode() == (*tag)->barcode() && mother->get_pid() == (*tag)->pdg_id())
          {
            decayIDs.push_back(g4particle->get_barcode());
          }
        }
      }

      std::map<int, std::pair<Jet::SRC, int>> fjMapMC;

      std::vector<fastjet::PseudoJet> particles;
      fastjet::PseudoJet fjMCParticle((*tag)->momentum().px(), (*tag)->momentum().py(), (*tag)->momentum().pz(), (*tag)->momentum().e());
      fjMapMC.insert(std::pair<int, std::pair<Jet::SRC, int>>(0, std::make_pair(Jet::SRC::VOID, (*tag)->barcode())));  // Maybe we could have a Jet::SRC::HF for HF-tagging? Or maybe Jet::SRC::TAG? Using VOID for the tag particle
      fjMCParticle.set_user_index(0);
      particles.push_back(fjMCParticle);

      int idpart = particles.size();

      for (HepMC::GenEvent::particle_const_iterator p = hepMCGenEvent->particles_begin(); p != hepMCGenEvent->particles_end(); ++p)
      {
        if (((*p)->momentum().eta() < m_track_mineta) || ((*p)->momentum().eta() > m_track_maxeta))
        {
          continue;
        }
        if ((*p)->status() > 1)
        {
          continue;
        }

        partPDG = database->GetParticle((*p)->pdg_id());
        double hepmcPartPt = std::sqrt(((*p)->momentum().px() * (*p)->momentum().px()) + ((*p)->momentum().py() * (*p)->momentum().py()));

        if (partPDG->Charge() != 0)
        {
          if ((hepmcPartPt < m_track_minpt) || (hepmcPartPt > m_track_maxpt))
          {
            continue;
          }
        }
        else
        {
          if ((hepmcPartPt < m_EMCal_cluster_minpt) || (hepmcPartPt > m_EMCal_cluster_maxpt))
          {
            continue;
          }
        }
        bool isTagDecay = false;
        for (auto mcid : decayIDs)
        {
          if (mcid == (*p)->barcode())
          {
            isTagDecay = true;
          }
        }
        if (isTagDecay)
        {
          continue;
        }

        fastjet::PseudoJet fjMCParticle2((*p)->momentum().px(), (*p)->momentum().py(), (*p)->momentum().pz(), (*p)->momentum().e());
        fjMapMC.insert(std::pair<int, std::pair<Jet::SRC, int>>(idpart, std::make_pair(Jet::SRC::PARTICLE, (*p)->barcode())));
        fjMCParticle2.set_user_index(idpart);
        particles.push_back(fjMCParticle2);
        idpart++;
      }

      std::unique_ptr<fastjet::JetDefinition> jetdef(new fastjet::JetDefinition(m_jetalgo, m_jetr, m_recomb_scheme, fastjet::Best));

      fastjet::ClusterSequence jetFinder(particles, *jetdef);
      std::vector<fastjet::PseudoJet> mcfastjets = jetFinder.inclusive_jets();

      mcTaggedJet = new Jetv1();

      for (auto &mcfastjet : mcfastjets)
      {
        bool isTaggedJet = false;
        std::vector<fastjet::PseudoJet> comps = mcfastjet.constituents();
        for (auto &comp : comps)
        {
          mcTaggedJet->insert_comp(fjMapMC[comp.user_index()].first, fjMapMC[comp.user_index()].second);
          if (comp.user_index() == 0)
          {
            mcTaggedJet->set_px(mcfastjet.px());
            mcTaggedJet->set_py(mcfastjet.py());
            mcTaggedJet->set_pz(mcfastjet.pz());
            mcTaggedJet->set_e(mcfastjet.e());
            mcTaggedJet->set_id(m_truth_jet_id);
            isTaggedJet = true;
          }
        }
        if (isTaggedJet)
        {
          break;
        }
        else
        {
          mcTaggedJet->clear_comp();
        }
      }
      m_truth_taggedJetMap->insert(mcTaggedJet);
      m_truth_jet_id++;
    }
  }
}


// Inspired by KFParticle_DST::createParticleNode
int ResonanceJetTagging::createJetNode(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);

  PHCompositeNode *lowerNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!lowerNode)
  {
    lowerNode = new PHCompositeNode("DST");
    topNode->addNode(lowerNode);
    std::cout << "DST node added" << std::endl;
  }

  std::string baseName;
  std::string jetNodeName;
  std::string jetNodeNameMC;

  if (m_jetcontainer_name.empty())
  {
    baseName = "taggedJet";
  }
  else
  {
    baseName = m_jetcontainer_name;
  }

  // Cant have forward slashes in DST or else you make a subdirectory on save!!!
  std::string undrscr = "_";
  std::string nothing = "";
  std::map<std::string, std::string> forbiddenStrings;
  forbiddenStrings["/"] = undrscr;
  forbiddenStrings["("] = undrscr;
  forbiddenStrings[")"] = nothing;
  forbiddenStrings["+"] = "plus";
  forbiddenStrings["-"] = "minus";
  forbiddenStrings["*"] = "star";
  for (auto const &[badString, goodString] : forbiddenStrings)
  {
    size_t pos;
    while ((pos = baseName.find(badString)) != std::string::npos)
    {
      baseName.replace(pos, 1, goodString);
    }
  }

  jetNodeName = baseName + "_Jet_Container";
  jetNodeNameMC = baseName + "_Truth_Jet_Container";

  m_taggedJetMap = new JetMapv1();
  PHIODataNode<PHObject> *jetNode = new PHIODataNode<PHObject>(m_taggedJetMap, jetNodeName.c_str(), "PHObject");
  lowerNode->addNode(jetNode);
  std::cout << jetNodeName << " node added" << std::endl;

  m_truth_taggedJetMap = new JetMapv1();
  PHIODataNode<PHObject> *jetNodeMC = new PHIODataNode<PHObject>(m_truth_taggedJetMap, jetNodeNameMC.c_str(), "PHObject");
  lowerNode->addNode(jetNodeMC);
  std::cout << jetNodeNameMC << " node added" << std::endl;

  return Fun4AllReturnCodes::EVENT_OK;
}
