#include "HFTrackEfficiency.h"

#include <decayfinder/DecayFinderContainerBase.h>

#include <trackbase/TrkrClusterHitAssoc.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrHitTruthAssoc.h>

#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrackMap_v2.h>
#include <trackbase_historic/TrackSeed.h>
#include <trackbase_historic/TrackSeedContainer.h>

#include <g4eval/SvtxClusterEval.h>
#include <g4eval/SvtxEvalStack.h>
#include <g4eval/SvtxTrackEval.h>
#include <g4eval/SvtxTruthEval.h>

#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4VtxPoint.h>

#include <phhepmc/PHHepMCGenEvent.h>
#include <phhepmc/PHHepMCGenEventMap.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>

#include <TDatabasePDG.h>
#include <TFile.h>
#include <TParticlePDG.h>
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
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <set>
#include <type_traits>

namespace
{
  // Speed of light in cm/ps, to turn a flight distance into a proper decay time
  constexpr double kCmPerPs = 0.0299792458;

  // PDG IDs (absolute value) of the particles that leave a track
  constexpr std::array<int, 5> kTrackablePDGs{11, 13, 211, 321, 2212};

  // Largest difference in a momentum component [GeV] that still counts as the same particle when there is no barcode
  // match, and the relative size of that difference for fast particles
  constexpr double kMomentumTolerance = 5e-3;

  // Number of vertices to follow back through a record, a guard against a malformed record
  constexpr int kMaxAncestryDepth = 200;

  // Vertices in the HepMC record closer than this are the same vertex, in the HepMC length unit
  constexpr double kVertexEpsilon = 1e-6;

  // Size of a vertex shift [cm] that is worth telling the user about
  constexpr double kVertexShiftWarning = 1e-3;

  // Rapidity of a particle moving along the beam with E <= |pz|, where it is not defined
  constexpr float kInvalidRapidity = -999.;

  constexpr float kNaN = std::numeric_limits<float>::quiet_NaN();

  bool isTrackable(const int pdgID)
  {
    return std::ranges::find(kTrackablePDGs, std::abs(pdgID)) != kTrackablePDGs.end();
  }

  std::string getParticleName(const int pdgID)
  {
    const TParticlePDG *particle = TDatabasePDG::Instance()->GetParticle(pdgID);
    return particle ? particle->GetName() : "unknown";
  }

  float getParticleMass(const int pdgID)
  {
    const TParticlePDG *particle = TDatabasePDG::Instance()->GetParticle(pdgID);
    return particle ? static_cast<float>(particle->Mass()) : kNaN;
  }

  // Add a branch with the leaf type of the variable
  template <typename T>
  void addBranch(TTree &tree, const std::string &name, T *address)
  {
    char leafType = 'F';
    if constexpr (std::is_same_v<T, bool>)
    {
      leafType = 'O';
    }
    else if constexpr (std::is_same_v<T, int>)
    {
      leafType = 'I';
    }
    else
    {
      static_assert(std::is_same_v<T, float>, "No leaf type for this branch type");
    }
    tree.Branch(name.c_str(), address, (name + "/" + leafType).c_str());
  }

  // Index of the detector in the cluster branches, -1 if the detector is not one of them
  int getDetectorIndex(const TrkrDefs::cluskey key)
  {
    switch (TrkrDefs::getTrkrId(key))
    {
    case TrkrDefs::mvtxId:
      return 0;
    case TrkrDefs::inttId:
      return 1;
    case TrkrDefs::tpcId:
      return 2;
    case TrkrDefs::micromegasId:
      return 3;
    default:
      return -1;
    }
  }

  // Does the reco track have the momentum of the true one to within the given fraction in all three components
  bool isMomentumMatch(const SvtxTrack &track, const CLHEP::HepLorentzVector &trueLV, const float fraction)
  {
    const auto delta_px = static_cast<float>((track.get_px() - trueLV.px()) / trueLV.px());
    const auto delta_py = static_cast<float>((track.get_py() - trueLV.py()) / trueLV.py());
    const auto delta_pz = static_cast<float>((track.get_pz() - trueLV.pz()) / trueLV.pz());
    return std::abs(delta_px) <= fraction && std::abs(delta_py) <= fraction && std::abs(delta_pz) <= fraction;
  }

  // Does the Geant4 particle have the PDG ID and barcode of an entry in the decay
  template <typename Entry>
  bool isDecayEntry(const PHG4Particle *particle, const Entry &entry)
  {
    return particle && particle->get_pid() == entry.second && particle->get_barcode() == entry.first.second;
  }

  // The parent of a Geant4 particle, nullptr if it has none
  PHG4Particle *getParent(PHG4TruthInfoContainer *truthInfo, const PHG4Particle *particle)
  {
    return particle->get_parent_id() != 0 ? truthInfo->GetParticle(particle->get_parent_id()) : nullptr;
  }

  // The child of an ancestor that a descendant comes from, which is the descendant itself for a direct child.
  // nullptr if the descendant does not come from the ancestor
  const PHG4Particle *getChildOf(PHG4TruthInfoContainer *truthInfo, const PHG4Particle *ancestor,
                                 const PHG4Particle *descendant)
  {
    const PHG4Particle *child = descendant;
    for (int depth = 0; child && depth < kMaxAncestryDepth; ++depth)
    {
      const PHG4Particle *parent = getParent(truthInfo, child);
      if (parent && parent->get_track_id() == ancestor->get_track_id())
      {
        return child;
      }
      child = parent;
    }
    return nullptr;
  }

  // Follow the first incoming particle of every vertex back to the start of the HepMC record
  HepMC::GenVertex *findRootVertex(const HepMC::GenParticle &particle)
  {
    HepMC::GenVertex *vtx = particle.production_vertex();
    for (int depth = 0; vtx && depth < kMaxAncestryDepth; ++depth)
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
  CLHEP::Hep3Vector getLabFrameVertex(const HepMC::GenVertex &vtx, const HepMC::GenEvent &evt,
                                      const CLHEP::HepLorentzRotation &evtGen2Lab,
                                      const HepMC::FourVector &collisionVertex)
  {
    const double lengthFactor = HepMC::Units::conversion_factor(evt.length_unit(), HepMC::Units::CM);
    const CLHEP::HepLorentzVector genVertex(vtx.position().x(), vtx.position().y(), vtx.position().z(),
                                            vtx.position().t());
    const CLHEP::HepLorentzVector labVertex = evtGen2Lab(genVertex);
    return CLHEP::Hep3Vector(labVertex.x() * lengthFactor + collisionVertex.x(),
                             labVertex.y() * lengthFactor + collisionVertex.y(),
                             labVertex.z() * lengthFactor + collisionVertex.z());
  }

  // Four momentum of a HepMC particle in GeV, whatever the unit of the event is, as HepMCNodeReader does
  HepMC::FourVector getMomentumInGeV(const HepMC::GenParticle &particle, const HepMC::GenEvent &evt)
  {
    const double momentumFactor = HepMC::Units::conversion_factor(evt.momentum_unit(), HepMC::Units::GEV);
    const HepMC::FourVector &momentum = particle.momentum();
    return HepMC::FourVector(momentum.px() * momentumFactor, momentum.py() * momentumFactor,
                             momentum.pz() * momentumFactor, momentum.e() * momentumFactor);
  }

  // Are two HepMC vertices at the same place
  bool isSameVertex(const HepMC::GenVertex &first, const HepMC::GenVertex &second)
  {
    return std::abs(first.position().x() - second.position().x()) < kVertexEpsilon &&
           std::abs(first.position().y() - second.position().y()) < kVertexEpsilon &&
           std::abs(first.position().z() - second.position().z()) < kVertexEpsilon;
  }
}  // namespace

struct HFTrackEfficiency::DecayTruth
{
  HepMC::GenEvent *event{nullptr};  // HepMC event of the mother, nullptr if the decay is only in the Geant4 record

  CLHEP::HepLorentzVector daughterSumLV;  // Sum of the true daughters that leave a track
  CLHEP::HepLorentzVector recoSumLV;      // Sum of the matched reco tracks, with the true masses
  std::set<SvtxTrack *> selectedTracks;   // The matched reco tracks

  // True vertices, all in cm and in the lab frame that Geant4 (and so the reconstruction) uses
  bool haveVertices{false};
  CLHEP::Hep3Vector primaryVtx;
  CLHEP::Hep3Vector productionVtx;
  CLHEP::Hep3Vector decayVtx;

  // The HepMC record is in the event generator frame. HepMCNodeReader rotates it to the lab frame (crossing angle) and
  // shifts it by the collision vertex, possibly plus a vertex offset/smearing of its own that is not visible from
  // here. The latter is recovered by comparing one HepMC daughter vertex to the Geant4 vertex of the same particle.
  CLHEP::HepLorentzRotation evtGen2Lab;
  double motherLabP{0.};  // Momentum of the mother in the lab frame [GeV], for the proper time
  CLHEP::Hep3Vector vertexOffset;
  bool vertexOffsetKnown{false};
};

struct HFTrackEfficiency::DaughterTruth
{
  PHG4Particle *g4{nullptr};  // Geant4 particle of the daughter, nullptr if there is none
  bool haveKinematics{false};
  CLHEP::HepLorentzVector lv;     // Four vector of the true daughter, generator frame if it came from HepMC
  CLHEP::HepLorentzVector labLV;  // Four vector of the true daughter in the lab frame, the frame of the reco tracks
};

//____________________________________________________________________________..
HFTrackEfficiency::HFTrackEfficiency(const std::string &name)
  : SubsysReco(name)
{
}

HFTrackEfficiency::~HFTrackEfficiency() = default;

//____________________________________________________________________________..
int HFTrackEfficiency::Init(PHCompositeNode *topNode)
{
  m_truthRecoMatchPercent /= 100.F;

  if (m_write_track_map)
  {
    createOutputTrackMap(topNode);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int HFTrackEfficiency::process_event(PHCompositeNode *topNode)
{
  int returnCode = getNodes(topNode);
  if (returnCode != Fun4AllReturnCodes::EVENT_OK)
  {
    return returnCode;
  }

  setupEvaluation(topNode);

  returnCode = setupDecay();
  if (returnCode != Fun4AllReturnCodes::EVENT_OK)
  {
    return returnCode;
  }

  const bool reconstructableDecay = processDecays(topNode);

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

//____________________________________________________________________________..
int HFTrackEfficiency::End(PHCompositeNode * /*topNode*/)
{
  PrintEff();

  if (m_write_nTuple && m_counter_allDecays != 0 && m_outfile)
  {
    m_outfile->Write();
    m_outfile->Close();
    m_outfile.reset();  // Closing the file deletes the tree
    m_tree = nullptr;
  }

  return Fun4AllReturnCodes::EVENT_OK;
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
    std::cout << "Efficiency: " << std::setw(4)
              << static_cast<float>(m_counter_acceptedDecays) * 100 / static_cast<float>(m_counter_allDecays) << "%"
              << std::endl;
  }
  if (m_counter_unmatchedTruth != 0)
  {
    std::cout << "Warning: " << m_counter_unmatchedTruth
              << " daughters had no Geant4 particle, these use the momentum matching to find a reco track" << std::endl;
  }
  std::cout << "----------------------------------------------------------------\n"
            << std::endl;
}

void HFTrackEfficiency::createOutputTrackMap(PHCompositeNode *topNode)
{
  PHNodeIterator nodeIter(topNode);

  auto *dstNode = dynamic_cast<PHCompositeNode *>(nodeIter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    dstNode = new PHCompositeNode("DST");
    topNode->addNode(dstNode);
    std::cout << "DST node added" << std::endl;
  }

  m_outputNodeName = m_output_track_map_node_name + "_SvtxTrackMap";
  m_output_trackMap = new SvtxTrackMap_v2();
  auto *outputTrackNode = new PHIODataNode<PHObject>(m_output_trackMap, m_outputNodeName, "PHObject");
  dstNode->addNode(outputTrackNode);
  std::cout << m_outputNodeName << " node added" << std::endl;
}

int HFTrackEfficiency::getNodes(PHCompositeNode *topNode)
{
  const std::string df_node_name = m_df_module_name + "_DecayMap";
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

  getClusterTruthNodes(topNode);

  return Fun4AllReturnCodes::EVENT_OK;
}

void HFTrackEfficiency::getClusterTruthNodes(PHCompositeNode *topNode)
{
  // SvtxClusterEval dereferences these without checking them, so they have to be there before it is asked for the
  // truth particles of a cluster. Same order as the detector index
  static constexpr std::array<const char *, s_nDetectors> kG4HitNodeNames{"G4HIT_MVTX", "G4HIT_INTT", "G4HIT_TPC",
                                                                          "G4HIT_MICROMEGAS"};

  m_haveClusterTruthAssoc = findNode::getClass<TrkrClusterHitAssoc>(topNode, "TRKR_CLUSTERHITASSOC") &&
                            findNode::getClass<TrkrHitTruthAssoc>(topNode, "TRKR_HITTRUTHASSOC");
  for (unsigned int iDet = 0; iDet < s_nDetectors; ++iDet)
  {
    m_haveDetectorG4Hits[iDet] = findNode::getClass<PHG4HitContainer>(topNode, kG4HitNodeNames[iDet]) != nullptr;
  }

  if (!m_haveClusterTruthAssoc && Verbosity() > 0)
  {
    std::cout << __FILE__ << ": Missing TRKR_CLUSTERHITASSOC or TRKR_HITTRUTHASSOC, the true cluster counts are -1"
              << std::endl;
  }
}

void HFTrackEfficiency::setupEvaluation(PHCompositeNode *topNode)
{
  if (!m_svtx_evalstack)
  {
    m_svtx_evalstack = std::make_unique<SvtxEvalStack>(topNode);
    m_trackeval = m_svtx_evalstack->get_track_eval();
    m_trackeval->set_track_nodename(m_input_track_map_node_name);
  }
  m_svtx_evalstack->next_event(topNode);
}

int HFTrackEfficiency::setupDecay()
{
  if (!m_decay_descriptor.empty() || m_decayMap->empty())
  {
    return Fun4AllReturnCodes::EVENT_OK;
  }

  getDecayDescriptor();
  getNDaughters();

  if (m_nDaughters > s_maxTracks)
  {
    std::cout << __FILE__ << ": The decay has " << m_nDaughters << " trackable daughters but at most " << s_maxTracks
              << " are supported" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  if (m_nDaughters != 2 && m_nDaughters != 3)
  {
    std::cout << __FILE__ << ": Warning, " << m_nDaughters
              << " trackable daughters. Decays found in the Geant4 record only support 2 or 3" << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

bool HFTrackEfficiency::processDecays(PHCompositeNode *topNode)
{
  bool reconstructableDecay = false;
  for (const auto &iter : *m_decayMap)
  {
    const Decay decay = iter.second;
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

  return reconstructableDecay;
}

bool HFTrackEfficiency::findTracks(PHCompositeNode *topNode, const Decay &decay)
{
  resetBranches();

  DecayTruth truth;
  if (m_geneventmap)
  {
    fillHepMCMother(decay, truth);
  }

  unsigned int nTrackable = 0;
  for (unsigned int i = 1; i < decay.size(); ++i)
  {
    if (!isTrackable(decay[i].second))
    {
      continue;
    }

    if (nTrackable >= s_maxTracks)
    {
      std::cout << __FILE__ << ": Decay has more than " << s_maxTracks << " trackable daughters, increase s_maxTracks"
                << std::endl;
      return false;
    }

    processDaughter(decay, i, nTrackable, truth);
    ++nTrackable;
  }

  m_true_mother_mass = static_cast<float>(truth.daughterSumLV.m());
  fillVertexBranches(decay, truth);

  return storeSelectedTracks(topNode, truth);
}

void HFTrackEfficiency::processDaughter(const Decay &decay, const unsigned int decayIndex,
                                        const unsigned int trackIndex, DecayTruth &truth)
{
  DaughterTruth daughter;
  findDaughterTruth(decay, decayIndex, trackIndex, truth, daughter);

  if (!daughter.g4)
  {
    ++m_counter_unmatchedTruth;
    if (Verbosity() > 0)
    {
      std::cout << __FILE__ << ": No Geant4 particle found for daughter " << trackIndex + 1 << " (PDG "
                << decay[decayIndex].second << ")" << std::endl;
    }
  }

  if (!daughter.haveKinematics)
  {
    return;  // Nothing to compare a reco track to
  }

  fillTrueDaughterBranches(trackIndex, daughter);

  SvtxTrack *track = findRecoTrack(trackIndex, daughter);
  if (!track)
  {
    return;
  }

  if (Verbosity() >= VERBOSITY_MORE)
  {
    track->identify();
  }

  truth.selectedTracks.insert(track);
  fillRecoTrackBranches(track, daughter.g4, trackIndex);

  CLHEP::HepLorentzVector daughterRecoLV;
  daughterRecoLV.setVectM(CLHEP::Hep3Vector(track->get_px(), track->get_py(), track->get_pz()),
                          getParticleMass(static_cast<int>(m_true_track_PID[trackIndex])));
  truth.recoSumLV += daughterRecoLV;
}

void HFTrackEfficiency::fillHepMCMother(const Decay &decay, DecayTruth &truth)
{
  m_genevt = m_geneventmap->get(decay[0].first.first);
  if (!m_genevt)
  {
    return;  // The decay is not in the HepMC record, the truth comes from Geant4
  }

  truth.event = m_genevt->getEvent();
  const HepMC::GenParticle *mother = truth.event ? truth.event->barcode_to_particle(decay[0].first.second) : nullptr;
  if (!mother)
  {
    truth.event = nullptr;
    return;
  }

  if (Verbosity() >= VERBOSITY_MORE)
  {
    mother->print();
  }

  truth.evtGen2Lab = m_genevt->get_LorentzRotation_EvtGen2Lab();
  fillHepMCMotherKinematics(*mother, truth);
  fillHepMCVertices(*mother, truth);
}

void HFTrackEfficiency::fillHepMCMotherKinematics(const HepMC::GenParticle &mother, DecayTruth &truth)
{
  // These are the generator frame values, the transformation to the lab frame (crossing angle rotation and boost) is
  // not applied. Compare with the Geant4 truth at the 1e-3 level
  const HepMC::FourVector momentum = getMomentumInGeV(mother, *truth.event);
  m_true_mother_pT = static_cast<float>(momentum.perp());
  // Must have an old HepMC build, no mag function
  m_true_mother_p = static_cast<float>(std::sqrt(momentum.px() * momentum.px() + momentum.py() * momentum.py() +
                                                 momentum.pz() * momentum.pz()));
  m_true_mother_eta = static_cast<float>(momentum.eta());
  m_true_mother_phi = static_cast<float>(momentum.phi());
  if (momentum.e() > std::abs(momentum.pz()))
  {
    m_true_mother_rapidity =
        static_cast<float>(0.5 * std::log((momentum.e() + momentum.pz()) / (momentum.e() - momentum.pz())));
  }
  else
  {
    m_true_mother_rapidity = kInvalidRapidity;
  }

  // The flight distance is in the lab frame so the proper time needs the lab frame momentum, the boost changes it
  const CLHEP::HepLorentzVector labMomentum =
      truth.evtGen2Lab(CLHEP::HepLorentzVector(momentum.px(), momentum.py(), momentum.pz(), momentum.e()));
  truth.motherLabP = labMomentum.vect().mag();
}

void HFTrackEfficiency::fillHepMCVertices(const HepMC::GenParticle &mother, DecayTruth &truth)
{
  const HepMC::GenVertex *motherProductionVtx = mother.production_vertex();
  const HepMC::GenVertex *motherDecayVtx = mother.end_vertex();
  const HepMC::GenVertex *rootVtx = findRootVertex(mother);

  if (!motherProductionVtx || !motherDecayVtx || !rootVtx)
  {
    if (Verbosity() > 0)
    {
      std::cout << __FILE__
                << ": The mother in the HepMC record has no production or decay vertex, no true vertices will be filled"
                << std::endl;
    }
    return;
  }

  const HepMC::FourVector &collisionVtx = m_genevt->get_collision_vertex();
  truth.primaryVtx = getLabFrameVertex(*rootVtx, *truth.event, truth.evtGen2Lab, collisionVtx);
  truth.productionVtx = getLabFrameVertex(*motherProductionVtx, *truth.event, truth.evtGen2Lab, collisionVtx);
  truth.decayVtx = getLabFrameVertex(*motherDecayVtx, *truth.event, truth.evtGen2Lab, collisionVtx);
  truth.haveVertices = true;

  // A prompt mother is produced where the record starts. Compared in the HepMC length unit, no reliance on the
  // record starting at the origin
  m_is_primary = isSameVertex(*motherProductionVtx, *rootVtx);
}

void HFTrackEfficiency::findDaughterTruth(const Decay &decay, const unsigned int decayIndex,
                                          const unsigned int trackIndex, DecayTruth &truth, DaughterTruth &daughter)
{
  const HepMC::GenParticle *daughterHepMC = nullptr;
  if (truth.event && decay[decayIndex].first.second > -1)
  {
    daughterHepMC = truth.event->barcode_to_particle(decay[decayIndex].first.second);
  }

  if (daughterHepMC)
  {
    fillHepMCDaughter(decay, decayIndex, trackIndex, *daughterHepMC, truth, daughter);
  }
  else if (m_truthInfo)
  {
    findGeant4Daughter(decay, decayIndex, trackIndex, truth, daughter);
  }
}

void HFTrackEfficiency::fillHepMCDaughter(const Decay &decay, const unsigned int decayIndex,
                                          const unsigned int trackIndex, const HepMC::GenParticle &daughterHepMC,
                                          DecayTruth &truth, DaughterTruth &daughter)
{
  if (Verbosity() >= VERBOSITY_MORE)
  {
    daughterHepMC.print();
  }

  const HepMC::FourVector momentum = getMomentumInGeV(daughterHepMC, *truth.event);
  daughter.lv.setVectM(CLHEP::Hep3Vector(momentum.px(), momentum.py(), momentum.pz()),
                       getParticleMass(decay[decayIndex].second));
  truth.daughterSumLV += daughter.lv;
  daughter.haveKinematics = true;

  m_true_track_PID[trackIndex] = static_cast<float>(daughterHepMC.pdg_id());

  // The Geant4 and reco momenta are in the lab frame so rotate the HepMC one before comparing with them
  daughter.labLV =
      truth.evtGen2Lab(CLHEP::HepLorentzVector(momentum.px(), momentum.py(), momentum.pz(), momentum.e()));

  // We need the G4 particle, not the HepMC one, to use the truth/reco matching
  daughter.g4 = findGeant4Primary(decay[decayIndex].first.first, decay[decayIndex].first.second,
                                  decay[decayIndex].second, daughter.labLV.vect());

  if (daughter.g4 && truth.haveVertices && !truth.vertexOffsetKnown)
  {
    updateVertexOffset(daughterHepMC, *daughter.g4, truth);
  }
}

void HFTrackEfficiency::updateVertexOffset(const HepMC::GenParticle &daughterHepMC, const PHG4Particle &daughterG4,
                                           DecayTruth &truth)
{
  const HepMC::GenVertex *daughterProductionVtx = daughterHepMC.production_vertex();
  const PHG4VtxPoint *daughterG4Vtx = m_truthInfo->GetVtx(daughterG4.get_vtx_id());
  if (!daughterProductionVtx || !daughterG4Vtx)
  {
    return;
  }

  const CLHEP::Hep3Vector expectedVtx =
      getLabFrameVertex(*daughterProductionVtx, *truth.event, truth.evtGen2Lab, m_genevt->get_collision_vertex());
  truth.vertexOffset.set(daughterG4Vtx->get_x() - expectedVtx.x(), daughterG4Vtx->get_y() - expectedVtx.y(),
                         daughterG4Vtx->get_z() - expectedVtx.z());
  truth.vertexOffsetKnown = true;

  if (Verbosity() > 0 && truth.vertexOffset.mag() > kVertexShiftWarning)
  {
    std::cout << __FILE__ << ": The Geant4 vertices are shifted by " << truth.vertexOffset
              << " cm relative to the HepMC record in the lab frame, correcting the true vertices" << std::endl;
  }
}

void HFTrackEfficiency::findGeant4Daughter(const Decay &decay, const unsigned int decayIndex,
                                           const unsigned int trackIndex, DecayTruth &truth, DaughterTruth &daughter)
{
  PHG4TruthInfoContainer::ConstRange range = m_truthInfo->GetParticleRange();
  for (PHG4TruthInfoContainer::ConstIterator iter = range.first; iter != range.second; ++iter)
  {
    PHG4Particle *candidate = iter->second;
    PHG4Particle *motherG4 = findGeant4Mother(decay, decayIndex, candidate);
    if (motherG4)
    {
      fillGeant4Truth(decay, motherG4, candidate, trackIndex, truth, daughter);
      return;
    }
  }
}

PHG4Particle *HFTrackEfficiency::findGeant4Mother(const Decay &decay, const unsigned int decayIndex,
                                                  const PHG4Particle *candidate) const
{
  if (!isDecayEntry(candidate, decay[decayIndex]))
  {
    return nullptr;
  }

  PHG4Particle *parent = getParent(m_truthInfo, candidate);
  if (!parent)
  {
    return nullptr;
  }

  // Two daughters, the parent is the mother
  if (m_nDaughters == 2)
  {
    return isDecayEntry(parent, decay[0]) ? parent : nullptr;
  }

  // Three daughters, direct decay [mother, d1, d2, d3], every daughter comes from the mother
  if (m_nDaughters == 3 && decay.size() == 4)
  {
    return isDecayEntry(parent, decay[0]) ? parent : nullptr;
  }

  // Three daughters, the decay is [mother, d1, d2, intermediate, d3]. d3 comes from the mother and d1 and d2 from the
  // intermediate, which comes from the mother
  if (m_nDaughters == 3 && decay.size() > 4)
  {
    if (decayIndex == 4)
    {
      return isDecayEntry(parent, decay[0]) ? parent : nullptr;
    }

    if (!isDecayEntry(parent, decay[3]))
    {
      return nullptr;
    }
    PHG4Particle *grandparent = getParent(m_truthInfo, parent);
    return isDecayEntry(grandparent, decay[0]) ? grandparent : nullptr;
  }

  return nullptr;
}

void HFTrackEfficiency::fillGeant4Truth(const Decay &decay, PHG4Particle *motherG4, PHG4Particle *daughterG4,
                                        const unsigned int trackIndex, DecayTruth &truth, DaughterTruth &daughter)
{
  if (Verbosity() >= VERBOSITY_MORE)
  {
    daughterG4->identify();
  }

  daughter.g4 = daughterG4;
  daughter.haveKinematics = true;

  m_is_primary = m_truthInfo->is_sPHENIX_primary(motherG4);

  const CLHEP::Hep3Vector mother3Vector(motherG4->get_px(), motherG4->get_py(), motherG4->get_pz());
  CLHEP::HepLorentzVector motherTrueLV;
  motherTrueLV.setVectM(mother3Vector, getParticleMass(decay[0].second));
  m_true_mother_pT = static_cast<float>(motherTrueLV.perp());
  m_true_mother_p = static_cast<float>(mother3Vector.mag());
  truth.motherLabP = mother3Vector.mag();
  m_true_mother_eta = static_cast<float>(motherTrueLV.pseudoRapidity());
  m_true_mother_phi = static_cast<float>(motherTrueLV.phi());
  m_true_mother_rapidity = static_cast<float>(motherTrueLV.rapidity());

  daughter.lv.setVectM(CLHEP::Hep3Vector(daughterG4->get_px(), daughterG4->get_py(), daughterG4->get_pz()),
                       getParticleMass(daughterG4->get_pid()));
  daughter.labLV = daughter.lv;  // The Geant4 truth is already in the lab frame
  truth.daughterSumLV += daughter.lv;

  // A HepMC mother already provided the vertices
  if (!truth.haveVertices)
  {
    truth.haveVertices =
        getGeant4Vertices(motherG4, daughterG4, truth.primaryVtx, truth.productionVtx, truth.decayVtx);
  }

  m_true_track_PID[trackIndex] = static_cast<float>(daughterG4->get_pid());
}

PHG4Particle *HFTrackEfficiency::findGeant4Primary(const int embeddingID, const int barcode, const int pdgID,
                                                   const CLHEP::Hep3Vector &labMomentum) const
{
  if (!m_truthInfo)
  {
    return nullptr;
  }

  // HepMCNodeReader copies the HepMC barcode to the Geant4 primary and registers it with the embedding ID of the HepMC
  // event so that is the unambiguous match. The barcode is only unique within one HepMC event, hence the embedding ID.
  // If that fails fall back on the closest momentum (in the lab frame, as the Geant4 particles are)
  PHG4Particle *closest = nullptr;
  double closestDelta = std::numeric_limits<double>::max();
  const double tolerance = std::max(kMomentumTolerance, kMomentumTolerance * labMomentum.mag());

  PHG4TruthInfoContainer::ConstRange range = m_truthInfo->GetPrimaryParticleRange();
  for (PHG4TruthInfoContainer::ConstIterator iter = range.first; iter != range.second; ++iter)
  {
    PHG4Particle *candidate = iter->second;
    // Only the primaries of the same HepMC event can match, for the barcode and for the momentum
    if (candidate->get_pid() != pdgID || m_truthInfo->isEmbeded(candidate->get_track_id()) != embeddingID)
    {
      continue;
    }

    if (candidate->get_barcode() == barcode)
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

bool HFTrackEfficiency::getGeant4Vertices(const PHG4Particle *motherG4, const PHG4Particle *daughterG4,
                                          CLHEP::Hep3Vector &primary, CLHEP::Hep3Vector &production,
                                          CLHEP::Hep3Vector &decay) const
{
  // Follow the parents of the mother back to the particle that started the chain
  const PHG4Particle *root = motherG4;
  for (int depth = 0; depth < kMaxAncestryDepth; ++depth)
  {
    const PHG4Particle *parent = getParent(m_truthInfo, root);
    if (!parent)
    {
      break;
    }
    root = parent;
  }

  // The mother decays where its own child is produced. For a daughter of an intermediate that is the production
  // vertex of the intermediate, not the vertex of the daughter, which is further along
  const PHG4Particle *child = getChildOf(m_truthInfo, motherG4, daughterG4);
  if (!child)
  {
    return false;  // The daughter does not descend from the mother
  }

  const PHG4VtxPoint *primaryVtx = m_truthInfo->GetVtx(root->get_vtx_id());
  const PHG4VtxPoint *productionVtx = m_truthInfo->GetVtx(motherG4->get_vtx_id());
  const PHG4VtxPoint *decayVtx = m_truthInfo->GetVtx(child->get_vtx_id());
  if (!primaryVtx || !productionVtx || !decayVtx)
  {
    return false;
  }

  primary.set(primaryVtx->get_x(), primaryVtx->get_y(), primaryVtx->get_z());
  production.set(productionVtx->get_x(), productionVtx->get_y(), productionVtx->get_z());
  decay.set(decayVtx->get_x(), decayVtx->get_y(), decayVtx->get_z());

  return true;
}

void HFTrackEfficiency::fillTrueDaughterBranches(const unsigned int trackIndex, const DaughterTruth &daughter)
{
  m_true_track_pT[trackIndex] = static_cast<float>(daughter.lv.perp());
  m_true_track_eta[trackIndex] = static_cast<float>(daughter.lv.pseudoRapidity());
  m_true_track_rapidity[trackIndex] = static_cast<float>(daughter.lv.rapidity());
  m_true_track_phi[trackIndex] = static_cast<float>(daughter.lv.phi());
  m_min_true_track_pT = std::min(m_true_track_pT[trackIndex], m_min_true_track_pT);
  m_max_true_track_pT = std::max(m_true_track_pT[trackIndex], m_max_true_track_pT);
}

SvtxTrack *HFTrackEfficiency::findRecoTrack(const unsigned int trackIndex, const DaughterTruth &daughter)
{
  if (m_trackeval && daughter.g4)
  {
    SvtxTrack *track = m_trackeval->best_track_from(daughter.g4);
    m_used_truth_reco_map[trackIndex] = (track != nullptr);
    return track;
  }

  for (const auto &iter : *m_input_trackMap)
  {
    SvtxTrack *track = iter.second;
    if (isMomentumMatch(*track, daughter.labLV, m_truthRecoMatchPercent))
    {
      return track;
    }
  }

  return nullptr;
}

void HFTrackEfficiency::fillRecoTrackBranches(SvtxTrack *track, PHG4Particle *truthParticle,
                                              const unsigned int trackIndex)
{
  m_reco_track_exists[trackIndex] = true;
  m_reco_track_pT[trackIndex] = track->get_pt();
  m_reco_track_eta[trackIndex] = track->get_eta();
  m_reco_track_phi[trackIndex] = track->get_phi();
  m_reco_track_chi2nDoF[trackIndex] = track->get_chisq() / static_cast<float>(track->get_ndf());

  countClusters(track, truthParticle, trackIndex);

  m_min_reco_track_pT = std::min(m_reco_track_pT[trackIndex], m_min_reco_track_pT);
  m_max_reco_track_pT = std::max(m_reco_track_pT[trackIndex], m_max_reco_track_pT);
}

std::set<TrkrDefs::cluskey> HFTrackEfficiency::collectClusterKeys(SvtxTrack *track) const
{
  // Same cluster keys that SvtxTrackEval uses, but also valid for tracks that only store the seed IDs
  std::set<TrkrDefs::cluskey> clusterKeys;
  auto collectKeys = [&clusterKeys](const TrackSeed *seed)
  {
    if (seed)
    {
      clusterKeys.insert(seed->begin_cluster_keys(), seed->end_cluster_keys());
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

  return clusterKeys;
}

void HFTrackEfficiency::countClusters(SvtxTrack *track, PHG4Particle *truthParticle, const unsigned int trackIndex)
{
  const std::set<TrkrDefs::cluskey> clusterKeys = collectClusterKeys(track);

  SvtxClusterEval *clusterEval = m_trackeval->get_cluster_eval();
  SvtxTruthEval *truthEval = m_trackeval->get_truth_eval();

  // The true cluster count is -1, rather than 0, for any detector whose cluster truth cannot be looked up
  m_reco_track_nclusters[trackIndex].fill(0);
  for (unsigned int iDet = 0; iDet < s_nDetectors; ++iDet)
  {
    const bool canMatchTruth = truthParticle && m_haveClusterTruthAssoc && m_haveDetectorG4Hits[iDet];
    m_reco_track_ntrue_clusters[trackIndex][iDet] = canMatchTruth ? 0 : -1;
  }

  for (const TrkrDefs::cluskey key : clusterKeys)
  {
    const int iDet = getDetectorIndex(key);
    if (iDet < 0)
    {
      continue;
    }

    ++m_reco_track_nclusters[trackIndex][iDet];

    if (m_reco_track_ntrue_clusters[trackIndex][iDet] < 0)
    {
      continue;
    }

    // A cluster is true if the truth daughter left energy in it, counted once however many particles share the cluster
    const auto contributors = clusterEval->all_truth_particles(key);
    const bool isTrue = std::ranges::any_of(contributors, [truthEval, truthParticle](PHG4Particle *contributor)
                                            { return truthEval->are_same_particle(contributor, truthParticle); });
    if (isTrue)
    {
      ++m_reco_track_ntrue_clusters[trackIndex][iDet];
    }
  }
}

void HFTrackEfficiency::fillVertexBranches(const Decay &decay, DecayTruth &truth)
{
  if (!truth.haveVertices)
  {
    return;
  }

  if (truth.vertexOffsetKnown)
  {
    truth.primaryVtx += truth.vertexOffset;
    truth.productionVtx += truth.vertexOffset;
    truth.decayVtx += truth.vertexOffset;
  }

  m_primary_vtx_x = static_cast<float>(truth.primaryVtx.x());
  m_primary_vtx_y = static_cast<float>(truth.primaryVtx.y());
  m_primary_vtx_z = static_cast<float>(truth.primaryVtx.z());
  m_production_vtx_x = static_cast<float>(truth.productionVtx.x());
  m_production_vtx_y = static_cast<float>(truth.productionVtx.y());
  m_production_vtx_z = static_cast<float>(truth.productionVtx.z());
  m_decay_vtx_x = static_cast<float>(truth.decayVtx.x());
  m_decay_vtx_y = static_cast<float>(truth.decayVtx.y());
  m_decay_vtx_z = static_cast<float>(truth.decayVtx.z());

  // Flight of the mother from where it was produced to where it decays. Not measured from the primary vertex,
  // that would include the flight of any parent (e.g. the B for a non-prompt D0)
  m_true_decay_length = static_cast<float>((truth.decayVtx - truth.productionVtx).mag());

  // Proper time of the mother, t = L m / (p c), with the flight and the momentum both in the lab frame
  if (truth.motherLabP > 0)
  {
    m_true_decay_time = static_cast<float>(m_true_decay_length * getParticleMass(decay[0].second) /
                                           (truth.motherLabP * kCmPerPs));
  }
}

bool HFTrackEfficiency::storeSelectedTracks(PHCompositeNode *topNode, const DecayTruth &truth)
{
  if (truth.selectedTracks.size() != m_nDaughters)
  {
    return false;
  }

  m_reco_mother_mass = static_cast<float>(truth.recoSumLV.m());
  m_reco_mother_pT = static_cast<float>(truth.recoSumLV.perp());

  if (m_write_track_map)
  {
    m_output_trackMap = findNode::getClass<SvtxTrackMap>(topNode, m_outputNodeName);
    if (m_output_trackMap)
    {
      for (SvtxTrack *track : truth.selectedTracks)
      {
        m_output_trackMap->insertWithKey(track, track->get_id());
      }
    }
  }

  return true;
}

void HFTrackEfficiency::initializeBranches()
{
  m_outfile = std::make_unique<TFile>(m_outfile_name.c_str(), "RECREATE");
  m_tree = new TTree("HFTrackEfficiency", "HFTrackEfficiency");
  m_tree->OptimizeBaskets();
  m_tree->SetAutoSave(-5e6);  // Save the output file every 5MB

  addMotherBranches();
  for (unsigned int iTrack = 0; iTrack < m_nDaughters; ++iTrack)
  {
    addDaughterBranches(iTrack);
  }
  addVertexBranches();
}

void HFTrackEfficiency::addMotherBranches()
{
  addBranch(*m_tree, "all_tracks_reconstructed", &m_all_tracks_reconstructed);
  addBranch(*m_tree, "is_primary", &m_is_primary);
  addBranch(*m_tree, "true_mother_mass", &m_true_mother_mass);
  addBranch(*m_tree, "reco_mother_mass", &m_reco_mother_mass);
  addBranch(*m_tree, "true_mother_pT", &m_true_mother_pT);
  addBranch(*m_tree, "reco_mother_pT", &m_reco_mother_pT);
  addBranch(*m_tree, "true_mother_p", &m_true_mother_p);
  addBranch(*m_tree, "true_mother_eta", &m_true_mother_eta);
  addBranch(*m_tree, "true_mother_rapidity", &m_true_mother_rapidity);
  addBranch(*m_tree, "true_mother_phi", &m_true_mother_phi);
  addBranch(*m_tree, "min_true_track_pT", &m_min_true_track_pT);
  addBranch(*m_tree, "min_reco_track_pT", &m_min_reco_track_pT);
  addBranch(*m_tree, "max_true_track_pT", &m_max_true_track_pT);
  addBranch(*m_tree, "max_reco_track_pT", &m_max_reco_track_pT);
}

void HFTrackEfficiency::addDaughterBranches(const unsigned int trackIndex)
{
  const std::string track = "track_" + std::to_string(trackIndex + 1);
  const auto reco = [&track](const std::string &what)
  { return "reco_" + track + "_" + what; };
  const auto truth = [&track](const std::string &what)
  { return "true_" + track + "_" + what; };

  addBranch(*m_tree, reco("exists"), &m_reco_track_exists[trackIndex]);
  addBranch(*m_tree, reco("used_truth_reco_map"), &m_used_truth_reco_map[trackIndex]);
  addBranch(*m_tree, truth("pT"), &m_true_track_pT[trackIndex]);
  addBranch(*m_tree, reco("pT"), &m_reco_track_pT[trackIndex]);
  addBranch(*m_tree, truth("eta"), &m_true_track_eta[trackIndex]);
  addBranch(*m_tree, reco("eta"), &m_reco_track_eta[trackIndex]);
  addBranch(*m_tree, truth("rapidity"), &m_true_track_rapidity[trackIndex]);
  addBranch(*m_tree, truth("phi"), &m_true_track_phi[trackIndex]);
  addBranch(*m_tree, reco("phi"), &m_reco_track_phi[trackIndex]);
  addBranch(*m_tree, truth("PID"), &m_true_track_PID[trackIndex]);
  addBranch(*m_tree, reco("chi2nDoF"), &m_reco_track_chi2nDoF[trackIndex]);

  const std::array<std::string, s_nDetectors> detectorNames = {"mvtx", "intt", "tpc", "tpot"};
  for (unsigned int iDet = 0; iDet < s_nDetectors; ++iDet)
  {
    addBranch(*m_tree, reco("nclusters_" + detectorNames[iDet]), &m_reco_track_nclusters[trackIndex][iDet]);
    addBranch(*m_tree, reco("ntrue_clusters_" + detectorNames[iDet]), &m_reco_track_ntrue_clusters[trackIndex][iDet]);
  }
}

void HFTrackEfficiency::addVertexBranches()
{
  addBranch(*m_tree, "true_primary_vertex_x", &m_primary_vtx_x);
  addBranch(*m_tree, "true_primary_vertex_y", &m_primary_vtx_y);
  addBranch(*m_tree, "true_primary_vertex_z", &m_primary_vtx_z);
  addBranch(*m_tree, "true_production_vertex_x", &m_production_vtx_x);
  addBranch(*m_tree, "true_production_vertex_y", &m_production_vtx_y);
  addBranch(*m_tree, "true_production_vertex_z", &m_production_vtx_z);
  addBranch(*m_tree, "true_decay_vertex_x", &m_decay_vtx_x);
  addBranch(*m_tree, "true_decay_vertex_y", &m_decay_vtx_y);
  addBranch(*m_tree, "true_decay_vertex_z", &m_decay_vtx_z);
  addBranch(*m_tree, "true_decay_length", &m_true_decay_length);
  addBranch(*m_tree, "true_decay_time", &m_true_decay_time);
}

void HFTrackEfficiency::resetBranches()
{
  m_all_tracks_reconstructed = false;
  m_is_primary = false;
  m_true_mother_mass = kNaN;
  m_reco_mother_mass = kNaN;
  m_true_mother_pT = kNaN;
  m_reco_mother_pT = kNaN;
  m_true_mother_p = kNaN;
  m_true_mother_eta = kNaN;
  m_true_mother_rapidity = kNaN;
  m_true_mother_phi = kNaN;
  m_min_true_track_pT = std::numeric_limits<float>::max();
  m_min_reco_track_pT = std::numeric_limits<float>::max();
  m_max_true_track_pT = std::numeric_limits<float>::lowest();
  m_max_reco_track_pT = std::numeric_limits<float>::lowest();

  m_reco_track_exists.fill(false);
  m_used_truth_reco_map.fill(false);
  m_true_track_pT.fill(kNaN);
  m_reco_track_pT.fill(kNaN);
  m_true_track_eta.fill(kNaN);
  m_reco_track_eta.fill(kNaN);
  m_true_track_rapidity.fill(kNaN);
  m_true_track_phi.fill(kNaN);
  m_reco_track_phi.fill(kNaN);
  m_true_track_PID.fill(kNaN);
  m_reco_track_chi2nDoF.fill(kNaN);
  for (unsigned int iTrack = 0; iTrack < s_maxTracks; ++iTrack)
  {
    m_reco_track_nclusters[iTrack].fill(-1);
    m_reco_track_ntrue_clusters[iTrack].fill(-1);
  }

  m_primary_vtx_x = kNaN;
  m_primary_vtx_y = kNaN;
  m_primary_vtx_z = kNaN;
  m_production_vtx_x = kNaN;
  m_production_vtx_y = kNaN;
  m_production_vtx_z = kNaN;
  m_decay_vtx_x = kNaN;
  m_decay_vtx_y = kNaN;
  m_decay_vtx_z = kNaN;
  m_true_decay_length = kNaN;
  m_true_decay_time = kNaN;
}

void HFTrackEfficiency::getDecayDescriptor()
{
  const Decay decay = m_decayMap->begin()->second;
  m_decay_descriptor = getParticleName(decay[0].second);
  m_decay_descriptor += " ->";
  for (unsigned int i = 1; i < decay.size(); ++i)
  {
    m_decay_descriptor += ' ';
    m_decay_descriptor += getParticleName(decay[i].second);
  }
}

void HFTrackEfficiency::getNDaughters()
{
  const Decay decay = m_decayMap->begin()->second;
  m_nDaughters = 0;
  for (unsigned int i = 1; i < decay.size(); ++i)
  {
    if (isTrackable(decay[i].second))
    {
      ++m_nDaughters;
    }
  }
}
