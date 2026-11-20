/*
 * Find decay topologies in HepMC
 * Cameron Dean
 * 04/06/2021
 */

#include "DecayFinder.h"

#include "DecayFinderContainerBase.h"  // for DecayFinderContainerBase::Iter
#include "DecayFinderContainer_v1.h"   // for DecayFinderContainer_v1

#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4VtxPoint.h>

#include <phhepmc/PHHepMCGenEvent.h>
#include <phhepmc/PHHepMCGenEventMap.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>    // for PHIODataNode
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/getClass.h>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <HepMC/GenEvent.h>
#include <HepMC/GenVertex.h>  // for GenVertex::particle_iterator
#pragma GCC diagnostic pop

#include <HepMC/GenParticle.h>
#include <HepMC/IteratorRange.h>
#include <HepMC/SimpleVector.h>

#include <TDatabasePDG.h>

#include <algorithm>  // for find, for_each, sort, transform
#include <array>
#include <cctype>    // for toupper
#include <cmath>     // for abs, sqrt, log, atan, tan
#include <cstddef>   // for ptrdiff_t, size_t
#include <cstdlib>   // for exit
#include <iostream>  // for operator<<, endl, basic_ostream
#include <iterator>  // for end, begin
#include <map>       // for map, map<>::mapped_type, _Rb...
#include <memory>    // for allocator_traits<>::value_type
#include <ranges>
#include <string_view>

namespace
{
  // Resonances that are ignored when searching for the decay. Each DecayFinder takes a copy and removes the particles
  // that are part of its own decay chain, see findDecay, so that finders with different decay descriptors do not
  // affect each other
  // This is an array of PID's which decay faster than we can see in the detector
  // This means if we find them in the HepMC record we need to skip over them
  // If your decay descriptor contains one of these IDs then it will be removed from the list before starting the search
  constexpr int defaultResonantPIDs[] = {
      111, 113, 213, 333, 310, 311, 313, 323, 413, 423, 513, 523, 441, 443, 100443, 9000111, 9000211, 100111, 100211,
      10111, 10211, 9010111, 9010211, 10113, 10213, 20113, 20213, 9000113, 9000213, 100113, 100213, 9010113, 9010213,
      9020113, 9020213, 30113, 30213, 9030113, 9030213, 9040113, 9040213, 115, 215, 10115, 10215, 9000115, 9000215,
      9010115, 9010215, 117, 217, 9000117, 9000217, 9010117, 9010217, 119, 219, 221, 331, 9000221, 9010221, 100221,
      10221, 9020221, 100331, 9030221, 10331, 9040221, 9050221, 9060221, 9070221, 9080221, 223, 10223, 20223, 10333,
      20333, 1000223, 9000223, 9010223, 30223, 100333, 225, 9000225, 335, 9010225, 9020225, 10225, 9030225, 10335,
      9040225, 9050225, 9060225, 9070225, 9080225, 9090225, 227, 337, 229, 9000229, 9010229, 9000311, 9000321, 10311,
      10321, 100311, 100321, 9010311, 9010321, 9020311, 9020321, 10313, 10323, 20313, 20323, 100313, 100323, 9000313,
      9000323, 30313, 30323, 315, 325, 9000315, 9000325, 10315, 10325, 20315, 20325, 9010315, 9010325, 9020315,
      9020325, 317, 327, 9010317, 9010327, 319, 329, 9000319, 9000329, 10411, 10421, 10413, 10423, 20413, 20423, 415,
      425, 431, 10431, 433, 10433, 20433, 435, 10511, 10521, 10513, 10523, 20513, 20523, 515, 525, 10531, 533, 10533,
      20533, 535, 10541, 543, 10543, 20543, 545, 10441, 100441, 10443, 20443, 30443, 9000443, 9010443, 9020443, 445,
      100445, 551, 10551, 100551, 110551, 200551, 210551, 10553, 20553, 30553, 110553, 120553, 130553, 210553,
      220553, 9000553, 9010553, 555, 10555, 20555, 100555, 110555, 120555, 200555, 557, 100557, 2224, 2214, 2114,
      1114, 3212, 3224, 3214, 3114, 3324, 3314, 4222, 4212, 4112, 4224, 4214, 4114, 4232, 4132, 4322, 4312, 4324,
      4314, 4332, 4334, 4412, 4422, 4414, 4424, 4432, 4434, 4444};

  // Neutral particles that are their own antiparticle, their PDG ID does not change sign in a charge conjugate decay
  constexpr auto kParticlesWithNoCC = std::to_array<int>(
      {111, 113, 115, 130, 220, 221, 223, 225, 310, 330, 331, 333, 335, 440, 441, 443, 445, 551, 553, 555, 10111,
       10113, 10221, 10223, 10331, 10333, 10441, 10443, 10551, 10553, 20113, 20223, 20333, 20443, 20553, 100443,
       100553});

  // These tracks require a + or - after their name for TDatabasePDG
  constexpr std::array<std::string_view, 4> kSpecialTracks{"e", "mu", "pi", "K"};

  // Pieces of the decay descriptor
  constexpr std::string_view kDecayArrow = "->";
  constexpr char kChargeIndicator = '^';
  constexpr char kStartIntermediate = '{';
  constexpr char kEndIntermediate = '}';

  // A decay chain this long means the search is stuck in a loop
  constexpr std::size_t kMaxDecayChainSize = 100;

  std::vector<int> absoluteValues(const std::vector<int> &values)
  {
    std::vector<int> absolutes;
    absolutes.reserve(values.size());
    for (const int value : values)
    {
      absolutes.push_back(std::abs(value));
    }
    return absolutes;
  }
}  // namespace

DecayFinder::DecayFinder()
  : SubsysReco("DECAYFINDER")
{
}

DecayFinder::DecayFinder(const std::string &name)
  : SubsysReco(name)
{
}

int DecayFinder::Init(PHCompositeNode *topNode)
{
  if (Verbosity() >= VERBOSITY_SOME)
  {
    std::cout << "DecayFinder name: " << Name() << std::endl;
    std::cout << "Decay descriptor: " << m_decayDescriptor << std::endl;
  }

  int canSearchDecay = parseDecayDescriptor();

  if (m_save_dst)
  {
    createDecayNode(topNode);
  }

  return canSearchDecay;
}

int DecayFinder::process_event(PHCompositeNode *topNode)
{
  int decayFound = findDecay(topNode);

  if (decayFound > 0 && m_save_dst && Verbosity() >= VERBOSITY_MORE)
  {
    printNode(topNode);
  }

  if (m_triggerOnDecay && decayFound < 1)
  {
    if (Verbosity() >= VERBOSITY_MORE)
    {
      std::cout << "The decay, " << m_decayDescriptor << " was not found or reconstructable in this event, skipping"
                << std::endl;
    }
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int DecayFinder::End(PHCompositeNode * /*topNode*/)
{
  if (Verbosity() >= VERBOSITY_SOME)
  {
    printInfo();
  }

  return 0;
}

/*
 *
 *This function tries to intepret the string you wrote
 *Will be semi-independent of the rest of this class
 *(semi- because we need this to set the charges and strings)
 *
 */
int DecayFinder::parseDecayDescriptor()
{
  // Everything below is built from the decay descriptor alone, so a second call must not add to the first
  // m_allowPhotons is not reset, it is a user setting
  m_getChargeConjugate = false;
  m_hasPhotonDaughter = false;
  m_mother_ID = 0;
  m_intermediates_ID.clear();
  m_daughters_ID.clear();
  m_nTracksFromMother = 0;
  m_nTracksFromIntermediates.clear();
  m_motherDecayProducts.clear();

  std::string descriptor = simplifyDecayDescriptor();
  std::vector<int> daughtersCharge;

  const bool motherParsed = parseMother(descriptor);
  const bool intermediatesParsed = parseIntermediates(descriptor, daughtersCharge);
  const bool tracksParsed = parseMotherTracks(descriptor, daughtersCharge);

  if (motherParsed && intermediatesParsed && tracksParsed)
  {
    // Only with every particle known, a particle that was not found has no charge in daughtersCharge
    assignCharges(daughtersCharge);

    if (Verbosity() >= VERBOSITY_MORE)
    {
      std::cout << "Your decay descriptor can be parsed" << std::endl;
    }
    return 0;
  }

  if (Verbosity() >= VERBOSITY_SOME)
  {
    std::cout << "Your decay descriptor cannot be parsed, " << Name() << " will not be registered" << std::endl;
  }
  return Fun4AllReturnCodes::DONOTREGISTERSUBSYSTEM;
}

std::string DecayFinder::simplifyDecayDescriptor()
{
  std::string descriptor = m_decayDescriptor;

  // Remove all white space before we begin
  std::erase(descriptor, ' ');

  // Check for charge conjugate requirement
  std::string checkForCC;
  if (descriptor.size() >= 4)  // Anything shorter cannot be "[...]CC" and substr would throw
  {
    checkForCC = descriptor.substr(0, 1) + descriptor.substr(descriptor.size() - 3, 3);
  }
  std::for_each(checkForCC.begin(), checkForCC.end(), [](char &c)
                { c = static_cast<char>(std::toupper(static_cast<unsigned char>(c))); });

  // Remove the CC check if needed
  if (checkForCC == "[]CC")
  {
    descriptor = descriptor.substr(1, descriptor.size() - 4);
    m_getChargeConjugate = true;
  }

  return descriptor;
}

bool DecayFinder::parseMother(std::string &descriptor)
{
  const size_t motherEndPoint = descriptor.find(kDecayArrow);
  const std::string mother = descriptor.substr(0, motherEndPoint);
  descriptor.erase(0, motherEndPoint + kDecayArrow.length());

  if (!findParticle(mother))
  {
    return false;
  }

  m_mother_ID = std::abs(get_pdgcode(mother));
  return true;
}

bool DecayFinder::parseIntermediates(std::string &descriptor, std::vector<int> &daughtersCharge)
{
  bool canBeParsed = true;

  size_t pos;
  while ((pos = descriptor.find(kStartIntermediate)) != std::string::npos)
  {
    const size_t endPoint = descriptor.find(kEndIntermediate, pos);
    if (endPoint == std::string::npos)
    {
      // Without this the intermediate is never removed from the descriptor and it is parsed forever
      if (Verbosity() >= VERBOSITY_SOME)
      {
        std::cout << "The intermediate decay in your decay descriptor has no closing " << kEndIntermediate << std::endl;
      }
      return false;
    }
    std::string intermediateDecay = descriptor.substr(pos + 1, endPoint - (pos + 1));

    const std::string intermediate = intermediateDecay.substr(0, intermediateDecay.find(kDecayArrow));
    if (findParticle(intermediate))
    {
      m_intermediates_ID.push_back(std::abs(get_pdgcode(intermediate)));
    }
    else
    {
      canBeParsed = false;
    }

    // Now find the daughters associated to this intermediate
    int nDaughters = 0;
    intermediateDecay.erase(0, intermediateDecay.find(kDecayArrow) + kDecayArrow.length());
    size_t daughterLocator;
    while ((daughterLocator = intermediateDecay.find(kChargeIndicator)) != std::string::npos)
    {
      int daughterID = 0;
      if (parseDaughter(intermediateDecay, daughterLocator, daughtersCharge, daughterID))
      {
        if (daughterID == 22)
        {
          m_hasPhotonDaughter = m_allowPhotons = true;
        }
      }
      else
      {
        canBeParsed = false;
      }
      ++nDaughters;
    }
    descriptor.erase(pos, endPoint + 1 - pos);
    m_nTracksFromIntermediates.push_back(nDaughters);
    m_nTracksFromMother += 1;
  }

  return canBeParsed;
}

bool DecayFinder::parseMotherTracks(std::string &descriptor, std::vector<int> &daughtersCharge)
{
  bool canBeParsed = true;

  size_t daughterLocator;
  while ((daughterLocator = descriptor.find(kChargeIndicator)) != std::string::npos)
  {
    int daughterID = 0;
    if (!parseDaughter(descriptor, daughterLocator, daughtersCharge, daughterID))
    {
      canBeParsed = false;
    }
    m_nTracksFromMother += 1;
  }

  return canBeParsed;
}

bool DecayFinder::parseDaughter(std::string &descriptor, const size_t chargeLocation,
                                std::vector<int> &daughtersCharge, int &daughterID)
{
  std::string daughter = descriptor.substr(0, chargeLocation);
  const std::string daughterChargeString = descriptor.substr(chargeLocation + 1, 1);
  if (std::find(kSpecialTracks.begin(), kSpecialTracks.end(), daughter) != kSpecialTracks.end())
  {
    daughter += daughterChargeString;
  }
  descriptor.erase(0, chargeLocation + 2);

  if (!findParticle(daughter))
  {
    return false;
  }

  daughterID = std::abs(get_pdgcode(daughter));
  m_daughters_ID.push_back(daughterID);
  daughtersCharge.push_back(get_charge(daughter));
  return true;
}

void DecayFinder::assignCharges(const std::vector<int> &daughtersCharge)
{
  std::vector<int> intermediates_charge;

  unsigned int trackStart = 0;
  unsigned int trackEnd = 0;
  for (unsigned int i = 0; i < m_intermediates_ID.size(); ++i)
  {
    trackStart = trackEnd;
    trackEnd = m_nTracksFromIntermediates[i] + trackStart;

    int vtxCharge = 0;

    for (unsigned int j = trackStart; j < trackEnd; ++j)
    {
      vtxCharge += daughtersCharge[j];
    }

    intermediates_charge.push_back(vtxCharge);
  }

  int mother_charge = 0;
  for (unsigned int i = 0; i < m_daughters_ID.size(); ++i)
  {
    mother_charge += daughtersCharge[i];
  }

  m_mother_ID = mother_charge == 0 ? m_mother_ID : mother_charge * m_mother_ID;
  for (unsigned int i = 0; i < m_intermediates_ID.size(); ++i)
  {
    m_intermediates_ID[i] =
        intermediates_charge[i] == 0 ? m_intermediates_ID[i] : intermediates_charge[i] * m_intermediates_ID[i];
    m_motherDecayProducts.push_back(m_intermediates_ID[i]);
  }
  for (unsigned int i = 0; i < m_daughters_ID.size(); ++i)
  {
    m_daughters_ID[i] = daughtersCharge[i] == 0 ? m_daughters_ID[i] : daughtersCharge[i] * m_daughters_ID[i];
    if (i >= trackEnd)
    {
      m_motherDecayProducts.push_back(m_daughters_ID[i]);
    }
  }
}

/*
 * Main body for searching your decay
 * Slight issue if you limit the generator decay
 * as decays wont enter the HepMC record
 * need a switch to go to Geant4 record
 */
int DecayFinder::findDecay(PHCompositeNode *topNode)
{
  resetResonances();
  const std::vector<int> positive_motherDecayProducts = absoluteValues(m_motherDecayProducts);

  getTruthNodes(topNode);

  // correct mother PID to search for match
  const int mother_id_to_match = m_getChargeConjugate ? std::abs(m_mother_ID) : m_mother_ID;

  if (m_truthinfo && !m_geneventmap)  // This should use the truth info container if we have no HepMC record
  {
    if (Verbosity() >= VERBOSITY_SOME)
    {
      std::cout << "You have a valid G4TruthInfo but no PHHepMCGenEventMap container" << std::endl;
      std::cout << "You are probably using a particle gun" << std::endl;
    }

    return searchGeant4Mothers(topNode, positive_motherDecayProducts, mother_id_to_match);
  }

  return searchHepMCMothers(topNode, positive_motherDecayProducts, mother_id_to_match);
}

void DecayFinder::resetResonances()
{
  // The mother and intermediates are not resonances to skip over. The IDs here can carry the sign of the charge so
  // compare absolute values
  m_resonantPIDs.assign(std::begin(defaultResonantPIDs), std::end(defaultResonantPIDs));
  auto removeResonance = [this](const int pdgID)
  {
    m_resonantPIDs.erase(std::remove(m_resonantPIDs.begin(), m_resonantPIDs.end(), std::abs(pdgID)),
                         m_resonantPIDs.end());
  };
  for (const int i : m_intermediates_ID)
  {
    removeResonance(i);
  }

  removeResonance(m_mother_ID);
}

void DecayFinder::getTruthNodes(PHCompositeNode *topNode)
{
  m_truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  if (!m_truthinfo)
  {
    std::cout << "DecayFinder: Missing node G4TruthInfo" << std::endl;
  }

  m_geneventmap = findNode::getClass<PHHepMCGenEventMap>(topNode, "PHHepMCGenEventMap");
  if (!m_geneventmap)
  {
    std::cout << "DecayFinder: Missing node PHHepMCGenEventMap" << std::endl;
  }

  if (!m_truthinfo && !m_geneventmap)
  {
    std::cout << "You have neither the PHHepMCGenEventMap or G4TruthInfo nodes" << std::endl;
    std::cout << "DecayFinder will crash, exiting now!" << std::endl;
    std::exit(1);
  }
}

int DecayFinder::searchGeant4Mothers(PHCompositeNode *topNode, const std::vector<int> &positiveProducts,
                                     const int motherIDToMatch)
{
  int reconstructableDecayWasFound = 0;

  PHG4TruthInfoContainer::ConstRange range = m_truthinfo->GetParticleRange();
  for (PHG4TruthInfoContainer::ConstIterator iter = range.first; iter != range.second; ++iter)
  {
    PHG4Particle *g4particle = iter->second;
    const int this_pid = m_getChargeConjugate ? std::abs(g4particle->get_pid()) : g4particle->get_pid();
    if (this_pid != motherIDToMatch)
    {
      continue;
    }

    if (Verbosity() >= VERBOSITY_MAX)
    {
      std::cout << "parent->pdg_id(): " << g4particle->get_pid() << std::endl;
    }

    SearchState state;
    m_decayChain.clear();
    m_decayChain.emplace_back(std::make_pair(g4particle->get_primary_id(), g4particle->get_barcode()),
                              g4particle->get_pid());

    searchGeant4Children(g4particle->get_track_id(), positiveProducts, state.breakLoop,
                         state.hasPhoton, state.hasPi0, state.failedPT, state.failedETA, state.products);

    if (evaluateCandidate(topNode, state))
    {
      ++reconstructableDecayWasFound;
    }
  }

  return reconstructableDecayWasFound;
}

int DecayFinder::searchHepMCMothers(PHCompositeNode *topNode, const std::vector<int> &positiveProducts,
                                    const int motherIDToMatch)
{
  int reconstructableDecayWasFound = 0;

  for (const auto &riter : std::ranges::reverse_view(*m_geneventmap))
  {
    m_genevt = riter.second;
    if (!m_genevt)
    {
      std::cout << "DecayFinder: Missing node PHHepMCGenEvent" << std::endl;
      continue;
    }

    reconstructableDecayWasFound += searchHepMCEvent(topNode, positiveProducts, motherIDToMatch);
  }

  return reconstructableDecayWasFound;
}

int DecayFinder::searchHepMCEvent(PHCompositeNode *topNode, const std::vector<int> &positiveProducts,
                                  const int motherIDToMatch)
{
  int reconstructableDecayWasFound = 0;

  HepMC::GenEvent *theEvent = m_genevt->getEvent();

  for (HepMC::GenEvent::particle_const_iterator p = theEvent->particles_begin(); p != theEvent->particles_end(); ++p)
  {
    const int this_pid = m_getChargeConjugate ? std::abs((*p)->pdg_id()) : (*p)->pdg_id();
    if (this_pid != motherIDToMatch)
    {
      continue;
    }

    if (Verbosity() >= VERBOSITY_MAX)
    {
      std::cout << "parent->pdg_id(): " << (*p)->pdg_id() << std::endl;
    }

    SearchState state;
    m_decayChain.clear();
    m_decayChain.emplace_back(std::make_pair(m_genevt->get_embedding_id(), (*p)->barcode()), (*p)->pdg_id());

    // Make sure that the mother has a decay in our record
    if (!(*p)->end_vertex())  // Mother has no end vertex, decay volume was limited
    {
      searchGeant4Record((*p)->barcode(), (*p)->pdg_id(), positiveProducts, state.breakLoop, state.hasPhoton,
                         state.hasPi0, state.failedPT, state.failedETA, state.products);
    }
    else
    {
      searchHepMCRecord((*p), positiveProducts, state.breakLoop, state.hasPhoton, state.hasPi0, state.failedPT,
                        state.failedETA, state.products);
    }

    if (evaluateCandidate(topNode, state))
    {
      ++reconstructableDecayWasFound;
    }
  }

  return reconstructableDecayWasFound;
}

bool DecayFinder::evaluateCandidate(PHCompositeNode *topNode, const SearchState &state)
{
  if (state.breakLoop || !compareDecays(m_motherDecayProducts, state.products))
  {
    return false;
  }

  m_counter += 1;
  if (state.failedPT || state.failedETA)
  {
    if (!state.failedETA)
    {
      m_nCandFail_pT += 1;
    }
    else if (!state.failedPT)
    {
      m_nCandFail_eta += 1;
    }
    else
    {
      m_nCandFail_pT_and_eta += 1;
    }
    return false;
  }

  m_nCandReconstructable += 1;
  if (m_save_dst)
  {
    fillDecayNode(topNode, m_decayChain);
  }
  countPhotonsAndPi0s(state.hasPhoton, state.hasPi0);

  return true;
}

void DecayFinder::countPhotonsAndPi0s(const bool hasPhoton, const bool hasPi0)
{
  if (hasPhoton && !hasPi0)
  {
    m_nCandHas_Photon += 1;
  }
  else if (!hasPhoton && hasPi0)
  {
    m_nCandHas_Pi0 += 1;
  }
  else if (hasPhoton && hasPi0)
  {
    m_nCandHas_Photon_and_Pi0 += 1;
  }
  else
  {
    m_nCandHas_noPhoton_and_noPi0 += 1;
  }
}

/*
 * Function to search HepMC record
 * Can switch to Geant4 search if needed
 */
void DecayFinder::searchHepMCRecord(HepMC::GenParticle *particle, const std::vector<int> &decayProducts,
                                    bool &breakLoop, bool &hasPhoton, bool &hasPi0, bool &failedPT, bool &failedETA,
                                    std::vector<int> &actualDecayProducts)
{
  for (HepMC::GenVertex::particle_iterator children = particle->end_vertex()->particles_begin(HepMC::children);
       children != particle->end_vertex()->particles_end(HepMC::children); ++children)
  {
    const int pdgID = (*children)->pdg_id();
    if (Verbosity() >= VERBOSITY_MAX)
    {
      std::cout << "--children->pdg_id(): " << pdgID << std::endl;
    }

    if (pdgID == 22)
    {
      hasPhoton = true;
    }
    if (pdgID == 111)
    {
      hasPi0 = true;
    }
    if ((!m_allowPhotons && pdgID == 22) || (!m_allowPi0 && pdgID == 111))
    {
      breakLoop = true;
      break;
    }

    if (std::find(decayProducts.begin(), decayProducts.end(), std::abs(pdgID)) != decayProducts.end())
    {
      // This is one of the children we are looking for
      if (Verbosity() >= VERBOSITY_MAX)
      {
        std::cout << "This is a child you were looking for" << std::endl;
      }
      searchWantedHepMCChild(*children, breakLoop, hasPhoton, hasPi0, failedPT, failedETA, actualDecayProducts);
    }
    else if (isIgnoredPhotonOrPi0(pdgID))
    {
      continue;  // An allowed photon or pi0 that is not part of the decay, as in the Geant4 search
    }
    else if (isResonance(pdgID))  // Now check if it's part of the resonance list
    {
      if (Verbosity() >= VERBOSITY_MAX)
      {
        std::cout << "This is a resonance to investigate further" << std::endl;
      }
      searchHepMCResonance(*children, decayProducts, breakLoop, hasPhoton, hasPi0, failedPT, failedETA,
                           actualDecayProducts);
    }
    else
    {
      breakLoop = true;  // This particle is not in the decay descriptor, stop
      break;
    }
  }
}

void DecayFinder::searchWantedHepMCChild(HepMC::GenParticle *child, bool &breakLoop, bool &hasPhoton, bool &hasPi0,
                                         bool &failedPT, bool &failedETA, std::vector<int> &actualDecayProducts)
{
  // Check if this is an internediate decay that didnt decay in the generator
  const std::vector<int> positive_intermediates_ID = getPositiveIntermediateIDs();
  const auto it =
      std::find(positive_intermediates_ID.begin(), positive_intermediates_ID.end(), std::abs(child->pdg_id()));

  if (!child->end_vertex() && it != positive_intermediates_ID.end())
  {
    const auto index = static_cast<unsigned int>(std::distance(positive_intermediates_ID.begin(), it));
    searchUndecayedIntermediate(child, index, breakLoop, hasPhoton, hasPi0, failedPT, failedETA,
                                actualDecayProducts);
  }
  else if (checkIfCorrectHepMCParticle(child, failedPT, failedETA))
  {
    recordHepMCParticle(*child, actualDecayProducts);
  }
}

void DecayFinder::searchUndecayedIntermediate(HepMC::GenParticle *intermediate, const unsigned int index,
                                              bool &breakLoop, bool &hasPhoton, bool &hasPi0, bool &failedPT,
                                              bool &failedETA, std::vector<int> &actualDecayProducts)
{
  const std::vector<int> requiredIntermediateDecayProducts = getIntermediateDaughterIDs(index);
  std::vector<int> actualIntermediateDecayProducts;
  const std::size_t chainSize = m_decayChain.size();

  // The tracks of an intermediate that is not needed must not fail the candidate, so keep their cuts apart
  bool probeFailedPT = false;
  bool probeFailedETA = false;
  searchGeant4Record(intermediate->barcode(), intermediate->pdg_id(),
                     absoluteValues(requiredIntermediateDecayProducts), breakLoop, hasPhoton, hasPi0, probeFailedPT,
                     probeFailedETA, actualIntermediateDecayProducts);

  if (compareDecays(requiredIntermediateDecayProducts, actualIntermediateDecayProducts))
  {
    failedPT = failedPT || probeFailedPT;
    failedETA = failedETA || probeFailedETA;
    recordHepMCParticle(*intermediate, actualDecayProducts);
  }
  else
  {
    restoreDecayChain(chainSize);  // Remove the products of the intermediate that was not needed
  }
}

void DecayFinder::searchHepMCResonance(HepMC::GenParticle *resonance, const std::vector<int> &decayProducts,
                                       bool &breakLoop, bool &hasPhoton, bool &hasPi0, bool &failedPT,
                                       bool &failedETA, std::vector<int> &actualDecayProducts)
{
  if (!resonance->end_vertex())
  {
    searchGeant4Record(resonance->barcode(), resonance->pdg_id(), decayProducts, breakLoop, hasPhoton, hasPi0,
                       failedPT, failedETA, actualDecayProducts);
    return;
  }

  for (HepMC::GenVertex::particle_iterator grandchildren = resonance->end_vertex()->particles_begin(HepMC::children);
       grandchildren != resonance->end_vertex()->particles_end(HepMC::children); ++grandchildren)
  {
    if (checkIfCorrectHepMCParticle(*grandchildren, failedPT, failedETA))
    {
      recordHepMCParticle(**grandchildren, actualDecayProducts);
    }
  }
}

void DecayFinder::recordHepMCParticle(const HepMC::GenParticle &particle, std::vector<int> &actualDecayProducts)
{
  actualDecayProducts.push_back(particle.pdg_id());
  m_decayChain.emplace_back(std::make_pair(m_genevt->get_embedding_id(), particle.barcode()), particle.pdg_id());
}

void DecayFinder::restoreDecayChain(const std::size_t size)
{
  if (m_decayChain.size() > size)
  {
    m_decayChain.erase(m_decayChain.begin() + static_cast<std::ptrdiff_t>(size), m_decayChain.end());
  }
}

/*
 *Function to search Geant4 record
 * Cannot switch to HepMC search (doesn't make sense to)
 */
void DecayFinder::searchGeant4Record(int barcode, int pid, const std::vector<int> &decayProducts, bool &breakLoop,
                                     bool &hasPhoton, bool &hasPi0, bool &failedPT, bool &failedETA,
                                     std::vector<int> &actualDecayProducts)
{
  // A HepMC only event can still have a decay that is only complete in the Geant4 record, this candidate is incomplete
  if (!m_truthinfo)
  {
    if (Verbosity() >= VERBOSITY_SOME)
    {
      std::cout << "DecayFinder: The decay is not complete in the HepMC record and there is no G4TruthInfo node"
                << std::endl;
    }
    breakLoop = true;
    return;
  }

  // The HepMC record only gives the barcode and PDG ID, from here the children are found from the Geant4 parent link
  PHG4TruthInfoContainer::ConstRange range = m_truthinfo->GetParticleRange();
  for (PHG4TruthInfoContainer::ConstIterator iter = range.first; iter != range.second; ++iter)
  {
    const PHG4Particle *parent = iter->second;
    if (parent->get_barcode() == barcode && std::abs(parent->get_pid()) == std::abs(pid))
    {
      searchGeant4Children(parent->get_track_id(), decayProducts, breakLoop, hasPhoton, hasPi0, failedPT, failedETA,
                           actualDecayProducts);
      if (breakLoop)
      {
        break;
      }
    }
  }
}

void DecayFinder::searchGeant4Children(const int parentTrackID, const std::vector<int> &decayProducts,
                                       bool &breakLoop, bool &hasPhoton, bool &hasPi0, bool &failedPT,
                                       bool &failedETA, std::vector<int> &actualDecayProducts)
{
  PHG4TruthInfoContainer::ConstRange range = m_truthinfo->GetParticleRange();
  for (PHG4TruthInfoContainer::ConstIterator iter = range.first; iter != range.second; ++iter)
  {
    if (m_decayChain.size() == kMaxDecayChainSize)
    {
      breakLoop = true;  // Stuck in loop. Sympton not cause!
      break;
    }

    // The track ID identifies the parent, the barcode does not for particles from a particle gun
    PHG4Particle *g4particle = iter->second;
    if (g4particle->get_parent_id() != parentTrackID)
    {
      continue;
    }

    if (processGeant4Child(g4particle, decayProducts, breakLoop, hasPhoton, hasPi0, failedPT, failedETA,
                           actualDecayProducts))
    {
      break;
    }
  }
}

bool DecayFinder::processGeant4Child(PHG4Particle *child, const std::vector<int> &decayProducts, bool &breakLoop,
                                     bool &hasPhoton, bool &hasPi0, bool &failedPT, bool &failedETA,
                                     std::vector<int> &actualDecayProducts)
{
  const int particleID = child->get_pid();
  if (Verbosity() >= VERBOSITY_MAX)
  {
    std::cout << "--children->pdg_id(): " << particleID << std::endl;
  }
  if (particleID == 22)
  {
    hasPhoton = true;
  }
  if (particleID == 111)
  {
    hasPi0 = true;
  }

  if (isRejectedPhotonOrPi0(particleID))
  {
    breakLoop = true;
    return true;
  }

  if (std::find(decayProducts.begin(), decayProducts.end(), std::abs(particleID)) != decayProducts.end())
  {
    // This is one of the children we are looking for
    if (checkIfCorrectGeant4Particle(child, hasPhoton, hasPi0, failedPT, failedETA))
    {
      if (Verbosity() >= VERBOSITY_MAX)
      {
        std::cout << "This is a child you were looking for" << std::endl;
      }
      actualDecayProducts.push_back(particleID);
      const int embedding_id = m_geneventmap ? m_genevt->get_embedding_id() : child->get_primary_id();
      m_decayChain.emplace_back(std::make_pair(embedding_id, child->get_barcode()), particleID);
    }
  }
  else if (isIgnoredPhotonOrPi0(particleID))
  {
    return false;
  }
  else if (isResonance(particleID))  // Now check if it's part of the other resonance list
  {
    if (Verbosity() >= VERBOSITY_MAX)
    {
      std::cout << "This is a resonance to investigate further" << std::endl;
    }
    searchGeant4Children(child->get_track_id(), decayProducts, breakLoop, hasPhoton, hasPi0, failedPT,
                         failedETA, actualDecayProducts);
  }
  else
  {
    breakLoop = true;  // This particle is not in the decay descriptor, stop
    return true;
  }

  return false;
}

/*
 * Checks if this prticle really matches your request
 */
bool DecayFinder::checkIfCorrectHepMCParticle(HepMC::GenParticle *particle, bool &trackFailedPT, bool &trackFailedETA)
{
  // Check if it is an intermediate or a final track
  const std::vector<int> positive_intermediates_ID = getPositiveIntermediateIDs();
  const auto it =
      std::find(positive_intermediates_ID.begin(), positive_intermediates_ID.end(), std::abs(particle->pdg_id()));

  if (it == positive_intermediates_ID.end())
  {
    if (Verbosity() >= VERBOSITY_MAX)
    {
      std::cout << "This is a final state track" << std::endl;
    }
    checkHepMCTrack(*particle, trackFailedPT, trackFailedETA);
    return true;
  }

  // Which intermediate decay list to we need
  const auto index = static_cast<unsigned int>(std::distance(positive_intermediates_ID.begin(), it));
  const std::vector<int> requiredIntermediateDecayProducts = getIntermediateDaughterIDs(index);
  std::vector<int> actualIntermediateDecayProducts;
  const std::size_t chainSize = m_decayChain.size();

  // The tracks of an intermediate that is not needed must not fail the candidate, so keep their cuts apart
  bool probeFailedPT = false;
  bool probeFailedETA = false;
  const bool vetoed = searchHepMCIntermediate(*particle, absoluteValues(requiredIntermediateDecayProducts),
                                              probeFailedPT, probeFailedETA, actualIntermediateDecayProducts);

  const bool acceptParticle =
      !vetoed && compareDecays(requiredIntermediateDecayProducts, actualIntermediateDecayProducts);
  if (acceptParticle)
  {
    trackFailedPT = trackFailedPT || probeFailedPT;
    trackFailedETA = trackFailedETA || probeFailedETA;
  }
  else
  {
    restoreDecayChain(chainSize);  // Remove the products of the intermediate that was not needed
  }

  return acceptParticle;
}

bool DecayFinder::searchHepMCIntermediate(const HepMC::GenParticle &intermediate,
                                          const std::vector<int> &positiveRequiredProducts, bool &failedPT,
                                          bool &failedETA, std::vector<int> &actualIntermediateDecayProducts)
{
  for (HepMC::GenVertex::particle_iterator grandchildren = intermediate.end_vertex()->particles_begin(HepMC::children);
       grandchildren != intermediate.end_vertex()->particles_end(HepMC::children); ++grandchildren)
  {
    const int pdgID = (*grandchildren)->pdg_id();
    if (Verbosity() >= VERBOSITY_MAX)
    {
      std::cout << "----grandchildren->pdg_id(): " << pdgID << std::endl;
    }

    if (isResonance(pdgID))
    {
      if (searchIntermediateResonance(**grandchildren, positiveRequiredProducts, failedPT, failedETA,
                                      actualIntermediateDecayProducts))
      {
        return true;
      }
    }
    else if (isIgnoredPhotonOrPi0(pdgID))
    {
      continue;
    }
    else if (isRejectedPhotonOrPi0(pdgID))
    {
      return true;  // A photon or pi0 that is not allowed rules out this decay, whatever was found before it
    }
    else
    {
      addHepMCIntermediateProduct(**grandchildren, failedPT, failedETA, actualIntermediateDecayProducts);
    }
  }

  return false;
}

bool DecayFinder::searchIntermediateResonance(const HepMC::GenParticle &resonance,
                                              const std::vector<int> &positiveRequiredProducts, bool &failedPT,
                                              bool &failedETA, std::vector<int> &actualIntermediateDecayProducts)
{
  HepMC::GenVertex *endVertex = resonance.end_vertex();
  if (!endVertex)
  {
    // The generator did not decay it and so Geant4 did. A photon or pi0 is dealt with as any other, otherwise the
    // decay is in the Geant4 record
    const int pdgID = resonance.pdg_id();
    if (isIgnoredPhotonOrPi0(pdgID))
    {
      return false;
    }
    if (isRejectedPhotonOrPi0(pdgID))
    {
      return true;
    }

    bool vetoed = false;
    bool hasPhoton = false;
    bool hasPi0 = false;
    searchGeant4Record(resonance.barcode(), pdgID, positiveRequiredProducts, vetoed, hasPhoton, hasPi0, failedPT,
                       failedETA, actualIntermediateDecayProducts);
    return vetoed;
  }

  for (HepMC::GenVertex::particle_iterator greatgrandchildren = endVertex->particles_begin(HepMC::children);
       greatgrandchildren != endVertex->particles_end(HepMC::children); ++greatgrandchildren)
  {
    const int greatgrandchildID = (*greatgrandchildren)->pdg_id();
    if (Verbosity() >= VERBOSITY_MAX)
    {
      std::cout << "--------greatgrandchildren->pdg_id(): " << greatgrandchildID << std::endl;
    }

    if (isIgnoredPhotonOrPi0(greatgrandchildID))
    {
      continue;
    }
    if (isRejectedPhotonOrPi0(greatgrandchildID))
    {
      return true;
    }

    addHepMCIntermediateProduct(**greatgrandchildren, failedPT, failedETA, actualIntermediateDecayProducts);
  }

  return false;
}

void DecayFinder::addHepMCIntermediateProduct(const HepMC::GenParticle &particle, bool &failedPT, bool &failedETA,
                                              std::vector<int> &actualIntermediateDecayProducts)
{
  recordHepMCParticle(particle, actualIntermediateDecayProducts);
  checkHepMCTrack(particle, failedPT, failedETA);
}

void DecayFinder::checkHepMCTrack(const HepMC::GenParticle &particle, bool &failedPT, bool &failedETA)
{
  const HepMC::FourVector &momentum = particle.momentum();

  const HepMC::GenVertex *thisVtx = particle.production_vertex();
  double vtxPos[3] = {thisVtx->point3d().x(), thisVtx->point3d().y(), thisVtx->point3d().z()};

  applyTrackRequirements(momentum.perp(), momentum.eta(), momentum.py(), vtxPos, failedPT, failedETA);
}

void DecayFinder::applyTrackRequirements(const double pt, const double eta, const double py, double vertex[3],
                                         bool &failedPT, bool &failedETA)
{
  if (m_recalcualteEtaRange)
  {
    recalculateEta(py, vertex);
  }

  if (pt < m_pt_req)
  {
    failedPT = true;
  }
  if (!isInRange(m_eta_low_req, eta, m_eta_high_req))
  {
    failedETA = true;
  }

  if (Verbosity() >= VERBOSITY_MAX)
  {
    std::cout << "pT = " << pt << ", eta = " << eta << std::endl;
    std::cout << "The track " << passOrFail(pt >= m_pt_req) << " the pT requirement, ";
    std::cout << "the track " << passOrFail(isInRange(m_eta_low_req, eta, m_eta_high_req))
              << " the eta requirement." << std::endl;
  }
}

/*
 * Checks if this prticle really matches your request
 */
bool DecayFinder::checkIfCorrectGeant4Particle(PHG4Particle *particle, bool &hasPhoton, bool &hasPi0,
                                               bool &trackFailedPT, bool &trackFailedETA)
{
  // Check if it is an intermediate or a final track
  const std::vector<int> positive_intermediates_ID = getPositiveIntermediateIDs();
  const auto it =
      std::find(positive_intermediates_ID.begin(), positive_intermediates_ID.end(), std::abs(particle->get_pid()));

  if (it == positive_intermediates_ID.end())
  {
    if (Verbosity() >= VERBOSITY_MAX)
    {
      std::cout << "This is a final state track" << std::endl;
    }
    const double px = particle->get_px();
    const double py = particle->get_py();
    const double pz = particle->get_pz();

    const double pt = std::sqrt(std::pow(px, 2) + std::pow(py, 2));
    const double p = std::sqrt(std::pow(px, 2) + std::pow(py, 2) + std::pow(pz, 2));
    const double eta = 0.5 * std::log((p + pz) / (p - pz));

    const PHG4VtxPoint *thisVtx = m_truthinfo->GetVtx(particle->get_vtx_id());
    double vtxPos[3] = {thisVtx->get_x(), thisVtx->get_y(), thisVtx->get_z()};

    applyTrackRequirements(pt, eta, py, vtxPos, trackFailedPT, trackFailedETA);
    return true;
  }

  // Which intermediate decay list to we need
  const auto index = static_cast<unsigned int>(std::distance(positive_intermediates_ID.begin(), it));
  const std::vector<int> requiredIntermediateDecayProducts = getIntermediateDaughterIDs(index);
  std::vector<int> actualIntermediateDecayProducts;

  const std::size_t chainSize = m_decayChain.size();

  // The tracks of an intermediate that is not needed must not fail the candidate, so keep their cuts apart
  bool probeFailedPT = false;
  bool probeFailedETA = false;
  bool fakeBreak = false;
  searchGeant4Children(particle->get_track_id(), absoluteValues(requiredIntermediateDecayProducts), fakeBreak,
                       hasPhoton, hasPi0, probeFailedPT, probeFailedETA, actualIntermediateDecayProducts);

  // A break means a photon or pi0 that is not allowed, or a particle that is not in the decay, was found, so the
  // products that were found before it do not make this the intermediate decay that is needed
  const bool acceptParticle =
      !fakeBreak && compareDecays(requiredIntermediateDecayProducts, actualIntermediateDecayProducts);
  if (acceptParticle)
  {
    trackFailedPT = trackFailedPT || probeFailedPT;
    trackFailedETA = trackFailedETA || probeFailedETA;
  }
  else
  {
    restoreDecayChain(chainSize);  // Remove the products of the intermediate that was not needed
  }

  return acceptParticle;
}

/*
 * Check if the decay matches what we need
 */
bool DecayFinder::compareDecays(std::vector<int> required, std::vector<int> actual)
{
  bool accept = false;
  multiplyVectorByScalarAndSort(required, +1);
  multiplyVectorByScalarAndSort(actual, +1);

  if (Verbosity() >= VERBOSITY_MAX)
  {
    printComparison(required, actual);
  }

  if (required == actual)
  {
    accept = true;
  }

  if (m_getChargeConjugate && !accept)
  {
    multiplyVectorByScalarAndSort(required, -1);
    if (required == actual)
    {
      accept = true;
    }

    if (Verbosity() >= VERBOSITY_MAX)
    {
      std::cout << "Checking CC state" << std::endl;
      std::cout << (required == actual ? "*\n* These vectors match\n*\n" : "*\n* These vectors DONT match\n*\n")
                << std::endl;
    }
  }
  return accept;
}

void DecayFinder::printComparison(const std::vector<int> &required, const std::vector<int> &actual) const
{
  std::cout << "Printing required decay products: ";
  for (const int i : required)
  {
    std::cout << i << ", ";
  }
  std::cout << std::endl;
  std::cout << "Printing actual decay products: ";
  for (const int i : actual)
  {
    std::cout << i << ", ";
  }
  std::cout << std::endl;
  std::cout << (required == actual ? "*\n* These vectors match\n*\n" : "*\n* These vectors DONT match\n*\n")
            << std::endl;
}

/*
 * Everything below this is helper functions
 */
bool DecayFinder::findParticle(const std::string &particle)
{
  bool particleFound = true;
  if (!TDatabasePDG::Instance()->GetParticle(particle.c_str()))
  {
    if (Verbosity() >= VERBOSITY_SOME)
    {
      std::cout << "The particle, " << particle << " is not in the TDatabasePDG particle list" << std::endl;
    }
    particleFound = false;
  }

  return particleFound;
}

std::vector<int> DecayFinder::getPositiveIntermediateIDs() const
{
  return absoluteValues(m_intermediates_ID);
}

std::vector<int> DecayFinder::getIntermediateDaughterIDs(const unsigned int index) const
{
  // The tracks of the intermediates are stored one intermediate after the other
  unsigned int trackStart = 0;
  for (unsigned int i = 0; i < index; ++i)
  {
    trackStart += m_nTracksFromIntermediates[i];
  }
  const unsigned int trackStop = trackStart + m_nTracksFromIntermediates[index];

  std::vector<int> daughterIDs;
  for (unsigned int i = trackStart; i < trackStop; ++i)
  {
    daughterIDs.push_back(m_daughters_ID[i]);
  }
  return daughterIDs;
}

bool DecayFinder::isResonance(const int pdgID) const
{
  return std::find(std::begin(m_resonantPIDs), std::end(m_resonantPIDs), std::abs(pdgID)) !=
         std::end(m_resonantPIDs);
}

bool DecayFinder::isIgnoredPhotonOrPi0(const int pdgID) const
{
  return (m_allowPhotons && !m_hasPhotonDaughter && pdgID == 22) || (m_allowPi0 && pdgID == 111);
}

bool DecayFinder::isRejectedPhotonOrPi0(const int pdgID) const
{
  return (!m_allowPhotons && pdgID == 22) || (!m_hasPhotonDaughter && pdgID == 22) || (!m_allowPi0 && pdgID == 111);
}

void DecayFinder::multiplyVectorByScalarAndSort(std::vector<int> &v, int k)
{
  // https://slaystudy.com/c-multiply-vector-by-scalar/
  std::transform(v.begin(), v.end(), v.begin(), [k](const int &c)
                 {
                   if (std::find(kParticlesWithNoCC.begin(), kParticlesWithNoCC.end(), c) != kParticlesWithNoCC.end())
                   {
                     return c;
                   }

                   return c * k; });
  std::sort(v.begin(), v.end());
}

int DecayFinder::get_pdgcode(const std::string &name)
{
  if (findParticle(name))
  {
    return TDatabasePDG::Instance()->GetParticle(name.c_str())->PdgCode();
  }

  return 0;
}

int DecayFinder::get_charge(const std::string &name)
{
  if (findParticle(name))
  {
    return static_cast<int>(TDatabasePDG::Instance()->GetParticle(name.c_str())->Charge() / 3);
  }

  return -99;
}

bool DecayFinder::isInRange(float min, float value, float max)
{
  return min <= value && value <= max;
}

void DecayFinder::calculateEffectiveTPCradius(double vertex[3], double &effective_top_r,
                                              double &effective_bottom_r) const
{
  double r_sq = std::pow(m_tpc_r, 2);
  double x_sq = std::pow(vertex[0], 2);

  effective_top_r = std::sqrt(r_sq - x_sq) - vertex[1];
  effective_bottom_r = std::sqrt(r_sq - x_sq) + vertex[1];
}

void DecayFinder::recalculateEta(double py, double vertex[3])
{
  calculateEffectiveTPCradius(vertex, m_effective_top_tpc_r, m_effective_bottom_tpc_r);

  double z_north = m_tpc_z - vertex[2];  // Distance between north end of TPC and vertex z position
  double z_south = m_tpc_z + vertex[2];  // Same but for south of TPC

  double top_north_theta = std::atan(m_effective_top_tpc_r / z_north);
  double bottom_north_theta = std::atan(m_effective_bottom_tpc_r / z_north);
  double top_south_theta = std::atan(m_effective_top_tpc_r / z_south);
  double bottom_south_theta = std::atan(m_effective_bottom_tpc_r / z_south);

  double top_north_eta = -1 * std::log(std::tan(top_north_theta / 2));
  double bottom_north_eta = -1 * std::log(std::tan(bottom_north_theta / 2));
  double top_south_eta = std::log(std::tan(top_south_theta / 2));
  double bottom_south_eta = std::log(std::tan(bottom_south_theta / 2));

  m_eta_high_req = py < 0 ? bottom_north_eta : top_north_eta;
  m_eta_low_req = py < 0 ? bottom_south_eta : top_south_eta;
}

int DecayFinder::createDecayNode(PHCompositeNode *topNode)
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

  if (m_container_name.empty())
  {
    // baseName = "decay";
    baseName = Name();
  }
  else
  {
    baseName = m_container_name;
  }

  // Cant have forward slashes in DST or else you make a subdirectory on save!!!
  size_t pos;
  std::string undrscr = "_";
  std::string nothing;
  std::map<std::string, std::string> forbiddenStrings;
  forbiddenStrings["/"] = undrscr;
  forbiddenStrings["("] = undrscr;
  forbiddenStrings[")"] = nothing;
  forbiddenStrings["+"] = "plus";
  forbiddenStrings["-"] = "minus";
  forbiddenStrings["*"] = "star";
  forbiddenStrings[" "] = nothing;
  for (auto const &[badString, goodString] : forbiddenStrings)
  {
    while ((pos = baseName.find(badString)) != std::string::npos)
    {
      baseName.replace(pos, 1, goodString);
    }
  }

  m_nodeName = baseName + "_DecayMap";

  m_decayMap = new DecayFinderContainer_v1();
  PHIODataNode<PHObject> *decayNode = new PHIODataNode<PHObject>(m_decayMap, m_nodeName, "PHObject");
  lowerNode->addNode(decayNode);
  std::cout << m_nodeName << " node added" << std::endl;

  return Fun4AllReturnCodes::EVENT_OK;
}

void DecayFinder::fillDecayNode(PHCompositeNode *topNode, Decay &decay)
{
  m_decayMap = findNode::getClass<DecayFinderContainer_v1>(topNode, m_nodeName);
  m_decayMap->insert(decay);
}

void DecayFinder::printInfo()
{
  std::cout << "\n---------------DecayFinder information---------------" << std::endl;
  std::cout << "Module name: " << Name() << std::endl;
  std::cout << "Decay descriptor: " << m_decayDescriptor << std::endl;
  std::cout << "Number of generated decays: " << m_counter << std::endl;
  std::cout << "  Number of decays that failed pT requirement: " << m_nCandFail_pT << std::endl;
  std::cout << "  Number of decays that failed eta requirement: " << m_nCandFail_eta << std::endl;
  std::cout << "  Number of decays that failed pT and eta requirements: " << m_nCandFail_pT_and_eta << std::endl;
  std::cout << "Number of decays that could be reconstructed: " << m_nCandReconstructable << std::endl;
  std::cout << "  Number of reconstructable mothers with associated photon: " << m_nCandHas_Photon << std::endl;
  std::cout << "  Number of reconstructable mothers with associated Pi0: " << m_nCandHas_Pi0 << std::endl;
  std::cout << "  Number of reconstructable mothers with associated photon and Pi0: " << m_nCandHas_Photon_and_Pi0
            << std::endl;
  std::cout << "  Number of reconstructable mothers with no associated photons or Pi0s: "
            << m_nCandHas_noPhoton_and_noPi0 << std::endl;
  std::cout << "-----------------------------------------------------\n";
  std::cout << std::endl;
}

void DecayFinder::printNode(PHCompositeNode *topNode)
{
  std::cout << "----------------";
  std::cout << " DecayFinderNode: " << m_nodeName << " information ";
  std::cout << "----------------" << std::endl;
  DecayFinderContainer_v1 *map = findNode::getClass<DecayFinderContainer_v1>(topNode, m_nodeName);
  for (auto &iter : *map)
  {
    Decay decay = iter.second;
    for (auto &i : decay)
    {
      std::cout << "Particle embedding ID: " << i.first.first << ", particle barcode: " << i.first.second
                << ", particle PDG ID: " << i.second << std::endl;
    }
  }
  std::cout << "--------------------------------------------------------------------------------------------------"
            << std::endl;
}

std::string DecayFinder::passOrFail(bool condition)
{
  return condition ? "passed" : "failed";
}
