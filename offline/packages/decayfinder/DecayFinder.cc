/*
 * Find decay topologies in HepMC
 * Cameron Dean
 * 04/06/2021
 */

//This is an array of PID's which decay faster than we can see in the detector
//This means if we find them in the HepMC record we need to skip over them
//If your decay descriptor contains one of these IDs then it will be removed from the list before starting the search

#include "DecayFinder.h"

#include "DecayFinderContainerBase.h"    // for DecayFinderContainerBase::Iter
#include "DecayFinderContainer_v1.h"     // for DecayFinderContainer_v1

#include <g4main/PHG4TruthInfoContainer.h>

#include <phhepmc/PHHepMCGenEvent.h>
#include <phhepmc/PHHepMCGenEventMap.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>          // for PHIODataNode
#include <phool/PHNodeIterator.h>        // for PHNodeIterator
#include <phool/getClass.h>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <HepMC/GenEvent.h>
#include <HepMC/GenVertex.h>             // for GenVertex::particle_iterator
#pragma GCC diagnostic pop

#include <HepMC/GenParticle.h>
#include <HepMC/IteratorRange.h>
#include <HepMC/SimpleVector.h>

#include <TDatabasePDG.h>

#include <algorithm>                     // for find, for_each, sort, transform
#include <ctype.h>                       // for toupper
#include <cstdlib>                      // for abs, size_t
#include <iostream>                      // for operator<<, endl, basic_ostream
#include <iterator>                      // for end, begin
#include <map>                           // for map, map<>::mapped_type, _Rb...
#include <memory>                        // for allocator_traits<>::value_type

int listOfResonantPIDs[] = {111, 113, 213, 333, 310, 311, 313, 323, 413, 423, 513, 523, 441, 443, 100443, 9000111, 9000211, 100111, 100211, 10111,
                            10211, 9010111, 9010211, 10113, 10213, 20113, 20213, 9000113, 9000213, 100113, 100213, 9010113, 9010213, 9020113, 9020213,
                            30113, 30213, 9030113, 9030213, 9040113, 9040213, 115, 215, 10115, 10215, 9000115, 9000215, 9010115, 9010215, 117, 217,
                            9000117, 9000217, 9010117, 9010217, 119, 219, 221, 331, 9000221, 9010221, 100221, 10221, 9020221, 100331, 9030221, 10331,
                            9040221, 9050221, 9060221, 9070221, 9080221, 223, 10223, 20223, 10333, 20333, 1000223, 9000223, 9010223, 30223, 100333, 225,
                            9000225, 335, 9010225, 9020225, 10225, 9030225, 10335, 9040225, 9050225, 9060225, 9070225, 9080225, 9090225, 227, 337, 229,
                            9000229, 9010229, 9000311, 9000321, 10311, 10321, 100311, 100321, 9010311, 9010321, 9020311, 9020321, 10313, 10323, 20313,
                            20323, 100313, 100323, 9000313, 9000323, 30313, 30323, 315, 325, 9000315, 9000325, 10315, 10325, 20315, 20325, 9010315,
                            9010325, 9020315, 9020325, 317, 327, 9010317, 9010327, 319, 329, 9000319, 9000329, 10411, 10421, 10413, 10423, 20413, 20423,
                            415, 425, 431, 10431, 433, 10433, 20433, 435, 10511, 10521, 10513, 10523, 20513, 20523, 515, 525, 10531, 533, 10533, 20533, 535,
                            10541, 543, 10543, 20543, 545, 10441, 100441, 10443, 20443, 30443, 9000443, 9010443, 9020443, 445, 100445, 551, 10551, 100551,
                            110551, 200551, 210551, 10553, 20553, 30553, 110553, 120553, 130553, 210553, 220553, 9000553, 9010553, 555, 10555, 20555, 100555,
                            110555, 120555, 200555, 557, 100557, 2224, 2214, 2114, 1114, 3212, 3224, 3214, 3114, 3324, 3314, 4222, 4212, 4112, 4224, 4214,
                            4114, 4232, 4132, 4322, 4312, 4324, 4314, 4332, 4334, 4412, 4422, 4414, 4424, 4432, 4434, 4444};

DecayFinder::DecayFinder()
  : SubsysReco("DECAYFINDER")
  , m_save_dst(false)
{
}

DecayFinder::DecayFinder(const std::string& name)
  : SubsysReco(name)
  , m_save_dst(false)
{
}

int DecayFinder::Init(PHCompositeNode* topNode)
{
  if (Verbosity() >= VERBOSITY_SOME)
  {
    std::cout << "DecayFinder name: " << Name() << std::endl;
    std::cout << "Decay descriptor: " << m_decayDescriptor << std::endl;
  }

  int canSearchDecay = parseDecayDescriptor();

  if (m_save_dst) createDecayNode(topNode);

  return canSearchDecay;
}

int DecayFinder::process_event(PHCompositeNode* topNode)
{
  bool decayFound = findDecay(topNode);

  if (decayFound && m_save_dst && Verbosity() >= VERBOSITY_MORE) printNode(topNode);

  if (m_triggerOnDecay && !decayFound)
  {
    if (Verbosity() >= VERBOSITY_MORE) std::cout << "The decay, " << m_decayDescriptor << " was not found or reconstructable in this event, skipping" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  else
  {
    return Fun4AllReturnCodes::EVENT_OK;
  }
}

int DecayFinder::End(PHCompositeNode* /*topNode*/)
{
  if (Verbosity() >= VERBOSITY_SOME) printInfo();

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
  bool ddCanBeParsed = true;

  size_t daughterLocator;

  std::string mother;
  std::string intermediate;
  std::string daughter;

  int mother_charge = 0;
  std::vector<int> intermediates_charge;
  std::vector<int> daughters_charge;

  std::string decayArrow = "->";
  std::string chargeIndicator = "^";
  std::string startIntermediate = "{";
  std::string endIntermediate = "}";

  //These tracks require a + or - after their name for TDatabasePDG
  std::string specialTracks[] = {"e", "mu", "pi", "K"};

  std::string manipulateDecayDescriptor = m_decayDescriptor;

  //Remove all white space before we begin
  size_t pos;
  while ((pos = manipulateDecayDescriptor.find(" ")) != std::string::npos) manipulateDecayDescriptor.replace(pos, 1, "");

  //Check for charge conjugate requirement
  std::string checkForCC = manipulateDecayDescriptor.substr(0, 1) + manipulateDecayDescriptor.substr(manipulateDecayDescriptor.size() - 3, 3);
  std::for_each(checkForCC.begin(), checkForCC.end(), [](char& c) { c = ::toupper(c); });

  //Remove the CC check if needed
  if (checkForCC == "[]CC")
  {
    manipulateDecayDescriptor = manipulateDecayDescriptor.substr(1, manipulateDecayDescriptor.size() - 4);
    m_getChargeConjugate = true;
  }

  //Try and find the initial particle
  size_t findMotherEndPoint = manipulateDecayDescriptor.find(decayArrow);
  mother = manipulateDecayDescriptor.substr(0, findMotherEndPoint);
  if (findParticle(mother))
    m_mother_ID = abs(get_pdgcode(mother));
  else
    ddCanBeParsed = false;
  manipulateDecayDescriptor.erase(0, findMotherEndPoint + decayArrow.length());

  //Try and find the intermediates
  while ((pos = manipulateDecayDescriptor.find(startIntermediate)) != std::string::npos)
  {
    size_t findIntermediateStartPoint = manipulateDecayDescriptor.find(startIntermediate, pos);
    size_t findIntermediateEndPoint = manipulateDecayDescriptor.find(endIntermediate, pos);
    std::string intermediateDecay = manipulateDecayDescriptor.substr(pos + 1, findIntermediateEndPoint - (pos + 1));

    intermediate = intermediateDecay.substr(0, intermediateDecay.find(decayArrow));
    if (findParticle(intermediate))
    {
      m_intermediates_ID.push_back(abs(get_pdgcode(intermediate)));
    }
    else
      ddCanBeParsed = false;

    //Now find the daughters associated to this intermediate
    int nDaughters = 0;
    intermediateDecay.erase(0, intermediateDecay.find(decayArrow) + decayArrow.length());
    while ((daughterLocator = intermediateDecay.find(chargeIndicator)) != std::string::npos)
    {
      daughter = intermediateDecay.substr(0, daughterLocator);
      std::string daughterChargeString = intermediateDecay.substr(daughterLocator + 1, 1);
      if (std::find(std::begin(specialTracks), std::end(specialTracks), daughter) != std::end(specialTracks))
      {
        daughter += daughterChargeString;
      }
      if (findParticle(daughter))
      {
        m_daughters_ID.push_back(abs(get_pdgcode(daughter)));
        daughters_charge.push_back(get_charge(daughter)); 
      }
      else
        ddCanBeParsed = false;
      intermediateDecay.erase(0, daughterLocator + 2);
      ++nDaughters;
    }
    manipulateDecayDescriptor.erase(findIntermediateStartPoint, findIntermediateEndPoint + 1 - findIntermediateStartPoint);
    m_nTracksFromIntermediates.push_back(nDaughters);
    m_nTracksFromMother += 1;
  }

  //Now find any remaining reconstructable tracks from the mother
  while ((daughterLocator = manipulateDecayDescriptor.find(chargeIndicator)) != std::string::npos)
  {
    daughter = manipulateDecayDescriptor.substr(0, daughterLocator);
    std::string daughterChargeString = manipulateDecayDescriptor.substr(daughterLocator + 1, 1);
    if (std::find(std::begin(specialTracks), std::end(specialTracks), daughter) != std::end(specialTracks))
    {
      daughter += daughterChargeString;
    }
    if (findParticle(daughter))
    {
      m_daughters_ID.push_back(abs(get_pdgcode(daughter)));
      daughters_charge.push_back(get_charge(daughter));
    }
    else
      ddCanBeParsed = false;
    manipulateDecayDescriptor.erase(0, daughterLocator + 2);
    m_nTracksFromMother += 1;
  }

  unsigned int trackStart = 0;
  unsigned int trackEnd = 0;
  for (unsigned int i = 0; i < m_intermediates_ID.size(); ++i)
  {
    trackStart = trackEnd;
    trackEnd = m_nTracksFromIntermediates[i] + trackStart;

    int vtxCharge = 0;

    for (unsigned int j = trackStart; j < trackEnd; ++j)
    {
      vtxCharge += daughters_charge[j];
    }

    intermediates_charge.push_back(vtxCharge);
  }

  for (unsigned int i = 0; i < m_daughters_ID.size(); ++i) mother_charge += daughters_charge[i];

  m_mother_ID = mother_charge == 0 ? m_mother_ID : mother_charge * m_mother_ID;
  for (unsigned int i = 0; i < m_intermediates_ID.size(); ++i)
  {
    m_intermediates_ID[i] = intermediates_charge[i] == 0 ? m_intermediates_ID[i] : intermediates_charge[i] * m_intermediates_ID[i];
    m_motherDecayProducts.push_back(m_intermediates_ID[i]);
  }
  for (unsigned int i = 0; i < m_daughters_ID.size(); ++i)
  {
    m_daughters_ID[i] = daughters_charge[i] == 0 ? m_daughters_ID[i] : daughters_charge[i] * m_daughters_ID[i];
    if (i >= trackEnd) m_motherDecayProducts.push_back(m_daughters_ID[i]);
  }

  if (ddCanBeParsed)
  {
    if (Verbosity() >= VERBOSITY_MORE) std::cout << "Your decay descriptor can be parsed" << std::endl;
    return 0;
  }
  else
  {
    if (Verbosity() >= VERBOSITY_SOME) std::cout << "Your decay descriptor cannot be parsed, " << Name() << " will not be registered" << std::endl;
    return Fun4AllReturnCodes::DONOTREGISTERSUBSYSTEM;
  }
}

/*
 * Main body for searching your decay
 * Slight issue if you limit the generator decay
 * as decays wont enter the HepMC record
 * need a switch to go to Geant4 record
 */
bool DecayFinder::findDecay(PHCompositeNode* topNode)
{
  bool decayWasFound = false;
  bool reconstructableDecayWasFound = false;
  bool aTrackFailedPT = false;
  bool aTrackFailedETA = false;
  bool aMotherHasPhoton = false;
  bool aMotherHasPi0 = false;
  std::vector<int> correctMotherProducts;
  std::vector<int> positive_motherDecayProducts;

  int n = sizeof(listOfResonantPIDs) / sizeof(listOfResonantPIDs[0]);
  for (unsigned int i = 0; i < m_intermediates_ID.size(); ++i)
    n = deleteElement(listOfResonantPIDs, n, m_intermediates_ID[i]);

  for (unsigned int i = 0; i < m_motherDecayProducts.size(); ++i)
    positive_motherDecayProducts.push_back(std::abs(m_motherDecayProducts[i]));

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
    exit(1);
  }

  if (m_truthinfo && !m_geneventmap) //This should use the truth info container if we have no HepMC record
  {
    if (Verbosity() >= VERBOSITY_SOME)
    {
    	std::cout << "You have a valid G4TruthInfo but no PHHepMCGenEventMap container" << std::endl;
    	std::cout << "You are probably using a particle gun" << std::endl;
    }
  
    PHG4TruthInfoContainer::ConstRange range = m_truthinfo->GetParticleRange();
    for(PHG4TruthInfoContainer::ConstIterator iter = range.first; iter != range.second; ++iter)
    {
      PHG4Particle* g4particle = iter->second;
      int this_pid = m_getChargeConjugate ? abs(g4particle->get_pid()) : g4particle->get_pid();
      if (this_pid == m_mother_ID)
      {
        if (Verbosity() >= VERBOSITY_MAX) std::cout << "parent->pdg_id(): " << g4particle->get_pid() << std::endl;
  
        bool breakOut = false;
        correctMotherProducts.clear();
        decayChain.clear();
        decayChain.push_back(std::make_pair(std::make_pair(g4particle->get_primary_id(), g4particle->get_barcode()), g4particle->get_pid()));
  
        searchGeant4Record(g4particle->get_barcode(), g4particle->get_pid(), positive_motherDecayProducts, breakOut, aMotherHasPhoton, aMotherHasPi0, aTrackFailedPT, aTrackFailedETA, correctMotherProducts);
  
        if (breakOut) continue;    
  
        decayWasFound = compareDecays(m_motherDecayProducts, correctMotherProducts);     
  
        if (decayWasFound)
        {
          m_counter += 1;
          if (aTrackFailedPT && !aTrackFailedETA) m_nCandFail_pT += 1;
          else if (!aTrackFailedPT && aTrackFailedETA) m_nCandFail_eta += 1;
          else if (aTrackFailedPT && aTrackFailedETA) m_nCandFail_pT_and_eta += 1;
          else
          {
            m_nCandReconstructable += 1;
            reconstructableDecayWasFound = true;
            if (m_save_dst) fillDecayNode(topNode, decayChain);
            if (aMotherHasPhoton && !aMotherHasPi0) m_nCandHas_Photon += 1;
            else if (!aMotherHasPhoton && aMotherHasPi0) m_nCandHas_Pi0 += 1;
            else if (aMotherHasPhoton && aMotherHasPi0) m_nCandHas_Photon_and_Pi0 += 1;
            else m_nCandHas_noPhoton_and_noPi0 += 1;
          }
        }
      }
    }
    return reconstructableDecayWasFound;
  }


  for (PHHepMCGenEventMap::ConstReverseIter riter = m_geneventmap->rbegin();
      riter !=  m_geneventmap->rend(); ++riter)
  {
    m_genevt = riter->second;
    if (!m_genevt)
    {
      std::cout << "DecayFinder: Missing node PHHepMCGenEvent" << std::endl;
      return 0;
    }

    HepMC::GenEvent* theEvent = m_genevt->getEvent();

    for (HepMC::GenEvent::particle_const_iterator p = theEvent->particles_begin(); p != theEvent->particles_end(); ++p)
    {
      int this_pid = m_getChargeConjugate ? abs((*p)->pdg_id()) : (*p)->pdg_id();
      if (this_pid == m_mother_ID)
      {
        if (Verbosity() >= VERBOSITY_MAX) std::cout << "parent->pdg_id(): " << (*p)->pdg_id() << std::endl;

        bool breakOut = false;
        correctMotherProducts.clear();
        decayChain.clear();
        decayChain.push_back(std::make_pair(std::make_pair(m_genevt->get_embedding_id(), (*p)->barcode()), (*p)->pdg_id()));

        // Make sure that the mother has a decay in our record
        if (!(*p)->end_vertex()) //Mother has no end vertex, decay volume was limited
        {
          searchGeant4Record((*p)->barcode(), (*p)->pdg_id(), positive_motherDecayProducts, breakOut, aMotherHasPhoton, aMotherHasPi0, aTrackFailedPT, aTrackFailedETA, correctMotherProducts);
        }
        else
        {
          searchHepMCRecord((*p), positive_motherDecayProducts, breakOut, aMotherHasPhoton, aMotherHasPi0, aTrackFailedPT, aTrackFailedETA, correctMotherProducts);
        }

        if (breakOut) continue;

        decayWasFound = compareDecays(m_motherDecayProducts, correctMotherProducts); 

        if (decayWasFound)
        {
          m_counter += 1;
          if (aTrackFailedPT && !aTrackFailedETA) m_nCandFail_pT += 1;
          else if (!aTrackFailedPT && aTrackFailedETA) m_nCandFail_eta += 1;
          else if (aTrackFailedPT && aTrackFailedETA) m_nCandFail_pT_and_eta += 1;
          else
          {
            m_nCandReconstructable += 1;
            reconstructableDecayWasFound = true;
            if (m_save_dst) fillDecayNode(topNode, decayChain);
            if (aMotherHasPhoton && !aMotherHasPi0) m_nCandHas_Photon += 1;
            else if (!aMotherHasPhoton && aMotherHasPi0) m_nCandHas_Pi0 += 1;
            else if (aMotherHasPhoton && aMotherHasPi0) m_nCandHas_Photon_and_Pi0 += 1;
            else m_nCandHas_noPhoton_and_noPi0 += 1;
          }
        }
      }
    }
  }

  return reconstructableDecayWasFound;
}

/*
 * Function to search HepMC record
 * Can switch to Geant4 search if needed
 */
void DecayFinder::searchHepMCRecord(HepMC::GenParticle* particle, std::vector<int> decayProducts, 
                                    bool &breakLoop, bool &hasPhoton, bool &hasPi0, bool &failedPT, bool &failedETA, 
                                    std::vector<int> &actualDecayProducts)
{
  for (HepMC::GenVertex::particle_iterator children = particle->end_vertex()->particles_begin(HepMC::children);
       children != particle->end_vertex()->particles_end(HepMC::children); ++children)
  {
    if (Verbosity() >= VERBOSITY_MAX) std::cout << "--children->pdg_id(): " << (*children)->pdg_id() << std::endl;

    if ((*children)->pdg_id() == 22) hasPhoton = true;
    if ((*children)->pdg_id() == 111) hasPi0 = true;
    if (!m_allowPhotons && (*children)->pdg_id() == 22)
    {
      breakLoop = true;
      break;
    }
    if (!m_allowPi0 && (*children)->pdg_id() == 111)
    {
      breakLoop = true;
      break;
    }

    //This is one of the children we are looking for
    if (std::find(decayProducts.begin(), decayProducts.end(),
                  std::abs((*children)->pdg_id())) != decayProducts.end())
    {
      if (Verbosity() >= VERBOSITY_MAX) std::cout << "This is a child you were looking for" << std::endl;
      //Check if this is an internediate decay that didnt decay in the generator
      std::vector<int> positive_intermediates_ID;
      for (unsigned int i = 0; i < m_intermediates_ID.size(); ++i) positive_intermediates_ID.push_back(abs(m_intermediates_ID[i]));

      if (!(*children)->end_vertex() && std::find(positive_intermediates_ID.begin(), positive_intermediates_ID.end(),
         abs((*children)->pdg_id())) != positive_intermediates_ID.end())
      {
        std::vector<int> requiredIntermediateDecayProducts;
        std::vector<int> positive_requiredIntermediateDecayProducts;
        std::vector<int> actualIntermediateDecayProducts;
        //Which intermediate decay list to we need
        auto it = std::find(positive_intermediates_ID.begin(), positive_intermediates_ID.end(), abs((*children)->pdg_id()));
        int index = std::distance(positive_intermediates_ID.begin(), it);

        unsigned int trackStart = 0, trackStop = 0;
        if (index == 0)
          trackStop = m_nTracksFromIntermediates[0];
        else
        {
          for (int i = 0; i < index; ++i) trackStart += m_nTracksFromIntermediates[i];
          trackStop = trackStart + m_nTracksFromIntermediates[index];
        }
        for (unsigned int i = trackStart; i < trackStop; ++i)
        {
          requiredIntermediateDecayProducts.push_back(m_daughters_ID[i]);
          positive_requiredIntermediateDecayProducts.push_back(abs(m_daughters_ID[i]));
        }

        searchGeant4Record((*children)->barcode(), (*children)->pdg_id(), positive_requiredIntermediateDecayProducts, 
                           breakLoop, hasPhoton, hasPi0, failedPT, failedETA, actualIntermediateDecayProducts);

        bool needThisParticle = compareDecays(requiredIntermediateDecayProducts, actualIntermediateDecayProducts);
        if (needThisParticle)
        {
          actualDecayProducts.push_back((*children)->pdg_id());
          decayChain.push_back(std::make_pair(std::make_pair(m_genevt->get_embedding_id(),(*children)->barcode()), (*children)->pdg_id()));
        }
      }
      else
      {
        bool needThisParticle = checkIfCorrectHepMCParticle((*children), failedPT, failedETA);
        if (needThisParticle)
        {
          actualDecayProducts.push_back((*children)->pdg_id());
          decayChain.push_back(std::make_pair(std::make_pair(m_genevt->get_embedding_id(),(*children)->barcode()), (*children)->pdg_id()));
        }
      }
    }  //Now check if it's part of the resonance list
    else if (std::find(std::begin(listOfResonantPIDs), std::end(listOfResonantPIDs),
                       std::abs((*children)->pdg_id())) != std::end(listOfResonantPIDs))
    {
      if (Verbosity() >= VERBOSITY_MAX) std::cout << "This is a resonance to investigate further" << std::endl;
      if (!(*children)->end_vertex())
      {
        searchGeant4Record((*children)->barcode(), (*children)->pdg_id(), decayProducts, 
                           breakLoop, hasPhoton, hasPi0, failedPT, failedETA, actualDecayProducts);
      }
      else
      {
        for (HepMC::GenVertex::particle_iterator grandchildren = (*children)->end_vertex()->particles_begin(HepMC::children);
             grandchildren != (*children)->end_vertex()->particles_end(HepMC::children); ++grandchildren)
        {
          bool needThisParticle = checkIfCorrectHepMCParticle((*grandchildren), failedPT, failedETA);
          if (needThisParticle)
          {
            actualDecayProducts.push_back((*grandchildren)->pdg_id());
            decayChain.push_back(std::make_pair(std::make_pair(m_genevt->get_embedding_id(),(*grandchildren)->barcode()), (*grandchildren)->pdg_id()));
          }
        }
      }
    }
    else
    {
      breakLoop = true;  //This particle is not in the decay descriptor, stop
      break;
    }
  }
}

/*
 *Function to search Geant4 record
 * Cannot switch to HepMC search (doesn't make sense to)
 */
void DecayFinder::searchGeant4Record(int barcode, int pid, std::vector<int> decayProducts, bool &breakLoop, bool &hasPhoton, bool &hasPi0, bool &failedPT, bool &failedETA, std::vector<int> &actualDecayProducts)
{
  PHG4TruthInfoContainer::ConstRange range = m_truthinfo->GetParticleRange();
  for(PHG4TruthInfoContainer::ConstIterator iter = range.first; iter != range.second; ++iter)
  {
    if (decayChain.size() == 100)
    {
      breakLoop = true;  //Stuck in loop. Sympton not cause!
      break;
    }
    
    PHG4Particle* g4particle = iter->second;
    PHG4Particle* mother = nullptr;
    if (g4particle->get_parent_id() != 0) mother = m_truthinfo->GetParticle(g4particle->get_parent_id());
    else continue;
    if (mother->get_barcode() == barcode && abs(mother->get_pid()) == abs(pid))
    {
      int particleID = g4particle->get_pid();
      if (Verbosity() >= VERBOSITY_MAX) std::cout << "--children->pdg_id(): " << particleID << std::endl;
      if (particleID == 22) hasPhoton = true;
      if (particleID == 111) hasPi0 = true;
      
      if (!m_allowPhotons && particleID == 22)
      {
        breakLoop = true;
        break;
      }
      if (!m_allowPi0 && particleID == 111)
      {
        breakLoop = true;
        break;
      }
      //This is one of the children we are looking for
      if (std::find(decayProducts.begin(), decayProducts.end(),
                    std::abs(particleID)) != decayProducts.end())
      {
        bool needThisParticle = checkIfCorrectGeant4Particle(g4particle, hasPhoton, hasPi0, failedPT, failedETA);
        if (needThisParticle)
        {
          if (Verbosity() >= VERBOSITY_MAX) std::cout << "This is a child you were looking for" << std::endl;
          actualDecayProducts.push_back(particleID);
          decayChain.push_back(std::make_pair(std::make_pair(m_genevt->get_embedding_id(),g4particle->get_barcode()), particleID));
        }
      }  //Now check if it's part of the other resonance list
      else if (m_allowPhotons && particleID == 22) continue;
      else if (m_allowPi0 && particleID == 111) continue;
      else if (std::find(std::begin(listOfResonantPIDs), std::end(listOfResonantPIDs),
               std::abs(particleID)) != std::end(listOfResonantPIDs))
      {
        if (Verbosity() >= VERBOSITY_MAX) std::cout << "This is a resonance to investigate further" << std::endl;
        searchGeant4Record(g4particle->get_barcode(), g4particle->get_pid(), decayProducts,
                           breakLoop, hasPhoton, hasPi0, failedPT, failedETA, actualDecayProducts);
      }
      else
      {
        breakLoop = true;  //This particle is not in the decay descriptor, stop
        break;
      }
    }
  }
}  

/*
 * Checks if this prticle really matches your request
 */
bool DecayFinder::checkIfCorrectHepMCParticle(HepMC::GenParticle* particle, bool& trackFailedPT, bool& trackFailedETA)
{
  bool acceptParticle = false;

  std::vector<int> positive_intermediates_ID;
  for (unsigned int i = 0; i < m_intermediates_ID.size(); ++i) positive_intermediates_ID.push_back(abs(m_intermediates_ID[i]));

  //Check if it is an intermediate or a final track
  if (std::find(positive_intermediates_ID.begin(), positive_intermediates_ID.end(),
                abs(particle->pdg_id())) != positive_intermediates_ID.end())
  {
    std::vector<int> requiredIntermediateDecayProducts;
    std::vector<int> actualIntermediateDecayProducts;
    //Which intermediate decay list to we need
    auto it = std::find(positive_intermediates_ID.begin(), positive_intermediates_ID.end(), abs(particle->pdg_id()));
    int index = std::distance(positive_intermediates_ID.begin(), it);

    unsigned int trackStart = 0, trackStop = 0;
    if (index == 0)
      trackStop = m_nTracksFromIntermediates[0];
    else
    {
      for (int i = 0; i < index; ++i) trackStart += m_nTracksFromIntermediates[i];
      trackStop = trackStart + m_nTracksFromIntermediates[index];
    }

    for (unsigned int i = trackStart; i < trackStop; ++i) requiredIntermediateDecayProducts.push_back(m_daughters_ID[i]);

    for (HepMC::GenVertex::particle_iterator grandchildren = particle->end_vertex()->particles_begin(HepMC::children);
         grandchildren != particle->end_vertex()->particles_end(HepMC::children); ++grandchildren)
    {
      if (Verbosity() >= VERBOSITY_MAX) std::cout << "----grandchildren->pdg_id(): " << (*grandchildren)->pdg_id() << std::endl;

      if (std::find(std::begin(listOfResonantPIDs), std::end(listOfResonantPIDs),
                    std::abs((*grandchildren)->pdg_id())) != std::end(listOfResonantPIDs))
      {
        for (HepMC::GenVertex::particle_iterator greatgrandchildren = (*grandchildren)->end_vertex()->particles_begin(HepMC::children);
             greatgrandchildren != (*grandchildren)->end_vertex()->particles_end(HepMC::children); ++greatgrandchildren)
        {
          if (Verbosity() >= VERBOSITY_MAX) std::cout << "--------greatgrandchildren->pdg_id(): " << (*greatgrandchildren)->pdg_id() << std::endl;

          if (m_allowPhotons && (*greatgrandchildren)->pdg_id() == 22) continue;
          else if (!m_allowPhotons && (*greatgrandchildren)->pdg_id() == 22) break;
          else if (m_allowPi0 && (*greatgrandchildren)->pdg_id() == 111) continue;
          else if (!m_allowPi0 && (*greatgrandchildren)->pdg_id() == 111) break;
          else
          {
            actualIntermediateDecayProducts.push_back((*greatgrandchildren)->pdg_id());
            HepMC::FourVector myFourVector = (*greatgrandchildren)->momentum();
            if (myFourVector.perp() < m_pt_req) trackFailedPT = true;
            if (!isInRange(m_eta_low_req, myFourVector.eta(), m_eta_high_req)) trackFailedETA = true;
          }
        }
      }
      else if (m_allowPhotons && (*grandchildren)->pdg_id() == 22) continue;
      else if (!m_allowPhotons && (*grandchildren)->pdg_id() == 22) break;
      else if (m_allowPi0 && (*grandchildren)->pdg_id() == 111) continue;
      else if (!m_allowPi0 && (*grandchildren)->pdg_id() == 111) break;
      else
      {
        actualIntermediateDecayProducts.push_back((*grandchildren)->pdg_id());
        HepMC::FourVector myFourVector = (*grandchildren)->momentum();
        if (myFourVector.perp() < m_pt_req) trackFailedPT = true;
        if (!isInRange(m_eta_low_req, myFourVector.eta(), m_eta_high_req)) trackFailedETA = true;
      }
    }

    acceptParticle = compareDecays(requiredIntermediateDecayProducts, actualIntermediateDecayProducts);

  }
  else if (particle->pdg_id() == 22) return false;
  else if (particle->pdg_id() == 111) return false;
  else
  {
    if (Verbosity() >= VERBOSITY_MAX) std::cout << "This is a final state track" << std::endl;
    HepMC::FourVector myFourVector = particle->momentum();
    if (myFourVector.perp() < m_pt_req) trackFailedPT = true;
    if (!isInRange(m_eta_low_req, myFourVector.eta(), m_eta_high_req)) trackFailedETA = true;
    acceptParticle = true;
  }

  return acceptParticle;
}

/*
 * Checks if this prticle really matches your request
 */
bool DecayFinder::checkIfCorrectGeant4Particle(PHG4Particle *particle, bool& hasPhoton, bool& hasPi0, bool& trackFailedPT, bool& trackFailedETA)
{
  bool acceptParticle = false;

  std::vector<int> positive_intermediates_ID;
  for (unsigned int i = 0; i < m_intermediates_ID.size(); ++i) positive_intermediates_ID.push_back(abs(m_intermediates_ID[i]));
  //Check if it is an intermediate or a final track
  if (std::find(positive_intermediates_ID.begin(), positive_intermediates_ID.end(),
                abs(particle->get_pid())) != positive_intermediates_ID.end())
  {
    std::vector<int> requiredIntermediateDecayProducts;
    std::vector<int> actualIntermediateDecayProducts;
    //Which intermediate decay list to we need
    auto it = std::find(positive_intermediates_ID.begin(), positive_intermediates_ID.end(), abs(particle->get_pid()));
    int index = std::distance(positive_intermediates_ID.begin(), it);

    unsigned int trackStart = 0, trackStop = 0;
    if (index == 0)
      trackStop = m_nTracksFromIntermediates[0];
    else
    {
      for (int i = 0; i < index; ++i) trackStart += m_nTracksFromIntermediates[i];
      trackStop = trackStart + m_nTracksFromIntermediates[index];
    }

    std::vector<int> positive_intermediateDecayProducts;
    for (unsigned int i = trackStart; i < trackStop; ++i)
    {
      positive_intermediateDecayProducts.push_back(abs(m_daughters_ID[i]));
      requiredIntermediateDecayProducts.push_back(m_daughters_ID[i]);
    }

    bool fakeBreak = false;
    searchGeant4Record(particle->get_barcode(), particle->get_pid(), positive_intermediateDecayProducts, fakeBreak, 
                       hasPhoton, hasPi0, trackFailedPT, trackFailedETA, actualIntermediateDecayProducts);

    acceptParticle = compareDecays(requiredIntermediateDecayProducts, actualIntermediateDecayProducts);

  }
  else if (particle->get_pid() == 22) return false;
  else if (particle->get_pid() == 111) return false;
  else
  {
    if (Verbosity() >= VERBOSITY_MAX) std::cout << "This is a final state track" << std::endl;
    double px = particle->get_px();
    double py = particle->get_py();
    double pz = particle->get_pz();

    double pt = std::sqrt(std::pow(px, 2) + std::pow(py, 2));
    double p = std::sqrt(std::pow(px, 2) + std::pow(py, 2) + std::pow(pz, 2));
    double eta = 0.5*std::log((p + pz)/(p - pz));

    if (pt < m_pt_req) trackFailedPT = true;
    if (!isInRange(m_eta_low_req, eta, m_eta_high_req)) trackFailedETA = true;

    acceptParticle = true;
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
    std::cout << "Printing required decay products: ";
    for (unsigned int i = 0; i < required.size(); ++i) std::cout << required[i] << ", ";
    std::cout << std::endl;
    std::cout << "Printing actual decay products: ";
    for (unsigned int i = 0; i < actual.size(); ++i) std::cout << actual[i] << ", ";
    std::cout << std::endl;
    if (required == actual) std::cout << "*\n* These vectors match\n*\n" << std::endl;
    else std::cout << "*\n* These vectors DONT match\n*\n" << std::endl;
  }

  if (required == actual) accept = true;

  if (m_getChargeConjugate && !accept)
  {
    multiplyVectorByScalarAndSort(required, -1);
    if (required == actual) accept = true;

    if (Verbosity() >= VERBOSITY_MAX)
    {
      std::cout << "Checking CC state" << std::endl;
      if (required == actual) std::cout << "*\n* These vectors match\n*\n" << std::endl;
      else std::cout << "*\n* These vectors DONT match\n*\n" << std::endl;
    }
  }
  return accept;
}

/*
 * Everything below this is helper functions
 */
int DecayFinder::deleteElement(int arr[], int n, int x)
{
  // https://www.geeksforgeeks.org/delete-an-element-from-array-using-two-traversals-and-one-traversal/
  // Search x in array
  int i;
  for (i = 0; i < n; i++)
    if (arr[i] == x)
      break;

  // If x found in array
  if (i < n)
  {
    // reduce size of array and move all
    // elements on space ahead
    n = n - 1;
    for (int j = i; j < n; j++)
      arr[j] = arr[j + 1];
  }

  return n;
}

bool DecayFinder::findParticle(const std::string particle)
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

void DecayFinder::multiplyVectorByScalarAndSort(std::vector<int>& v, int k)
{
  //https://slaystudy.com/c-multiply-vector-by-scalar/
  std::transform(v.begin(), v.end(), v.begin(), [k](const int& c) { return c * k; });
  std::sort(v.begin(), v.end());
}

int DecayFinder::get_pdgcode(const std::string name)
{
  if (findParticle(name))
    return TDatabasePDG::Instance()->GetParticle(name.c_str())->PdgCode();
  else
    return 0;
}

int DecayFinder::get_charge(const std::string name)
{
  if (findParticle(name))
    return TDatabasePDG::Instance()->GetParticle(name.c_str())->Charge()/3;
  else
    return -99;
}

bool DecayFinder::isInRange(float min, float value, float max)
{
  return min <= value && value <= max;
}

int DecayFinder::createDecayNode(PHCompositeNode* topNode)
{
  PHNodeIterator iter(topNode);

  PHCompositeNode* lowerNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!lowerNode)
  {
    lowerNode = new PHCompositeNode("DST");
    topNode->addNode(lowerNode);
    std::cout << "DST node added" << std::endl;
  }

  std::string baseName;

  if (m_container_name.empty())
    //baseName = "decay";
    baseName = Name();
  else
    baseName = m_container_name;

  //Cant have forward slashes in DST or else you make a subdirectory on save!!!
  size_t pos;
  std::string undrscr = "_";
  std::string nothing = "";
  std::map<std::string, std::string> forbiddenStrings;
  forbiddenStrings["/"] = undrscr;
  forbiddenStrings["("] = undrscr;
  forbiddenStrings[")"] = nothing;
  forbiddenStrings["+"] = "plus";
  forbiddenStrings["-"] = "minus";
  forbiddenStrings["*"] = "star";
  forbiddenStrings[" "] = nothing;
  for (auto const& [badString, goodString] : forbiddenStrings)
  {
    while ((pos = baseName.find(badString)) != std::string::npos) baseName.replace(pos, 1, goodString);
  }

  m_nodeName = baseName + "_DecayMap";

  m_decayMap = new DecayFinderContainer_v1();
  PHIODataNode<PHObject>* decayNode = new PHIODataNode<PHObject>(m_decayMap, m_nodeName.c_str(), "PHObject");
  lowerNode->addNode(decayNode);
  std::cout << m_nodeName << " node added" << std::endl;

  return Fun4AllReturnCodes::EVENT_OK;
}

void DecayFinder::fillDecayNode(PHCompositeNode* topNode, Decay decay)
{
  m_decayMap = findNode::getClass<DecayFinderContainer_v1>(topNode, m_nodeName.c_str());
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
  std::cout << "  Number of reconstructable mothers with associated photon and Pi0: " << m_nCandHas_Photon_and_Pi0 << std::endl;
  std::cout << "  Number of reconstructable mothers with no associated photons or Pi0s: " << m_nCandHas_noPhoton_and_noPi0 << std::endl;
  std::cout << "-----------------------------------------------------\n";
  std::cout << std::endl;
}

void DecayFinder::printNode(PHCompositeNode* topNode)
{
  std::cout << "----------------";
  std::cout << " DecayFinderNode: " << m_nodeName << " information ";
  std::cout << "----------------" << std::endl;
  DecayFinderContainer_v1* map = findNode::getClass<DecayFinderContainer_v1>(topNode, m_nodeName.c_str());
  for (DecayFinderContainer_v1::Iter iter = map->begin(); iter != map->end(); ++iter)
  {
    Decay decay = iter->second;
    for (unsigned int i = 0; i < decay.size(); ++i)
    {
      std::cout << "Particle embedding ID: " << decay[i].first.first << ", particle barcode: " << decay[i].first.second << ", particle PDG ID: " << decay[i].second << std::endl;
    }
  }
  std::cout << "--------------------------------------------------------------------------------------------------" << std::endl;
}
