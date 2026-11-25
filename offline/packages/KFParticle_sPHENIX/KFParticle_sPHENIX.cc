/*
 * This file is part of KFParticle package
 * Copyright ( C ) 2007-2019 FIAS Frankfurt Institute for Advanced Studies
 *               2007-2019 Goethe University of Frankfurt
 *               2007-2019 Ivan Kisel <I.Kisel@compeng.uni-frankfurt.de>
 *               2007-2019 Maksym Zyzak
 *
 * KFParticle is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * ( at your option ) any later version.
 *
 * KFParticle is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "KFParticle_sPHENIX.h"

#include <centrality/CentralityInfo.h>

#include <globalvertex/MbdVertex.h>
#include <globalvertex/MbdVertexMap.h>
#include <globalvertex/SvtxVertexMap.h>
#include <trackbase_historic/SvtxTrackMap.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/getClass.h>
#include <ffaobjects/EventHeader.h>
#include <ffarawobjects/Gl1Packet.h>

#include <TEntryList.h>
#include <TFile.h>
#include <TLeaf.h>
#include <TSystem.h>
#include <TTree.h>  // for getting the TTree from the file

#include <KFPVertex.h>
#include <KFParticle.h>           // for KFParticle
#include <fun4all/Fun4AllBase.h>  // for Fun4AllBase::VERBOSITY...
#include <fun4all/SubsysReco.h>   // for SubsysReco

#include <ffamodules/CDBInterface.h>  // for accessing the field map file from the CDB
#include <cctype>                     // for toupper
#include <cmath>                      // for sqrt
#include <cstdlib>                    // for size_t, exit
#include <filesystem>
#include <iostream>  // for operator<<, endl, basi...
#include <map>       // for map

class PHCompositeNode;

//int candidateCounter = 0;

/// KFParticle constructor
KFParticle_sPHENIX::KFParticle_sPHENIX() : KFParticle_sPHENIX("KFPARTICLE"){}

KFParticle_sPHENIX::KFParticle_sPHENIX(const std::string &name)
  : SubsysReco(name)
  , KFParticle_nTuple(m_config)  // binds only; m_config is constructed before any of it is read
  , KFParticle_DST(m_config)  // same
  , KFParticle_eventReconstruction(m_config)  // same
  , m_outfile(nullptr)
{
  setToolsInstance(static_cast<KFParticle_Tools *>(this));
}

int KFParticle_sPHENIX::Init(PHCompositeNode *topNode)
{
  m_config.verbosity = Verbosity();
  m_qa_tree_path = Name() + "/DecayTree";  // one directory per module so several KFParticle modules do not clash
  
  if (m_config.saveOutput && !m_config.lightweightNTupleToQA && Verbosity() >= VERBOSITY_SOME)
  {
    std::cout << "Output nTuple: " << m_config.outputName << std::endl;
  }

  if (m_config.saveDST)
  {
    createParticleNode(topNode);
  }

  int returnCode = 0;
  if (!m_config.decayDescriptor.empty())
  {
    returnCode = parseDecayDescriptor();
  }

  if (m_config.usePID)
  {
    init_dEdx_fits();
  }

  if (m_config.bunchCrossingMode == KFParticle_Config::BunchCrossingMode::Mixed && Verbosity() >= VERBOSITY_MORE)
  {
    std::cout << "KFParticle: You have forced to make a mixed-event sample. This overwrites any request to match bunch-crossings between tracks and primary vertices" << std::endl;
  }

  return returnCode;
}

int KFParticle_sPHENIX::InitRun(PHCompositeNode *topNode)
{
  assert(topNode);

  getField();

  if (m_config.useCentrality && !findNode::getClass<CentralityInfo>(topNode, "CentralityInfo"))
  {
    std::cout << "KFParticle - [WARNING] - can't find CentralityInfo node CentralityInfo" << std::endl;
  }

  return 0;
}

int KFParticle_sPHENIX::process_event(PHCompositeNode *topNode)
{ 
  std::vector<KFParticle> mother;
  std::vector<KFParticle> vertex_kfparticle;
  std::vector<std::vector<KFParticle>> daughters;
  std::vector<std::vector<KFParticle>> intermediates;
  int nPVs;

  // Adding BCO Matching  
  auto* evtHeader = findNode::getClass<EventHeader>(topNode, "EventHeader"); // event header node
  auto* gl1packet = findNode::getClass<Gl1Packet>(topNode, "GL1RAWHIT");     // gl1 packet node

  if (!gl1packet)
  {
    gl1packet = findNode::getClass<Gl1Packet>(topNode, "GL1Packet");
  }

  if (evtHeader && gl1packet)
  {
    const int64_t run = evtHeader->get_RunNumber();
    const int64_t evn = evtHeader->get_EvtSequence();
    m_this_event_bco = static_cast<uint64_t>(gl1packet->lValue(0, "BCO"));

    if (Verbosity() >= VERBOSITY_A_LOT)
    {
      std::cout << "Event start | run: " << run << " event: " << evn << " this_event_bco: " << m_this_event_bco << std::endl;
    }

    if (run != m_prev_runNumber || evn != m_prev_eventNumber)
    {

      if (Verbosity() >= VERBOSITY_A_LOT)
      {
        std::cout << "New event detected" << std::endl;
      	std::cout << "Previous event BCO: " << m_prev_event_bco << std::endl;
      }
 
      m_last_event_bco = (run == m_prev_runNumber) ? m_prev_event_bco : -1;
      m_prev_event_bco = m_this_event_bco;

      m_prev_runNumber = run;
      m_prev_eventNumber = evn;

      if (Verbosity() >= VERBOSITY_A_LOT)
      {
        std::cout << "Updated values | last_event_bco: " << m_last_event_bco
                  << " stored_prev_event_bco: " << m_prev_event_bco
                  << std::endl;
      }
    }
  }
  else
  {
	
    if (Verbosity() >= VERBOSITY_MORE)
    {
      std::cout << "KFParticle: EventHeader or GL1 packet not found" << std::endl;
    }
    m_this_event_bco = 0;
    m_last_event_bco = 0;
    m_prev_event_bco = 0;
    m_prev_runNumber = 0;
    m_prev_eventNumber = 0;
  }
  // End BCO matching here

  SvtxTrackMap *check_trackmap = findNode::getClass<SvtxTrackMap>(topNode, m_config.trackMapNodeName);

  if (!check_trackmap || check_trackmap->empty())
  {
    if (Verbosity() >= VERBOSITY_SOME)
    {
      std::cout << "KFParticle: Event skipped as there are no tracks" << std::endl;
    }
    return Fun4AllReturnCodes::EVENT_OK;
  }
  const int multiplicity = check_trackmap->size();

  if (!m_config.useFakePrimaryVertex)
  {
    const bool noVertices = m_config.useMbdVertex
                                ? findNode::getClass<MbdVertexMap>(topNode, "MbdVertexMap")->empty()
                                : findNode::getClass<SvtxVertexMap>(topNode, m_config.vertexMapNodeName)->empty();
    if (noVertices)
    {
      if (Verbosity() >= VERBOSITY_SOME)
      {
        std::cout << "KFParticle: Event skipped as there are no vertices" << std::endl;
      }
      return Fun4AllReturnCodes::EVENT_OK;
    }
  }
  

  createDecay(topNode, mother, vertex_kfparticle, daughters, intermediates, nPVs);
  if (!m_config.hasIntermediates)
  {
    intermediates = daughters;
  }
  if (!m_config.constrainToPrimaryVertex)
  {
    vertex_kfparticle = mother;
  }

  for (unsigned int i = 0; i < mother.size(); ++i)
  {
    
    if (m_config.saveOutput && getCandidateCounter() == 0)
    {
      if (!m_config.lightweightNTupleToQA)
      {
        m_outfile = new TFile(m_config.outputName.c_str(), "RECREATE");
      }
      initializeBranches(topNode);
    }

    incrementCandidateCounter();

    if (m_config.saveOutput)
    {
      set_event_bcos(m_this_event_bco, m_last_event_bco); //filling nTuple for BCO Matching
      fillBranch(topNode, mother[i], vertex_kfparticle[i], daughters[i], intermediates[i]);
    }
    if (m_config.saveDST)
    {
      fillParticleNode(topNode, mother[i], vertex_kfparticle[i], daughters[i], intermediates[i]);
    }

    if (Verbosity() >= VERBOSITY_SOME)
    {
      printParticles(mother[i], vertex_kfparticle[i], daughters[i], intermediates[i], nPVs, multiplicity);
    }
    if (Verbosity() >= VERBOSITY_MORE && m_config.saveDST)
    {
      printNode(topNode);
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int KFParticle_sPHENIX::End(PHCompositeNode * /*topNode*/)
{
  std::cout << "KFParticle_sPHENIX object " << Name() << " finished. Number of candidates: " << getCandidateCounter() << std::endl;

  if (m_config.saveOutput && !m_config.lightweightNTupleToQA && getCandidateCounter() != 0)
  {
    m_outfile->Write();
    m_outfile->Close();
    delete m_outfile;
  }

  return 0;
}

void KFParticle_sPHENIX::printParticles(const KFParticle &motherParticle,
                                        const KFParticle &chosenVertex,
                                        const std::vector<KFParticle> &daughterParticles,
                                        const std::vector<KFParticle> &intermediateParticles,
                                        const int numPVs, const int numTracks) const
{
  std::cout << "\n---------------KFParticle candidate information---------------" << std::endl;

  std::cout << "Mother information:" << std::endl;
  KFParticle_util::identify(motherParticle);

  if (m_config.hasIntermediates)
  {
    std::cout << "Intermediate state information:" << std::endl;
    for (const auto &intermediateParticle : intermediateParticles)
    {
      KFParticle_util::identify(intermediateParticle);
    }
  }

  std::cout << "Final track information:" << std::endl;
  for (const auto &daughterParticle : daughterParticles)
  {
    KFParticle_util::identify(daughterParticle);
  }

  if (m_config.constrainToPrimaryVertex)
  {
    std::cout << "Primary vertex information:" << std::endl;
    std::cout << "(x,y,z) = (" << chosenVertex.GetX() << " +/- " << std::sqrt(chosenVertex.GetCovariance(0, 0)) << ", ";
    std::cout << chosenVertex.GetY() << " +/- " << std::sqrt(chosenVertex.GetCovariance(1, 1)) << ", ";
    std::cout << chosenVertex.GetZ() << " +/- " << std::sqrt(chosenVertex.GetCovariance(2, 2)) << ") cm\n"
              << std::endl;
  }

  std::cout << "The number of primary vertices is: " << numPVs << std::endl;
  std::cout << "The number of tracks in the event is: " << numTracks << std::endl;

  std::cout << "------------------------------------------------------------\n"
            << std::endl;
}

int KFParticle_sPHENIX::parseDecayDescriptor()
{
  std::string mother;
  std::string intermediate;

  std::vector<std::pair<std::string, int>> intermediate_list;
  std::vector<std::string> intermediates_name;
  std::vector<int> intermediates_charge;

  std::vector<std::pair<std::string, int>> daughter_list;
  std::vector<std::string> daughters_name;
  std::vector<int> daughters_charge;

  int nTracks = 0;
  std::vector<int> m_nTracksFromIntermediates;

  std::string decayArrow = "->";
  std::string chargeIndicator = "^";
  std::string startIntermediate = "{";
  std::string endIntermediate = "}";
  static const std::map<std::string, int> chargeLookup = {{"+", +1}, {"-", -1}, {"0", 0}};

  // These tracks require a + or - after their name for TDatabasePDG
  std::string specialTracks[] = {"e", "mu", "pi", "K"};

  std::string manipulateDecayDescriptor = m_config.decayDescriptor;

  // Remove all white space before we begin
  size_t pos;
  while ((pos = manipulateDecayDescriptor.find(' ')) != std::string::npos)
  {
    manipulateDecayDescriptor.replace(pos, 1, "");
  }

  // Check for charge conjugate requirement
  std::string checkForCC = manipulateDecayDescriptor.substr(0, 1) + manipulateDecayDescriptor.substr(manipulateDecayDescriptor.size() - 3, 3);
  std::for_each(checkForCC.begin(), checkForCC.end(), [](char &c)
                { c = ::toupper(c); });

  // Remove the CC check if needed
  if (checkForCC == "[]CC")
  {
    manipulateDecayDescriptor = manipulateDecayDescriptor.substr(1, manipulateDecayDescriptor.size() - 4);
    getChargeConjugate();
  }

  // Find the initial particle
  size_t findMotherEndPoint = manipulateDecayDescriptor.find(decayArrow);
  mother = manipulateDecayDescriptor.substr(0, findMotherEndPoint);
  bool ddCanBeParsed = KFParticle_util::findParticle(mother);
  manipulateDecayDescriptor.erase(0, findMotherEndPoint + decayArrow.length());

  // Consume every "name^charge" daughter in decay, filling daughters_name/daughters_charge; returns how many were read
  auto parseDaughters = [&](std::string &decay)
  {
    int nFound = 0;
    size_t daughterLocator;
    while ((daughterLocator = decay.find(chargeIndicator)) != std::string::npos)
    {
      std::string daughter = decay.substr(0, daughterLocator);
      std::string daughterChargeString = decay.substr(daughterLocator + 1, 1);
      if (std::find(std::begin(specialTracks), std::end(specialTracks), daughter) != std::end(specialTracks))
      {
        daughter += daughterChargeString;
      }
      if (KFParticle_util::findParticle(daughter))
      {
        daughters_name.emplace_back(daughter.c_str());
        auto charge = chargeLookup.find(daughterChargeString);
        if (charge != chargeLookup.end())
        {
          daughters_charge.push_back(charge->second);
        }
        else
        {
          if (Verbosity() >= VERBOSITY_MORE)
          {
            std::cout << "The charge of " << daughterChargeString << " was not known" << std::endl;
          }
          ddCanBeParsed = false;
        }
      }
      else
      {
        ddCanBeParsed = false;
      }
      decay.erase(0, daughterLocator + 2);
      ++nFound;
    }
    return nFound;
  };

  // Try and find the intermediates
  while ((pos = manipulateDecayDescriptor.find(startIntermediate)) != std::string::npos)
  {
    size_t findIntermediateStartPoint = manipulateDecayDescriptor.find(startIntermediate, pos);
    size_t findIntermediateEndPoint = manipulateDecayDescriptor.find(endIntermediate, pos);
    std::string intermediateDecay = manipulateDecayDescriptor.substr(pos + 1, findIntermediateEndPoint - (pos + 1));

    intermediate = intermediateDecay.substr(0, intermediateDecay.find(decayArrow));
    if (KFParticle_util::findParticle(intermediate))
    {
      intermediates_name.emplace_back(intermediate.c_str());
    }
    else
    {
      ddCanBeParsed = false;
    }

    // Now find the daughters associated to this intermediate
    intermediateDecay.erase(0, intermediateDecay.find(decayArrow) + decayArrow.length());
    const int nDaughters = parseDaughters(intermediateDecay);
    manipulateDecayDescriptor.erase(findIntermediateStartPoint, findIntermediateEndPoint + 1 - findIntermediateStartPoint);
    m_nTracksFromIntermediates.push_back(nDaughters);
    nTracks += nDaughters;
  }

  // Now find any remaining reconstructable tracks from the mother
  nTracks += parseDaughters(manipulateDecayDescriptor);

  int trackEnd = 0;
  for (unsigned int i = 0; i < intermediates_name.size(); ++i)
  {
    int trackStart = trackEnd;
    trackEnd = m_nTracksFromIntermediates[i] + trackStart;

    int vtxCharge = 0;

    for (int j = trackStart; j < trackEnd; ++j)
    {
      vtxCharge += daughters_charge[j];
    }

    intermediates_charge.push_back(vtxCharge);

    intermediate_list.emplace_back(intermediates_name[i], intermediates_charge[i]);
  }

  daughter_list.reserve(nTracks);
  for (int i = 0; i < nTracks; ++i)
  {
    daughter_list.emplace_back(daughters_name[i], daughters_charge[i]);
  }

  setMotherName(mother);
  setNumberOfTracks(nTracks);
  setDaughters(daughter_list);

  if (!intermediates_name.empty())
  {
    hasIntermediateStates();
    setIntermediateStates(intermediate_list);
    setNumberOfIntermediateStates(intermediates_name.size());
    setNumberTracksFromIntermeditateState(m_nTracksFromIntermediates);
  }

  if (ddCanBeParsed)
  {
    if (Verbosity() >= VERBOSITY_MORE)
    {
      std::cout << "Your decay descriptor can be parsed" << std::endl;
    }
    return 0;
  }
  
  if (Verbosity() >= VERBOSITY_SOME)
  {
    std::cout << "KFParticle: Your decay descriptor, " << Name() << " cannot be parsed"
	      << "\nExiting!" << std::endl;
  }
  return Fun4AllReturnCodes::ABORTRUN;
 
}

void KFParticle_sPHENIX::getField()
{
  //This sweeps the sPHENIX magnetic field map from some point radially then grabs the first event that passes the selection
  m_config.magFieldFile = std::filesystem::exists(m_config.magFieldFile) ? m_config.magFieldFile : CDBInterface::instance()->getUrl(m_config.magFieldFile);

  if (Verbosity() > 0)
  {
    std::cout << PHWHERE << ": using fieldmap : " << m_config.magFieldFile << std::endl;
  }

  TFile *fin = new TFile(m_config.magFieldFile.c_str());
  TTree *fieldmap = (TTree *) fin->Get("fieldmap");

  float Bz = 0.;
  unsigned int r = 0.;
  float z = 0.;

  double arc = M_PI/2;
  unsigned int n = 0;

  while (Bz == 0)
  {
    if (n == 4)
    {
      ++r;
    }

    if (r == 3) //Dont go too far out radially
    {
      ++z;
    }

    n = n & 0x3U; //Constrains n from 0 to 3
    r = r & 0x2U;

    double x = r*std::cos(n*arc);
    double y = r*std::sin(n*arc);

    std::string sweep = "x == " + std::to_string(x) + " && y == " + std::to_string(y) + " && z == " + std::to_string(z);

    fieldmap->Draw(">>elist", sweep.c_str(), "entrylist");
    TEntryList *elist = (TEntryList*)gDirectory->Get("elist");
    if (elist->GetEntry(0))
    {
      TLeaf *fieldValue = fieldmap->GetLeaf("bz");
      fieldValue->GetBranch()->GetEntry(elist->GetEntry(0));
      Bz = fieldValue->GetValue();
    }

    ++n;

    if (r == 0) // No point in rescanning (0,0)
    {
      ++r;
      n = 0;
    }
  }
  // The actual unit of KFParticle is in kilo Gauss (kG), which is equivalent to 0.1 T, instead of Tesla (T). The positive value indicates the B field is in the +z direction
  Bz *= 10;  // Factor of 10 to convert the B field unit from kG to T
  KFParticle::SetField((double) Bz);

  fieldmap->Delete();
  fin->Close();
}
