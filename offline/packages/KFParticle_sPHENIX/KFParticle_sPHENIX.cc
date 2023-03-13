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

#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxVertexMap.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/getClass.h>

#include <TFile.h>

#include <KFParticle.h>           // for KFParticle
#include <fun4all/Fun4AllBase.h>  // for Fun4AllBase::VERBOSITY...
#include <fun4all/SubsysReco.h>   // for SubsysReco

#include <cctype>    // for toupper
#include <cmath>     // for sqrt
#include <cstdlib>   // for size_t, exit
#include <iostream>  // for operator<<, endl, basi...
#include <map>       // for map
#include <tuple>     // for tie, tuple

class PHCompositeNode;

namespace TMVA
{
  class Reader;
}

int candidateCounter = 0;

/// KFParticle constructor
KFParticle_sPHENIX::KFParticle_sPHENIX()
  : SubsysReco("KFPARTICLE")
  , m_has_intermediates_sPHENIX(false)
  , m_constrain_to_vertex_sPHENIX(false)
  , m_require_mva(false)
  , m_save_dst(0)
  , m_save_output(1)
  , m_outfile_name("outputData.root")
  , m_outfile(nullptr)
{
}

KFParticle_sPHENIX::KFParticle_sPHENIX(const std::string &name)
  : SubsysReco(name)
  , m_has_intermediates_sPHENIX(false)
  , m_constrain_to_vertex_sPHENIX(false)
  , m_require_mva(false)
  , m_save_dst(0)
  , m_save_output(1)
  , m_outfile_name("outputData.root")
  , m_outfile(nullptr)
{
}

int KFParticle_sPHENIX::Init(PHCompositeNode *topNode)
{
  if (m_save_output && Verbosity() >= VERBOSITY_SOME) std::cout << "Output nTuple: " << m_outfile_name << std::endl;

  if (m_save_dst) createParticleNode(topNode);

  if (m_require_mva)
  {
    TMVA::Reader *reader;
    std::vector<Float_t> MVA_parValues;
    tie(reader, MVA_parValues) = initMVA();
  }

  int returnCode = 0;
  if (!m_decayDescriptor.empty()) returnCode = parseDecayDescriptor();

  return returnCode;
}

int KFParticle_sPHENIX::process_event(PHCompositeNode *topNode)
{
  std::vector<KFParticle> mother, vertex;
  std::vector<std::vector<KFParticle>> daughters, intermediates;
  int nPVs, multiplicity;

  if (!m_use_fake_pv)
  {
    SvtxVertexMap *check_vertexmap = findNode::getClass<SvtxVertexMap>(topNode, m_vtx_map_node_name);
    if (check_vertexmap->size() == 0)
    {
      if (Verbosity() >= VERBOSITY_SOME) std::cout << "KFParticle: Event skipped as there are no vertices" << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }

  SvtxTrackMap *check_trackmap = findNode::getClass<SvtxTrackMap>(topNode, m_trk_map_node_name);
  if (check_trackmap->size() == 0)
  {
    if (Verbosity() >= VERBOSITY_SOME) std::cout << "KFParticle: Event skipped as there are no tracks" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  createDecay(topNode, mother, vertex, daughters, intermediates, nPVs, multiplicity);

  if (!m_has_intermediates_sPHENIX) intermediates = daughters;
  if (!m_constrain_to_vertex_sPHENIX) vertex = mother;

  if (mother.size() != 0)
  {
    for (unsigned int i = 0; i < mother.size(); ++i)
    {
      if (m_save_output && candidateCounter == 0)
      {
        m_outfile = new TFile(m_outfile_name.c_str(), "RECREATE");
        initializeBranches();
      }

      candidateCounter += 1;

      if (m_save_output) fillBranch(topNode, mother[i], vertex[i], daughters[i], intermediates[i], nPVs, multiplicity);
      if (m_save_dst) fillParticleNode(topNode, mother[i], daughters[i], intermediates[i]);

      if (Verbosity() >= VERBOSITY_SOME)
      {
        printParticles(mother[i], vertex[i], daughters[i], intermediates[i], nPVs, multiplicity);
      }
      if (Verbosity() >= VERBOSITY_MORE)
      {
        if (m_save_dst) printNode(topNode);
      }
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int KFParticle_sPHENIX::End(PHCompositeNode * /*topNode*/)
{
  std::cout << "KFParticle_sPHENIX object " << Name() << " finished. Number of candidates: " << candidateCounter << std::endl;

  if (m_save_output && candidateCounter != 0)
  {
    m_outfile->Write();
    m_outfile->Close();
    delete m_outfile;
  }

  return 0;
}

void KFParticle_sPHENIX::printParticles(const KFParticle motherParticle,
                                        const KFParticle chosenVertex,
                                        const std::vector<KFParticle> &daughterParticles,
                                        const std::vector<KFParticle> &intermediateParticles,
                                        const int numPVs, const int numTracks)
{
  std::cout << "\n---------------KFParticle candidate information---------------" << std::endl;

  std::cout << "Mother information:" << std::endl;
  identify(motherParticle);

  if (m_has_intermediates_sPHENIX)
  {
    std::cout << "Intermediate state information:" << std::endl;
    for (unsigned int i = 0; i < intermediateParticles.size(); i++)
    {
      identify(intermediateParticles[i]);
    }
  }

  std::cout << "Final track information:" << std::endl;
  for (unsigned int i = 0; i < daughterParticles.size(); i++)
  {
    identify(daughterParticles[i]);
  }

  if (m_constrain_to_vertex_sPHENIX)
  {
    std::cout << "Primary vertex information:" << std::endl;
    std::cout << "(x,y,z) = (" << chosenVertex.GetX() << " +/- " << sqrt(chosenVertex.GetCovariance(0, 0)) << ", ";
    std::cout << chosenVertex.GetY() << " +/- " << sqrt(chosenVertex.GetCovariance(1, 1)) << ", ";
    std::cout << chosenVertex.GetZ() << " +/- " << sqrt(chosenVertex.GetCovariance(2, 2)) << ") cm\n"
              << std::endl;
  }

  std::cout << "The number of primary vertices is: " << numPVs << std::endl;
  std::cout << "The number of tracks in the event is: " << numTracks << std::endl;

  std::cout << "------------------------------------------------------------\n"
            << std::endl;
}

int KFParticle_sPHENIX::parseDecayDescriptor()
{
  bool ddCanBeParsed = true;

  size_t daughterLocator;

  std::string mother;
  std::string intermediate;
  std::string daughter;

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

  //These tracks require a + or - after their name for TDatabasePDG
  std::string specialTracks[] = {"e", "mu", "pi", "K"};

  std::string manipulateDecayDescriptor = m_decayDescriptor;

  //Remove all white space before we begin
  size_t pos;
  while ((pos = manipulateDecayDescriptor.find(" ")) != std::string::npos) manipulateDecayDescriptor.replace(pos, 1, "");

  //Check for charge conjugate requirement
  std::string checkForCC = manipulateDecayDescriptor.substr(0, 1) + manipulateDecayDescriptor.substr(manipulateDecayDescriptor.size() - 3, 3);
  std::for_each(checkForCC.begin(), checkForCC.end(), [](char &c) { c = ::toupper(c); });

  //Remove the CC check if needed
  if (checkForCC == "[]CC")
  {
    manipulateDecayDescriptor = manipulateDecayDescriptor.substr(1, manipulateDecayDescriptor.size() - 4);
    getChargeConjugate(true);
  }

  //Find the initial particle
  size_t findMotherEndPoint = manipulateDecayDescriptor.find(decayArrow);
  mother = manipulateDecayDescriptor.substr(0, findMotherEndPoint);
  if (!findParticle(mother)) ddCanBeParsed = false;
  manipulateDecayDescriptor.erase(0, findMotherEndPoint + decayArrow.length());

  //Try and find the intermediates
  while ((pos = manipulateDecayDescriptor.find(startIntermediate)) != std::string::npos)
  {
    size_t findIntermediateStartPoint = manipulateDecayDescriptor.find(startIntermediate, pos);
    size_t findIntermediateEndPoint = manipulateDecayDescriptor.find(endIntermediate, pos);
    std::string intermediateDecay = manipulateDecayDescriptor.substr(pos + 1, findIntermediateEndPoint - (pos + 1));

    intermediate = intermediateDecay.substr(0, intermediateDecay.find(decayArrow));
    if (findParticle(intermediate))
      intermediates_name.push_back(intermediate.c_str());
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
        daughters_name.push_back(daughter.c_str());

        if (daughterChargeString == "+")
        {
          daughters_charge.push_back(+1);
        }
        else if (daughterChargeString == "-")
        {
          daughters_charge.push_back(-1);
        }
        else if (daughterChargeString == "0")
        {
          daughters_charge.push_back(0);
        }
        else
        {
          if (Verbosity() >= VERBOSITY_MORE) std::cout << "The charge of " << daughterChargeString << " was not known" << std::endl;
          ddCanBeParsed = false;
        }
      }
      else
        ddCanBeParsed = false;
      intermediateDecay.erase(0, daughterLocator + 2);
      ++nDaughters;
    }
    manipulateDecayDescriptor.erase(findIntermediateStartPoint, findIntermediateEndPoint + 1 - findIntermediateStartPoint);
    m_nTracksFromIntermediates.push_back(nDaughters);
    nTracks += nDaughters;
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
      daughters_name.push_back(daughter.c_str());
      if (daughterChargeString == "+")
      {
        daughters_charge.push_back(+1);
      }
      else if (daughterChargeString == "-")
      {
        daughters_charge.push_back(-1);
      }
      else if (daughterChargeString == "0")
      {
        daughters_charge.push_back(0);
      }
      else
      {
        if (Verbosity() >= VERBOSITY_MORE) std::cout << "The charge of " << daughterChargeString << " was not known" << std::endl;
        ddCanBeParsed = false;
      }
    }
    else
      ddCanBeParsed = false;
    manipulateDecayDescriptor.erase(0, daughterLocator + 2);
    nTracks += 1;
  }

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

    intermediate_list.push_back(std::make_pair(intermediates_name[i], intermediates_charge[i]));
  }

  for (int i = 0; i < nTracks; ++i)
  {
    daughter_list.push_back(std::make_pair(daughters_name[i], daughters_charge[i]));
  }

  setMotherName(mother);
  setNumberOfTracks(nTracks);
  setDaughters(daughter_list);

  if (intermediates_name.size() > 0)
  {
    hasIntermediateStates(true);
    setIntermediateStates(intermediate_list);
    setNumberOfIntermediateStates(intermediates_name.size());
    setNumberTracksFromIntermeditateState(m_nTracksFromIntermediates);
  }

  if (ddCanBeParsed)
  {
    if (Verbosity() >= VERBOSITY_MORE) std::cout << "Your decay descriptor can be parsed" << std::endl;
    return 0;
  }
  else
  {
    if (Verbosity() >= VERBOSITY_SOME) std::cout << "KFParticle: Your decay descriptor, " << Name() << " cannot be parsed"
                                                 << "\nExiting!" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
}
