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

/*****************/
/* Cameron Dean  */
/*   LANL 2020   */
/* cdean@bnl.gov */
/*****************/

#include "KFParticle_eventReconstruction.h"
#include "KFParticle_Tools.h"
#include "KFParticle_util.h"

//sPHENIX stuff
#include <trackbase_historic/SvtxTrack.h>

// KFParticle stuff
#include <KFParticle.h>

#include <algorithm>
#include <cassert>
#include <cstdlib>   // for std::abs(int)
#include <iterator>  // for begin, distance, end
#include <memory>   // for allocator_traits<>::value_type
#include <numeric>  // for iota
#include <set>
#include <string>   // for string
#include <tuple>    // for tie, tuple

#include <iostream>

namespace
{
  // pi+ pdgid = 211, pi- pdgid = -211
  // e+ pdgid = -11, e- pdgid = 11
  // mu+ pdgid = -13, mu- pdgid = 13
  // tau+ pdgid = -15, tau- pdgid = 15
  const std::set<int> chargedLeptonPDGIDs = {11, 13, 15};
}  // namespace

namespace
{
  // Copy of track with the mass and PDG ID for hypothesis pid. A charged track gets q * pid
  // (sign flipped for charged leptons, so e-/mu- keep a positive PDG code); a neutral track gets pid.
  KFParticle withAssignedPID(KFParticle &track, int pid)
  {
    double mass = -1;
    int pdg = 0;
    if ((Int_t) track.GetQ() != 0)
    {
      mass = KFParticle_util::getParticleMass((Int_t) track.GetQ() * pid);
      pdg = (Int_t) track.GetQ() * pid;

      if (chargedLeptonPDGIDs.contains(std::abs(pdg)))
      {
        pdg *= -1;
      }
    }
    else
    {
      mass = KFParticle_util::getParticleMass(pid);
      pdg = pid;
    }
    KFParticle particle;
    particle.Create(track.Parameters(), track.CovarianceMatrix(), (Int_t) track.GetQ(), mass);
    particle.NDF() = track.GetNDF();
    particle.Chi2() = track.GetChi2();
    particle.SetId(track.Id());
    particle.SetPDG(pdg);
    return particle;
  }
}  // namespace

/// KFParticle constructor
KFParticle_eventReconstruction::KFParticle_eventReconstruction(const KFParticle_Config& config)
  : KFParticle_Tools(config)
  , m_config(config)
{
}

void KFParticle_eventReconstruction::createDecay(PHCompositeNode* topNode, std::vector<KFParticle>& selectedMother, std::vector<KFParticle>& selectedVertex,
                                                 std::vector<std::vector<KFParticle>>& selectedDaughters,
                                                 std::vector<std::vector<KFParticle>>& selectedIntermediates,
                                                 int& nPVs)
{
  std::vector<KFParticle> primaryVertices;
  if (m_config.useFakePrimaryVertex)
  {
    primaryVertices.push_back(createFakePV());
  }
  else
  {
    primaryVertices = makeAllPrimaryVertices(topNode, m_config.vertexMapNodeName);
  }

  std::vector<KFParticle> daughterParticles = makeAllDaughterParticles(topNode);

  nPVs = primaryVertices.size();

  std::vector<int> goodTrackIndex(daughterParticles.size());
  std::iota(goodTrackIndex.begin(), goodTrackIndex.end(), 0);  // every track index, 0..n-1

  if (m_config.verbosity >= 10)
  {
    KFParticle_util::printSelectionCheck("Number of daughters passing state selection", daughterParticles.size());
    KFParticle_util::printSelectionCheck("Number of daughters passing track selection", goodTrackIndex.size());
    KFParticle_util::printSelectionCheck("Number of PVs passing selection", primaryVertices.size());
  }

  if (!m_config.hasIntermediates)
  {
    buildBasicChain(selectedMother, selectedVertex, selectedDaughters, daughterParticles, goodTrackIndex, primaryVertices, topNode);
  }
  else
  {
    buildChain(selectedMother, selectedVertex, selectedDaughters, selectedIntermediates, daughterParticles, goodTrackIndex, primaryVertices, topNode);
  }
}

/*
 *  This function is used to build a basic n-body decay without any intermediate particles such as D's or J/psi's
 */
void KFParticle_eventReconstruction::buildBasicChain(std::vector<KFParticle>& selectedMotherBasic,
                                                     std::vector<KFParticle>& selectedVertexBasic,
                                                     std::vector<std::vector<KFParticle>>& selectedDaughtersBasic,
                                                     const std::vector<KFParticle>& daughterParticlesBasic,
                                                     const std::vector<int>& goodTrackIndexBasic,
                                                     const std::vector<KFParticle>& primaryVerticesBasic, PHCompositeNode* topNode)
{
  std::vector<std::vector<int>> goodTracksThatMeet = findTwoProngs(daughterParticlesBasic, goodTrackIndexBasic, m_config.numTracks, primaryVerticesBasic);
  for (int p = 3; p < m_config.numTracks + 1; ++p)
  {
    goodTracksThatMeet = findNProngs(daughterParticlesBasic, goodTrackIndexBasic, goodTracksThatMeet, m_config.numTracks, p, primaryVerticesBasic);
  }

  if (m_config.verbosity >= 10)
  { 
    KFParticle_util::printSelectionCheck("Number of SVs passing selection", goodTracksThatMeet.size());
  }

  getCandidateDecay(selectedMotherBasic, selectedVertexBasic, selectedDaughtersBasic, daughterParticlesBasic,
                    goodTracksThatMeet, primaryVerticesBasic, 0, m_config.numTracks, false, 0, true, topNode);
}

/*
 *  This function is used to build a more complicated decay with intermediate particles such as D's or J/psi's
 */
void KFParticle_eventReconstruction::buildChain(std::vector<KFParticle>& selectedMotherAdv,
                                                std::vector<KFParticle>& selectedVertexAdv,
                                                std::vector<std::vector<KFParticle>>& selectedDaughtersAdv,
                                                std::vector<std::vector<KFParticle>>& selectedIntermediatesAdv,
                                                const std::vector<KFParticle>& daughterParticlesAdv,
                                                const std::vector<int>& goodTrackIndexAdv,
                                                const std::vector<KFParticle>& primaryVerticesAdv, PHCompositeNode* topNode)
{
  int track_start = 0;
  int track_stop = m_config.numTracksFromIntermediate[0];
  std::vector<KFParticle> goodCandidates;
  std::vector<KFParticle> goodVertex;
  std::vector<std::vector<KFParticle>> goodDaughters(m_config.numTracks);
  std::vector<std::vector<KFParticle>> goodIntermediates(m_config.numIntermediateStates);
  std::vector<std::vector<KFParticle>> potentialIntermediates(m_config.numIntermediateStates);
  std::vector<std::vector<std::vector<KFParticle>>> potentialDaughters(m_config.numIntermediateStates);
  for (int i = 0; i < m_config.numIntermediateStates; ++i)
  {
    std::vector<KFParticle> vertices;
    std::vector<std::vector<int>> goodTracksThatMeet = findTwoProngs(daughterParticlesAdv, goodTrackIndexAdv, m_config.numTracksFromIntermediate[i], primaryVerticesAdv);
    for (int p = 3; p <= m_config.numTracksFromIntermediate[i]; ++p)
    {
      goodTracksThatMeet = findNProngs(daughterParticlesAdv,
                                       goodTrackIndexAdv,
                                       goodTracksThatMeet,
                                       m_config.numTracksFromIntermediate[i], p, primaryVerticesAdv);
    }

    if (m_config.verbosity >= 10)
    { 
      KFParticle_util::printSelectionCheck("Number of SVs passing selection", goodTracksThatMeet.size());
    }

    getCandidateDecay(potentialIntermediates[i], vertices, potentialDaughters[i], daughterParticlesAdv,
                      goodTracksThatMeet, primaryVerticesAdv, track_start, track_stop, true, i, m_config.constrainIntermediateMasses, topNode);
    if (i + 1 >= m_config.numIntermediateStates)
    {
      break;
    }
    track_start = track_stop;
    track_stop += m_config.numTracksFromIntermediate[i + 1];
    if (track_stop > m_config.numTracks)
    {
      break;
    }
  }
  int num_tracks_used_by_intermediates = 0;
  for (int i = 0; i < m_config.numIntermediateStates; ++i)
  {
    num_tracks_used_by_intermediates += m_config.numTracksFromIntermediate[i];
  }

  int num_remaining_tracks = m_config.numTracks - num_tracks_used_by_intermediates;
  // Number of potential intermediates found
  const unsigned int num_pot_inter_a = potentialIntermediates[0].size();
  const unsigned int num_pot_inter_b = m_config.numIntermediateStates < 2 ? 1 : potentialIntermediates[1].size();  // Ensure the code inside the loop below is executed
  const unsigned int num_pot_inter_c = m_config.numIntermediateStates < 3 ? 1 : potentialIntermediates[2].size();
  const unsigned int num_pot_inter_d = m_config.numIntermediateStates < 4 ? 1 : potentialIntermediates[3].size();

  for (unsigned int a = 0; a < num_pot_inter_a; ++a)
  {
    for (unsigned int b = 0; b < num_pot_inter_b; ++b)
    {
      for (unsigned int c = 0; c < num_pot_inter_c; ++c)
      {
        for (unsigned int d = 0; d < num_pot_inter_d; ++d)
        {
          const unsigned int matchIterators[4] = {a, b, c, d};

          std::vector<KFParticle> finalTracks = potentialDaughters[0][a];
          for (int j = 1; j < m_config.numIntermediateStates; ++j)
          {
            finalTracks.insert(finalTracks.end(), potentialDaughters[j][matchIterators[j]].begin(), potentialDaughters[j][matchIterators[j]].end());
          }

          // Skip combinations where intermediates share a track
          std::vector<int> trackIds;
          trackIds.reserve(finalTracks.size());
          for (const auto& track : finalTracks)
          {
            trackIds.push_back(track.Id());
          }
          std::sort(trackIds.begin(), trackIds.end());
          if (std::adjacent_find(trackIds.begin(), trackIds.end()) != trackIds.end())
          {
            continue;
          }

          int num_mother_decay_products = m_config.numIntermediateStates + num_remaining_tracks;
          assert(num_mother_decay_products > 0);
          std::vector<KFParticle> motherDecayProducts(num_mother_decay_products);
          for (int i = 0; i < m_config.numIntermediateStates; ++i)
          {
            motherDecayProducts[i] = potentialIntermediates[i][matchIterators[i]];
          }

          // If there are daughter tracks coming from the mother not an intermediate, need to ensure that the intermeditate decay tracks aren't used again
          std::vector<int> goodTrackIndexAdv_withoutIntermediates = goodTrackIndexAdv;
          for (int m = 0; m < num_tracks_used_by_intermediates; ++m)
          {
            int trackID_to_remove = finalTracks[m].Id();
            int trackElement_to_remove = -1;

            auto it = std::find_if(daughterParticlesAdv.begin(), daughterParticlesAdv.end(),
                                   [&trackID_to_remove](const KFParticle& obj)
                                   { return obj.Id() == trackID_to_remove; });

            if (it != daughterParticlesAdv.end())
            {
              trackElement_to_remove = std::distance(daughterParticlesAdv.begin(), it);
            }

            goodTrackIndexAdv_withoutIntermediates.erase(remove(goodTrackIndexAdv_withoutIntermediates.begin(),
                                                                goodTrackIndexAdv_withoutIntermediates.end(), trackElement_to_remove),
                                                         goodTrackIndexAdv_withoutIntermediates.end());
          }

          float required_unique_vertexID = 0;
          for (int n = 0; n < m_config.numIntermediateStates; ++n)
          {
            required_unique_vertexID += m_config.intermediateCharges[n] * KFParticle_util::getParticleMass(m_config.intermediateNames[n]);
          }

          std::vector<std::vector<int>> uniqueCombinations;
          std::vector<std::vector<int>> listOfTracksToAppend;

          if (num_remaining_tracks != 0)
          {
            for (int i = num_tracks_used_by_intermediates; i < m_config.numTracks; ++i)
            {
              required_unique_vertexID += m_config.daughterCharges[i] * KFParticle_util::getParticleMass(m_config.daughterNames[i]);
            }

            uniqueCombinations = findUniqueDaughterCombinations(num_tracks_used_by_intermediates, m_config.numTracks);  // Unique comb of remaining trackIDs

            listOfTracksToAppend = appendTracksToIntermediates(motherDecayProducts.data(), daughterParticlesAdv, goodTrackIndexAdv_withoutIntermediates, num_remaining_tracks, primaryVerticesAdv);

            for (auto& uniqueCombination : uniqueCombinations)
            {
              for (const auto& element_of_intermediate : m_config.intermediateNames)
              {
                uniqueCombination.insert(begin(uniqueCombination), KFParticle_util::getParticleID(element_of_intermediate));
              }
            }
          }
          else
          {
            std::vector<int> intermediate_id;
            intermediate_id.reserve(m_config.intermediateNames.size());
            for (const auto& element_of_intermediate : m_config.intermediateNames)
            {
              intermediate_id.push_back(KFParticle_util::getParticleID(element_of_intermediate));
            }
            uniqueCombinations.push_back(intermediate_id);
            listOfTracksToAppend.push_back({0});
          }

          for (auto& n_tracks : listOfTracksToAppend)
          {
            for (int n_trackID = 0; n_trackID < num_remaining_tracks; ++n_trackID)
            {
              int mDP_trackElem = m_config.numIntermediateStates + n_trackID;
              int dP_trackElem = n_tracks[n_trackID];
              motherDecayProducts[mDP_trackElem] = daughterParticlesAdv[dP_trackElem];
            }

            for (auto& uniqueCombination : uniqueCombinations)
            {
              for (const auto& i_pv : primaryVerticesAdv)
              {
                auto [candidate, isGood] = getCombination(motherDecayProducts.data(), uniqueCombination.data(), i_pv,
                                                          m_config.constrainToPrimaryVertex, false, 0, num_mother_decay_products, m_config.constrainIntermediateMasses, required_unique_vertexID, topNode);

                if (isGood)
                {
                  // The charge and species of every final-state track must match the decay descriptor. The mother
                  // check only sees the intermediate states (e.g. a neutral Lambda0), so on its own it cannot tell
                  // Xi- -> Lambda0 (p+ pi-) pi- from Xi- -> Lambda0 (p- pi+) pi-
                  std::vector<std::pair<int, int>> assigned;  // (charge, |PDG ID|) of each final-state track
                  assigned.reserve(m_config.numTracks);
                  for (int k = 0; k < num_tracks_used_by_intermediates; ++k)
                  {
                    assigned.emplace_back((int) finalTracks[k].GetQ(), std::abs(finalTracks[k].GetPDG()));
                  }
                  for (int k = 0; k < num_remaining_tracks; ++k)
                  {
                    const int slot = m_config.numIntermediateStates + k;
                    assigned.emplace_back((int) motherDecayProducts[slot].GetQ(), std::abs(uniqueCombination[slot]));
                  }
                  isGood = finalStateMatchesDecay(assigned);
                }

                if (isGood)
                {
                  goodCandidates.push_back(candidate);
                  if (m_config.constrainToPrimaryVertex)
                  {
                    goodVertex.push_back(i_pv);
                  }
                  for (int k = 0; k < m_config.numIntermediateStates; ++k)
                  {
                    goodIntermediates[k].push_back(motherDecayProducts[k]);
                  }
                  for (int k = 0; k < num_tracks_used_by_intermediates; ++k)
                  {
                    goodDaughters[k].push_back(finalTracks[k]);
                  }
                  for (int k = 0; k < num_remaining_tracks; ++k)
                  {  // Need to deal with track mass and PID assignment for extra tracks
                    int trackArrayID = k + m_config.numIntermediateStates;
                    goodDaughters[k + num_tracks_used_by_intermediates].push_back(withAssignedPID(motherDecayProducts[trackArrayID], uniqueCombination[trackArrayID]));
                  }
                }
              }
            }

            if (!goodCandidates.empty())
            {
              int bestCombinationIndex = selectBestCombination(m_config.constrainToPrimaryVertex, false, goodCandidates, goodVertex);

              selectedMotherAdv.push_back(goodCandidates[bestCombinationIndex]);
              if (m_config.constrainToPrimaryVertex)
              {
                selectedVertexAdv.push_back(goodVertex[bestCombinationIndex]);
              }
              std::vector<KFParticle> intermediates;
              intermediates.reserve(m_config.numIntermediateStates);
              for (int i = 0; i < m_config.numIntermediateStates; ++i)
              {
                intermediates.push_back(goodIntermediates[i][bestCombinationIndex]);
              }
              selectedIntermediatesAdv.push_back(intermediates);
              std::vector<KFParticle> particles;
              particles.reserve(m_config.numTracks);
              for (int i = 0; i < m_config.numTracks; ++i)
              {
                particles.push_back(goodDaughters[i][bestCombinationIndex]);
              }
              selectedDaughtersAdv.push_back(particles);
            }
            goodCandidates.clear();
            goodVertex.clear();
            for (auto& v : goodIntermediates)
            {
              v.clear();
            }
            for (auto& v : goodDaughters)
            {
              v.clear();
            }
          }
        }  // Close forth intermediate
      }    // Close third intermediate
    }      // Close second intermediate
  }        // Close first intermediate
}

void KFParticle_eventReconstruction::getCandidateDecay(std::vector<KFParticle>& selectedMotherCand,
                                                       std::vector<KFParticle>& selectedVertexCand,
                                                       std::vector<std::vector<KFParticle>>& selectedDaughtersCand,
                                                       const std::vector<KFParticle>& daughterParticlesCand,
                                                       const std::vector<std::vector<int>>& goodTracksThatMeetCand,
                                                       const std::vector<KFParticle>& primaryVerticesCand,
                                                       int n_track_start, int n_track_stop,
                                                       bool isIntermediate, int intermediateNumber, bool constrainMass, PHCompositeNode* topNode)
{
  int nTracks = n_track_stop - n_track_start;
  std::vector<std::vector<int>> uniqueCombinations = findUniqueDaughterCombinations(n_track_start, n_track_stop);
  std::vector<KFParticle> goodCandidates;
  std::vector<KFParticle> goodVertex;
  std::vector<std::vector<KFParticle>> goodDaughters(nTracks);
  bool fixToPV = m_config.constrainToPrimaryVertex && !isIntermediate;

  float required_unique_vertexID = 0;
  for (int i = n_track_start; i < n_track_stop; ++i)
  {
    required_unique_vertexID += m_config.daughterCharges[i] * KFParticle_util::getParticleMass(m_config.daughterNames[i]);
  }

  for (const auto& i_comb : goodTracksThatMeetCand)  // Loop over all good track combinations
  {
    std::vector<KFParticle> daughterTracks(nTracks);

    for (int i_track = 0; i_track < nTracks; ++i_track)
    {
      daughterTracks[i_track] = daughterParticlesCand[i_comb[i_track]];
    }  // Build array of the good tracks in that combination

    for (auto& uniqueCombination : uniqueCombinations)  // Loop over unique track PID assignments
    {
      for (const auto& pv : primaryVerticesCand)  // Loop over all PVs in the event
      {
        int* PDGIDofFirstParticleInCombination = uniqueCombination.data();
        auto [candidate, isGood] = getCombination(daughterTracks.data(), PDGIDofFirstParticleInCombination, pv, m_config.constrainToPrimaryVertex,
                                                  isIntermediate, intermediateNumber, nTracks, constrainMass, required_unique_vertexID, topNode);
        if (isIntermediate && isGood)
        {
          float min_PV_dca = 0;
          float min_PV_dca_stddev = 0;
          float min_PV_dca_xy = 0;
          float min_PV_dca_stddev_xy = 0;
          calcMinPV_DCA(candidate, primaryVerticesCand, min_PV_dca, min_PV_dca_stddev);
          calcMinPV_DCA(candidate, primaryVerticesCand, min_PV_dca_xy, min_PV_dca_stddev_xy, false);
          const int n = intermediateNumber;
          isGood = KFParticle_util::isInRange(m_config.intermediateMinPVDCA[n], min_PV_dca, m_config.intermediateMaxPVDCA[n]) &&
                   KFParticle_util::isInRange(m_config.intermediateMinPVDCAStdDev[n], min_PV_dca_stddev, m_config.intermediateMaxPVDCAStdDev[n]) &&
                   KFParticle_util::isInRange(m_config.intermediateMinPVDCAXY[n], min_PV_dca_xy, m_config.intermediateMaxPVDCAXY[n]) &&
                   KFParticle_util::isInRange(m_config.intermediateMinPVDCAStdDevXY[n], min_PV_dca_stddev_xy, m_config.intermediateMaxPVDCAStdDevXY[n]);
        }

        if (isGood)
        {
          goodCandidates.push_back(candidate);
          goodVertex.push_back(pv);
          for (int i = 0; i < nTracks; ++i)
          {
            goodDaughters[i].push_back(withAssignedPID(daughterTracks[i], PDGIDofFirstParticleInCombination[i]));
          }
        }
      }
    }

    if (!goodCandidates.empty())
    {
      int bestCombinationIndex = selectBestCombination(fixToPV, isIntermediate, goodCandidates, goodVertex);

      selectedMotherCand.push_back(goodCandidates[bestCombinationIndex]);
      if (fixToPV)
      {
        selectedVertexCand.push_back(goodVertex[bestCombinationIndex]);
      }
      std::vector<KFParticle> particles;
      particles.reserve(nTracks);
      for (int i = 0; i < nTracks; ++i)
      {
        particles.push_back(goodDaughters[i][bestCombinationIndex]);
      }
      selectedDaughtersCand.push_back(particles);

      goodCandidates.clear();
      goodVertex.clear();
      for (auto& v : goodDaughters)
      {
        v.clear();
      }
    }
  }
}

int KFParticle_eventReconstruction::selectBestCombination(bool PVconstraint, bool isAnInterMother,
                                                          std::vector<KFParticle> possibleCandidates,
                                                          std::vector<KFParticle> possibleVertex) const
{
  const bool byMassError = !PVconstraint || isAnInterMother || m_config.selectMotherByMassError;
  auto score = [&](unsigned int i) -> float
  {
    if (byMassError)
    {
      return possibleCandidates[i].GetErrMass();
    }
    return m_config.use2DMatching ? possibleCandidates[i].GetDeviationFromVertexXY(possibleVertex[i])
                                   : possibleCandidates[i].GetDeviationFromVertex(possibleVertex[i]);
  };

  int bestCombinationIndex = 0;
  float bestScore = score(0);
  for (unsigned int i = 1; i < possibleCandidates.size(); ++i)
  {
    const float thisScore = score(i);
    if (thisScore < bestScore)
    {
      bestScore = thisScore;
      bestCombinationIndex = i;
    }
  }
  return bestCombinationIndex;
}

KFParticle KFParticle_eventReconstruction::createFakePV()
{
  float f_vertexParameters[6] = {0};

  float f_vertexCovariance[21] = {0};

  KFParticle kfp_vertex;
  kfp_vertex.Create(f_vertexParameters, f_vertexCovariance, 0, -1);
  kfp_vertex.NDF() = 0;
  kfp_vertex.Chi2() = 0;
  kfp_vertex.SetId(0);
  return kfp_vertex;
}
