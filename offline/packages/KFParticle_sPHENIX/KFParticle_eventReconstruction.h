/**
 * This file is part of KFParticle package
 * Copyright (C) 2007-2019 FIAS Frankfurt Institute for Advanced Studies
 *               2007-2019 Goethe University of Frankfurt
 *               2007-2019 Ivan Kisel <I.Kisel@compeng.uni-frankfurt.de>
 *               2007-2019 Maksym Zyzak
 *
 * KFParticle is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * KFParticle is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef KFPARTICLESPHENIX_KFPARTICLEEVENTRECONSTRUCTION_H
#define KFPARTICLESPHENIX_KFPARTICLEEVENTRECONSTRUCTION_H

#include "KFParticle_Tools.h"

#include <KFParticle.h>

#include <vector>

class PHCompositeNode;

class KFParticle_eventReconstruction : public KFParticle_Tools
{
 public:
  KFParticle_eventReconstruction();

  virtual ~KFParticle_eventReconstruction() {}

  /**
   * Starts the reconstruction chain
   *
   * @param selectedMother Input a vector and it will be filled with any mother candidates
   * @param selectedVertex Input a vector and it will be filled with the production vertex associated to your mother
   * @param selectedDaughters Input a vector and it will be filled with any tracks associated to your mother
   * @param selectedIntermediates Input a vector and it will be filled with any intermediate states associated to your mother
   */
  void createDecay(PHCompositeNode* topNode, std::vector<KFParticle>& selectedMother, std::vector<KFParticle>& selectedVertex,
                   std::vector<std::vector<KFParticle>>& selectedDaughters,
                   std::vector<std::vector<KFParticle>>& selectedIntermediates,
                   int& nPVs, int& multiplicity);

  ///Used to reconstruct simple decays with no intermediate states
  void buildBasicChain(std::vector<KFParticle>& selectedMotherBasic,
                       std::vector<KFParticle>& selectedVertexBasic,
                       std::vector<std::vector<KFParticle>>& selectedDaughtersBasic,
                       const std::vector<KFParticle>& daughterParticlesBasic,
                       const std::vector<int>& goodTrackIndexBasic,
                       const std::vector<KFParticle>& primaryVerticesBasic);

  ///Used to reconstruct more complicated decays with up to four intermediate states
  void buildChain(std::vector<KFParticle>& selectedMotherAdv,
                  std::vector<KFParticle>& selectedVertexAdv,
                  std::vector<std::vector<KFParticle>>& selectedDaughtersAdv,
                  std::vector<std::vector<KFParticle>>& selectedIntermediatesAdv,
                  const std::vector<KFParticle>& daughterParticlesAdv,
                  const std::vector<int>& goodTrackIndexAdv,
                  const std::vector<KFParticle>& primaryVerticesAdv);

  ///Basic building block for event reconstruction and selection
  void getCandidateDecay(std::vector<KFParticle>& selectedMotherCand,
                         std::vector<KFParticle>& selectedVertexCand,
                         std::vector<std::vector<KFParticle>>& selectedDaughtersCand,
                         std::vector<KFParticle> daughterParticlesCand,
                         std::vector<std::vector<int>> goodTracksThatMeetCand,
                         std::vector<KFParticle> primaryVerticesCand,
                         int n_track_start, int n_track_stop,
                         bool isIntermediate, int intermediateNumber, bool constrainMass);

  ///Method to chose best candidate from a selection of common SV's
  int selectBestCombination(bool PVconstraint, bool isAnInterMother,
                            std::vector<KFParticle> possibleCandidates,
                            std::vector<KFParticle> possibleVertex);

  KFParticle createFakePV();

 protected:
  bool m_constrain_to_vertex;
  bool m_constrain_int_mass;
  bool m_use_fake_pv;

  //private:
};

#endif  //KFPARTICLESPHENIX_KFPARTICLEEVENTRECONSTRUCTION_H
