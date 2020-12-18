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

  virtual ~KFParticle_eventReconstruction(){}

  /**
   * Starts the reconstruction chain
   *
   * @param selectedMother Input a vector and it will be filled with any mother candidates
   * @param selectedVertex Input a vector and it will be filled with the production vertex associated to your mother
   * @param selectedDaughters Input a vector and it will be filled with any tracks associated to your mother
   * @param selectedIntermediates Input a vector and it will be filled with any intermediate states associated to your mother
   */
  void createDecay(PHCompositeNode* topNode, vector<KFParticle>& selectedMother, vector<KFParticle>& selectedVertex,
                   vector<vector<KFParticle>>& selectedDaughters,
                   vector<vector<KFParticle>>& selectedIntermediates,
                   int& nPVs, int& multiplicity);
 
  ///Used to reconstruct simple decays with no intermediate states
  void buildBasicChain(vector<KFParticle>& selectedMotherBasic,
                       vector<KFParticle>& selectedVertexBasic,
                       vector<vector<KFParticle>>& selectedDaughtersBasic,
                       const vector<KFParticle> daughterParticlesBasic,
                       const vector<int> goodTrackIndexBasic,
                       const vector<KFParticle> primaryVerticesBasic);

  ///Used to reconstruct more complicated decays with up to four intermediate states
  void buildChain(vector<KFParticle>& selectedMotherAdv,
                  vector<KFParticle>& selectedVertexAdv,
                  vector<vector<KFParticle>>& selectedDaughtersAdv,
                  vector<vector<KFParticle>>& selectedIntermediatesAdv,
                  vector<KFParticle> daughterParticlesAdv,
                  const vector<int> goodTrackIndexAdv,
                  vector<KFParticle> primaryVerticesAdv);

  ///Basic building block for event reconstruction and selection
  void getCandidateDecay(vector<KFParticle>& selectedMotherCand,
                         vector<KFParticle>& selectedVertexCand,
                         vector<vector<KFParticle>>& selectedDaughtersCand,
                         vector<KFParticle> daughterParticlesCand,
                         vector<vector<int>> goodTracksThatMeetCand,
                         vector<KFParticle> primaryVerticesCand,
                         int n_track_start, int n_track_stop,
                         bool isIntermediate, int intermediateNumber, bool constrainMass);

 protected:
  static const int max_tracks = 99;
  bool m_has_intermediates;
  int m_num_tracks;
  string m_daughter_name_evt[max_tracks];
  int m_daughter_charge_evt[max_tracks];
  int m_intermediate_charge[max_tracks];
  bool m_constrain_to_vertex;
  bool m_constrain_int_mass;
};

#endif  //KFPARTICLESPHENIX_KFPARTICLEEVENTRECONSTRUCTION_H
