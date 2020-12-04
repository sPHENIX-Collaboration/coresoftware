/*
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

#ifndef KFParticle_eventReconstruction_H
#define KFParticle_eventReconstruction_H

//sPHENIX stuff
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxVertexMap.h>
#include <trackbase_historic/SvtxVertex.h>
#include <phool/getClass.h>
#include <KFParticle_Tools.h>

//ROOT stuff
#include "TMatrixD.h"

//C++ stuff
#include <iostream>
#include <map>
#include <vector>
#include <iomanip>
#include <cmath>

using namespace std;


class PHCompositeNode;
class KFParticle;

class KFParticle_eventReconstruction : public KFParticle_Tools 
{
 public:
  
  KFParticle_eventReconstruction();

  ~KFParticle_eventReconstruction();

  void createDecay(PHCompositeNode *topNode, vector<KFParticle>& selectedMother, vector<KFParticle>& selectedVertex,
                                             vector<vector<KFParticle>>& selectedDaughters,
                                             vector<vector<KFParticle>>& selectedIntermediates,
                                             int& nPVs, int& multiplicity);

  void buildBasicChain(vector<KFParticle>& selectedMother, 
                       vector<KFParticle>& selectedVertex,
                       vector<vector<KFParticle>>& selectedDaughters, 
                       vector<KFParticle> daughterParticles,
                       vector<int>  goodTrackIndex,
                       vector<KFParticle> primaryVertices);

  void buildChain(vector<KFParticle>& selectedMother,
                  vector<KFParticle>& selectedVertex,
                  vector<vector<KFParticle>>& selectedDaughters,
                  vector<vector<KFParticle>>& selectedIntermediates,
                  vector<KFParticle> daughterParticles,
                  vector<int> goodTrackIndex,
                  vector<KFParticle> primaryVertices);

  void getCandidateDecay(vector<KFParticle>& selectedMother,
                         vector<KFParticle>& selectedVertex,
                         vector<vector<KFParticle>>& selectedDaughters,
                         vector<KFParticle> daughterParticles,
                         vector<vector<int>> goodTracksThatMeet,
                         vector<KFParticle> primaryVertices,
                         int n_track_start, int n_track_stop,
                         bool isIntermediate, int intermediateNumber, bool constrainMass);

 protected:
   
  bool m_has_intermediates;
  int m_num_tracks;
  string m_daughter_name_evt[99];
  int m_daughter_charge_evt[99];
  int m_intermediate_charge[99];
  bool m_constrain_to_vertex;
  bool m_constrain_int_mass;

};

#endif //KFParticle_eventReconstruction_H

