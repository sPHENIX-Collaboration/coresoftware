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

#ifndef KFParticle_Tools_H
#define KFParticle_Tools_H

//sPHENIX stuff
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxVertexMap.h>
#include <trackbase_historic/SvtxVertex.h>
#include <phool/getClass.h>
#include <KFParticle_MVA.h>

//ROOT stuff
#include "TMatrixD.h"

//C++ stuff
#include <iostream>
#include <map>
#include <vector>
#include <iomanip>
#include <cmath>

class PHCompositeNode;
class SvtxVertexMap;
class SvtxTrackMap;
class SvtxVertex;
class SvtxTrack;

class KFPTrack;
class KFParticle;
class KFParticleSIMD;
class KFPVertex;
class KFParticleDatabase;

class KFParticle_Tools : protected KFParticle_MVA 
{
 public:
  
  KFParticle_Tools();

  ~KFParticle_Tools();

  KFPVertex makeVertex(PHCompositeNode *topNode);

  std::vector<KFPVertex> makeAllPrimaryVertices(PHCompositeNode *topNode);

  KFPTrack makeTrack(PHCompositeNode *topNode);

  KFParticle makeParticle(PHCompositeNode *topNode, int massHypothesis);

  std::vector<KFParticle> makeAllDaughterParticles(PHCompositeNode *topNode);

  void createDecay(PHCompositeNode *topNode, std::vector<KFParticle>& selectedMother, std::vector<KFParticle>& selectedVertex,
                                             std::vector<KFParticle>& selectedDaughters_1, std::vector<KFParticle>& selectedDaughters_2, 
                                             std::vector<KFParticle>& selectedDaughters_3, std::vector<KFParticle>& selectedDaughters_4,
                                             int& nPVs, int& multiplicity);

  bool isGoodTrack(KFParticle particle, std::vector<KFPVertex> primaryVertices);

  std::vector<int> findAllGoodTracks(std::vector<KFParticle> daughterParticles, std::vector<KFPVertex> primaryVertices);

  std::vector<std::vector<int>> findTwoProngs(std::vector<KFParticle> daughterParticles, std::vector<int> goodTrackIndex);

  std::vector<std::vector<int>> findThreeProngs(std::vector<KFParticle> daughterParticles, std::vector<int> goodTrackIndex, std::vector<std::vector<int>> goodTracksThatMeet);

  std::vector<std::vector<int>> findFourProngs(std::vector<KFParticle> daughterParticles, std::vector<int> goodTrackIndex, std::vector<std::vector<int>> goodTracksThatMeet);

  float eventDIRA(KFParticle particle, KFPVertex vertex);

  float flightDistanceChi2(KFParticle particle, KFPVertex vertex);

  std::tuple<KFParticle, bool> getCombination(KFParticle vDaughters[], std::string daughterOrder[], KFPVertex vertex);

  std::vector<std::vector<std::string>> findUniqueDaughterCombinations();

 
 protected:
  
  int m_num_tracks;
  std::string m_daughter_one;
  std::string m_daughter_two;
  std::string m_daughter_three;
  std::string m_daughter_four;
  int m_daughter_one_charge;
  int m_daughter_two_charge;
  int m_daughter_three_charge;
  int m_daughter_four_charge;
  float m_min_mass;
  float m_max_mass;
  float m_min_lifetime;
  float m_max_lifetime;
  float m_track_pt;
  float m_track_ptchi2;
  float m_track_ipchi2;
  float m_track_chi2ndof;
  float m_comb_DCA;
  float m_vertex_chi2ndof;
  float m_fdchi2;
  float m_dira_min;
  float m_dira_max;
  float m_mother_pt;
  int   m_mother_charge;
  float m_mva_cut_value;
  SvtxVertexMap *m_dst_vertexmap;
  SvtxTrackMap *m_dst_trackmap;
  SvtxVertex *m_dst_vertex;
  SvtxTrack *m_dst_track;


 private:
  
  bool chargeChecker( KFParticle vDaughters[] );
 
  float returnPDGMass( const int pdgIndex);

  void removeDuplicates(std::vector<int> &v);
  void removeDuplicates(std::vector<std::vector<int>> &v);
  void removeDuplicates(std::vector<std::vector<std::string>> &v);
};

#endif //KFParticle_Tools_H

