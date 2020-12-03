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
#include <KFParticle_particleList.h>

//ROOT stuff
#include "TMatrixD.h"

//C++ stuff
#include <iostream>
#include <map>
#include <vector>
#include <iomanip>
#include <cmath>


class SvtxEvalStack;

class PHCompositeNode;
class SvtxVertexMap;
class SvtxTrackMap;
class SvtxVertex;
class SvtxTrack;

class KFPTrack;
class KFParticle;
class KFParticleSIMD;
class KFParticle;
class KFParticleDatabase;

class KFParticle_Tools : public KFParticle_particleList, protected KFParticle_MVA 
{
 public:
  
  KFParticle_Tools();

  ~KFParticle_Tools();

  KFParticle makeVertex(PHCompositeNode *topNode);

  std::vector<KFParticle> makeAllPrimaryVertices(PHCompositeNode *topNode);

  KFParticle makeParticle(PHCompositeNode *topNode);

  std::vector<KFParticle> makeAllDaughterParticles(PHCompositeNode *topNode);

  void createDecay(PHCompositeNode *topNode, std::vector<KFParticle>& selectedMother, std::vector<KFParticle>& selectedVertex,
                                             std::vector<std::vector<KFParticle>>& selectedDaughters,
                                             std::vector<std::vector<KFParticle>>& selectedIntermediates,
                                             int& nPVs, int& multiplicity);

  void buildBasicChain(std::vector<KFParticle>& selectedMother, 
                       std::vector<KFParticle>& selectedVertex,
                       std::vector<std::vector<KFParticle>>& selectedDaughters, 
                       std::vector<KFParticle> daughterParticles,
                       std::vector<int>  goodTrackIndex,
                       std::vector<KFParticle> primaryVertices);

  void buildChain(std::vector<KFParticle>& selectedMother,
                  std::vector<KFParticle>& selectedVertex,
                  std::vector<std::vector<KFParticle>>& selectedDaughters,
                  std::vector<std::vector<KFParticle>>& selectedIntermediates,
                  std::vector<KFParticle> daughterParticles,
                  std::vector<int> goodTrackIndex,
                  std::vector<KFParticle> primaryVertices);

  void getCandidateDecay(std::vector<KFParticle>& selectedMother,
                         std::vector<KFParticle>& selectedVertex,
                         std::vector<std::vector<KFParticle>>& selectedDaughters,
                         std::vector<KFParticle> daughterParticles,
                         std::vector<std::vector<int>> goodTracksThatMeet,
                         std::vector<KFParticle> primaryVertices,
                         int n_track_start, int n_track_stop,
                         bool isIntermediate, int intermediateNumber, bool constrainMass);

  int getTracksFromVertex( PHCompositeNode *topNode, KFParticle vertex );

  const bool isGoodTrack(KFParticle particle, std::vector<KFParticle> primaryVertices);

  std::vector<int> findAllGoodTracks(std::vector<KFParticle> daughterParticles, std::vector<KFParticle> primaryVertices);

  std::vector<std::vector<int>> findTwoProngs(std::vector<KFParticle> daughterParticles, std::vector<int> goodTrackIndex, int nTracks);

  std::vector<std::vector<int>> findNProngs( std::vector<KFParticle> daughterParticles, 
                                             std::vector<int> goodTrackIndex, 
                                             std::vector<std::vector<int>> goodTracksThatMeet, 
                                             int nRequiredTracks, unsigned int nProngs );

  std::vector<std::vector<int>> appendTracksToIntermediates( KFParticle intermediateResonances[], std::vector<KFParticle> daughterParticles, std::vector<int> goodTrackIndex, int num_remaining_tracks);

  float eventDIRA(KFParticle particle, KFParticle vertex);

  float flightDistanceChi2(KFParticle particle, KFParticle vertex);

  std::tuple<KFParticle, bool> buildMother( KFParticle vDaughters[], std::string daughterOrder[], bool isIntermediate, int intermediateNumber, int nTracks, bool constrainMass, float required_vertexID );

  void constrainToVertex(KFParticle& particle, bool& goodCandidate, KFParticle& vertex );

  std::tuple<KFParticle, bool> getCombination(KFParticle vDaughters[], std::string daughterOrder[], KFParticle vertex, 
                                              bool constrain_to_vertex, bool isIntermediate, int intermediateNumber, int nTracks, bool constrainMass, float required_vertexID );

  std::vector<std::vector<std::string>> findUniqueDaughterCombinations( int start, int end );

  double calculateEllipsoidRadius( int posOrNeg, double sigma_ii, double sigma_jj, double sigma_ij );
  
  float calculateEllipsoidVolume( KFParticle particle );
 
 protected:
   
  std::string m_mother_name_Tools;
  bool m_has_intermediates;
  int m_num_tracks;
  int m_num_intermediate_states;
  int m_num_tracks_from_intermediate[99];
  std::string m_daughter_name[99];
  int m_daughter_charge[99];
  std::string m_intermediate_name[99];
  int m_intermediate_charge[99];
  float m_min_mass;
  float m_max_mass;
  std::pair<float, float> m_intermediate_mass_range[99];
  float m_intermediate_min_pt[99];
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
  float m_mother_ipchi2;
  float m_mva_cut_value;
  bool m_constrain_to_vertex;
  bool m_constrain_int_mass;
  bool m_get_charge_conjugate;
  std::string m_vtx_map_node_name;
  SvtxVertexMap *m_dst_vertexmap;
  SvtxTrackMap *m_dst_trackmap;
  SvtxVertex *m_dst_vertex;
  SvtxTrack *m_dst_track;


 private:
 
    SvtxEvalStack *m_svtx_evalstack;
  void removeDuplicates(std::vector<double> &v);
  void removeDuplicates(std::vector<int> &v);
  void removeDuplicates(std::vector<std::vector<int>> &v);
  void removeDuplicates(std::vector<std::vector<std::string>> &v);
};

#endif //KFParticle_Tools_H

