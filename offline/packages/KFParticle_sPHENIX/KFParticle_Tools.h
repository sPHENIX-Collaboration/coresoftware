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

#ifndef KFPARTICLESPHENIX_KFPARTICLETOOLS_H
#define KFPARTICLESPHENIX_KFPARTICLETOOLS_H

#include "KFParticle_MVA.h"
#include "KFParticle_particleList.h"

#include <KFParticle.h>

#include <map>
#include <vector>

class PHCompositeNode;

class SvtxVertexMap;
class SvtxTrackMap;
class SvtxVertex;
class SvtxTrack;

class KFParticle_Tools : public KFParticle_particleList, protected KFParticle_MVA
{
 public:
  KFParticle_Tools();

  virtual ~KFParticle_Tools() {}

  KFParticle makeVertex(PHCompositeNode *topNode);

  std::vector<KFParticle> makeAllPrimaryVertices(PHCompositeNode *topNode);

  KFParticle makeParticle(PHCompositeNode *topNode);

  std::vector<KFParticle> makeAllDaughterParticles(PHCompositeNode *topNode);

  int getTracksFromVertex(PHCompositeNode *topNode, KFParticle vertex);

  const bool isGoodTrack(KFParticle particle, std::vector<KFParticle> primaryVertices);

  std::vector<int> findAllGoodTracks(std::vector<KFParticle> daughterParticles, std::vector<KFParticle> primaryVertices);

  std::vector<std::vector<int>> findTwoProngs(std::vector<KFParticle> daughterParticles, std::vector<int> goodTrackIndex, int nTracks);

  std::vector<std::vector<int>> findNProngs(std::vector<KFParticle> daughterParticles,
                                  std::vector<int> goodTrackIndex,
                                  std::vector<std::vector<int>> goodTracksThatMeet,
                                  int nRequiredTracks, unsigned int nProngs);

  std::vector<std::vector<int>> appendTracksToIntermediates(KFParticle intermediateResonances[], std::vector<KFParticle> daughterParticles, std::vector<int> goodTrackIndex, int num_remaining_tracks);

  ///Calculates the cosine of the angle betweent the flight direction and momentum
  float eventDIRA(KFParticle particle, KFParticle vertex);

  float flightDistanceChi2(KFParticle particle, KFParticle vertex);

  std::tuple<KFParticle, bool> buildMother(KFParticle vDaughters[], std::string daughterOrder[], bool isIntermediate, int intermediateNumber, int nTracks, bool constrainMass, float required_vertexID);

  void constrainToVertex(KFParticle &particle, bool &goodCandidate, KFParticle &vertex);

  std::tuple<KFParticle, bool> getCombination(KFParticle vDaughters[], std::string daughterOrder[], KFParticle vertex,
                                         bool constrain_to_vertex, bool isIntermediate, int intermediateNumber, int nTracks, bool constrainMass, float required_vertexID);

  std::vector<std::vector<std::string>> findUniqueDaughterCombinations(int start, int end);

  double calculateEllipsoidRadius(int posOrNeg, double sigma_ii, double sigma_jj, double sigma_ij);

  float calculateEllipsoidVolume(KFParticle particle);

  void identify(KFParticle particle);

 protected:
  static const int max_tracks = 99;
  std::string m_mother_name_Tools;
  int m_num_intermediate_states = -1;
  int m_num_tracks_from_intermediate[max_tracks] = {0};
  std::string m_daughter_name[max_tracks];
  int m_daughter_charge[max_tracks] = {0};
  int m_num_tracks = -1;
  std::string m_intermediate_name[max_tracks];
  float m_min_mass = -1;

  float m_max_mass = -1;

  std::pair<float, float> m_intermediate_mass_range[max_tracks];
  float m_intermediate_min_pt[max_tracks] = {0};
  float m_min_lifetime = -1;

  float m_max_lifetime = -1;

  float m_track_pt = -1;

  float m_track_ptchi2 = -1;

  float m_track_ipchi2 = -1;

  float m_track_chi2ndof = -1;

  float m_comb_DCA = -1;

  float m_vertex_chi2ndof = -1;

  float m_fdchi2 = -1;

  float m_dira_min = -1;

  float m_dira_max = -1;

  float m_mother_pt = -1;

  float m_mother_ipchi2 = -1;

  float m_mva_cut_value = -1;

  bool m_get_charge_conjugate = -1;

  std::string m_vtx_map_node_name;
  std::string m_trk_map_node_name;
  SvtxVertexMap *m_dst_vertexmap = nullptr;
  SvtxTrackMap *m_dst_trackmap = nullptr;
  SvtxVertex *m_dst_vertex = nullptr;
  SvtxTrack *m_dst_track = nullptr;

 private:
  void removeDuplicates(std::vector<double> &v);
  void removeDuplicates(std::vector<int> &v);
  void removeDuplicates(std::vector<std::vector<int>> &v);
  void removeDuplicates(std::vector<std::vector<std::string>> &v);
};

#endif  //KFPARTICLESPHENIX_KFPARTICLETOOLS_H
