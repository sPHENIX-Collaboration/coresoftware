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

#include "KFParticle_Config.h"

#include <globalvertex/MbdVertex.h>
#include <globalvertex/MbdVertexMap.h>
#include <trackbase_historic/TrackSeed.h>
#include <trackbase_historic/TrackSeedContainer.h>

#include <dedxfitter/bethe_bloch.h> // included here so inline functions are defined on user end

#include <KFParticle.h>

#include <TF1.h>

#include <limits>
#include <map>      // for map
#include <string>   // for string
#include <tuple>    // for tuple
#include <utility>  // for pair
#include <vector>

class PHCompositeNode;

class SvtxVertexMap;
class SvtxTrackMap;
class SvtxVertex;
class SvtxTrack;
class SvtxTrack_v5;
class GlobalVertexMap;
class GlobalVertex;
class TrkrClusterContainer;
class PHG4TpcGeomContainer;

class KFParticle_Tools
{
 public:
  /// Standalone use (QA, helper objects): reads a default configuration
  KFParticle_Tools();

  /// Used by KFParticle_sPHENIX. The config is owned there and is not yet constructed when this runs:
  /// store the reference only, never read through it here.
  explicit KFParticle_Tools(const KFParticle_Config &config);

  virtual ~KFParticle_Tools() = default;

  KFParticle makeVertex(PHCompositeNode *topNode);

  std::vector<KFParticle> makeAllPrimaryVertices(PHCompositeNode *topNode, const std::string &vertexMapName);

  KFParticle makeParticle(PHCompositeNode *topNode);

  std::vector<KFParticle> makeAllDaughterParticles(PHCompositeNode *topNode);

  /*const*/ bool isGoodTrack(const KFParticle &particle, const std::vector<KFParticle> &primaryVertices);

  bool isTPConlyAtZeroCrossing(const KFParticle &track);

  int calcMinPV_DCA(const KFParticle &track, const std::vector<KFParticle> &PVs, float &minimumPV_DCA, float &minimumPV_DCA_stddev, bool do3D = true);

  std::vector<int> findAllGoodTracks(const std::vector<KFParticle> &daughterParticles);//, const std::vector<KFParticle> &primaryVertices);

  std::vector<std::vector<int>> findTwoProngs(const std::vector<KFParticle> &daughterParticles, std::vector<int> goodTrackIndex, int nTracks, const std::vector<KFParticle> &primaryVertices);

  std::vector<std::vector<int>> findNProngs(const std::vector<KFParticle> &daughterParticles,
                                            const std::vector<int> &goodTrackIndex,
                                            std::vector<std::vector<int>> goodTracksThatMeet,
                                            int nRequiredTracks, unsigned int nProngs, const std::vector<KFParticle> &primaryVertices);

  std::vector<std::vector<int>> appendTracksToIntermediates(KFParticle intermediateResonances[], const std::vector<KFParticle> &daughterParticles, const std::vector<int> &goodTrackIndex, int num_remaining_tracks, const std::vector<KFParticle> &primaryVertices);

  std::tuple<KFParticle, bool> buildMother(KFParticle vDaughters[], int daughterOrder[], bool isIntermediate, int intermediateNumber, int nTracks, bool constrainMass, float required_vertexID, PHCompositeNode *topNode);

  void constrainToVertex(KFParticle &particle, bool &goodCandidate, KFParticle &vertex);

  std::tuple<KFParticle, bool> getCombination(KFParticle vDaughters[], int daughterOrder[], KFParticle vertex,
                                              bool constrain_to_vertex, bool isIntermediate, int intermediateNumber, int nTracks, bool constrainMass, float required_vertexID, PHCompositeNode *topNode);

  std::vector<std::vector<int>> findUniqueDaughterCombinations(int start, int end);

  bool isInRange(float min, float value, float max);

  bool finalStateMatchesDecay(const std::vector<std::pair<int, int>> &assigned);

  float get_dEdx(PHCompositeNode *topNode, const KFParticle &daughter);

  void init_dEdx_fits();

  double get_dEdx_fitValue(float momentum, int PID);

  bool vertexToleranceCheck(float vertexMap_pos[3], float kfp_vertex_pos[3]);

  bool checkTrackAndVertexMatch(KFParticle vDaughters[], int nTracks, const KFParticle &vertex);

  std::vector<int> getParticleCrossings(const KFParticle &particle);

 protected:
  TF1 *f_pion_plus{nullptr};
  TF1 *f_kaon_plus{nullptr};
  TF1 *f_proton_plus{nullptr};
  TF1 *f_pion_minus{nullptr};
  TF1 *f_kaon_minus{nullptr};
  TF1 *f_proton_minus{nullptr};

  std::map<int, TF1 *> pidMap;

  std::map<int, std::vector<int>> m_intermediate_crossings;

  int m_next_intermediate_id{-2};

  GlobalVertexMap *m_dst_globalvertexmap{nullptr};
  GlobalVertex *m_dst_globalvertex{nullptr};
  MbdVertexMap *m_dst_mbdvertexmap{nullptr};
  MbdVertex *m_dst_mbdvertex{nullptr};
  SvtxTrackMap *m_dst_trackmap{nullptr};
  SvtxTrack *m_dst_track{nullptr};
  SvtxVertexMap *m_dst_vertexmap{nullptr};
  SvtxVertex *m_dst_vertex{nullptr};
  TrkrClusterContainer *m_cluster_map{nullptr};
  PHG4TpcGeomContainer *m_geom_container{nullptr};

  void removeDuplicates(std::vector<double> &v);
  void removeDuplicates(std::vector<int> &v);
  void removeDuplicates(std::vector<std::vector<int>> &v);
  void removeDuplicates(std::vector<std::vector<std::string>> &v);

  void printSelectionCheck(const std::string &parameter, float min, float val, float max);
  void printSelectionCheck(const std::string &start, const std::string &accept, const std::string &reject, const std::string &end, bool equality);
  void printSelectionCheck(const std::string &info, unsigned int value);
  std::string accept_colour = "32";
  std::string reject_colour = "31";

 private:
  const KFParticle_Config &m_config;
};

#endif  // KFPARTICLESPHENIX_KFPARTICLETOOLS_H
