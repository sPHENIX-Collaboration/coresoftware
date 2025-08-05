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

#include <globalvertex/MbdVertex.h>
#include <globalvertex/MbdVertexMap.h>

#include <KFParticle.h>

#include <TF1.h>

#include <limits>
#include <string>   // for string
#include <tuple>    // for tuple
#include <utility>  // for pair
#include <vector>

class PHCompositeNode;

class SvtxVertexMap;
class SvtxTrackMap;
class SvtxVertex;
class SvtxTrack;
class GlobalVertexMap;
class GlobalVertex;
class TrkrClusterContainer;
class PHG4TpcCylinderGeomContainer;

class KFParticle_Tools : protected KFParticle_MVA
{
 public:
  KFParticle_Tools();

  ~KFParticle_Tools() override = default;

  KFParticle makeVertex(PHCompositeNode *topNode);

  std::vector<KFParticle> makeAllPrimaryVertices(PHCompositeNode *topNode, const std::string &vertexMapName);

  KFParticle makeParticle(PHCompositeNode *topNode);

  std::vector<KFParticle> makeAllDaughterParticles(PHCompositeNode *topNode);

  void getTracksFromBC(PHCompositeNode *topNode, const int &bunch_crossing, const std::string &vertexMapName, int &nTracks, int &nPVs);

  int getTracksFromVertex(PHCompositeNode *topNode, const KFParticle &vertex, const std::string &vertexMapName);

  /*const*/ bool isGoodTrack(const KFParticle &particle, const std::vector<KFParticle> &primaryVertices);

  int calcMinIP(const KFParticle &track, const std::vector<KFParticle> &PVs, float &minimumIP, float &minimumIPchi2, bool do3D = true);

  std::vector<int> findAllGoodTracks(const std::vector<KFParticle> &daughterParticles, const std::vector<KFParticle> &primaryVertices);

  std::vector<std::vector<int>> findTwoProngs(std::vector<KFParticle> daughterParticles, std::vector<int> goodTrackIndex, int nTracks);

  std::vector<std::vector<int>> findNProngs(std::vector<KFParticle> daughterParticles,
                                            const std::vector<int> &goodTrackIndex,
                                            std::vector<std::vector<int>> goodTracksThatMeet,
                                            int nRequiredTracks, unsigned int nProngs);

  std::vector<std::vector<int>> appendTracksToIntermediates(KFParticle intermediateResonances[], const std::vector<KFParticle> &daughterParticles, const std::vector<int> &goodTrackIndex, int num_remaining_tracks);

  /// Calculates the cosine of the angle betweent the flight direction and momentum
  float eventDIRA(const KFParticle &particle, const KFParticle &vertex, bool do3D = true);

  float flightDistanceChi2(const KFParticle &particle, const KFParticle &vertex);

  std::tuple<KFParticle, bool> buildMother(KFParticle vDaughters[], int daughterOrder[], bool isIntermediate, int intermediateNumber, int nTracks, bool constrainMass, float required_vertexID, PHCompositeNode* topNode);

  void constrainToVertex(KFParticle &particle, bool &goodCandidate, KFParticle &vertex);

  std::tuple<KFParticle, bool> getCombination(KFParticle vDaughters[], int daughterOrder[], KFParticle vertex,
                                              bool constrain_to_vertex, bool isIntermediate, int intermediateNumber, int nTracks, bool constrainMass, float required_vertexID, PHCompositeNode* topNode);

  std::vector<std::vector<int>> findUniqueDaughterCombinations(int start, int end);

  double calculateEllipsoidRadius(int posOrNeg, double sigma_ii, double sigma_jj, double sigma_ij);

  float calculateEllipsoidVolume(const KFParticle &particle);

  float calculateJT(const KFParticle &mother, const KFParticle &daughter);

  bool isInRange(float min, float value, float max);

  bool findParticle(const std::string &particle);

  int getParticleID(const std::string &particle);

  float getParticleMass(const std::string &particle);
  float getParticleMass(const int PDGID);

  void identify(const KFParticle &particle);

  float get_dEdx(PHCompositeNode *topNode, const KFParticle &daughter);

  void init_dEdx_fits();

  double get_dEdx_fitValue(float momentum, int PID);

  bool checkTrackAndVertexMatch(KFParticle vDaughters[], int nTracks, KFParticle vertex);

  void set_dont_use_global_vertex(bool set_variable){ m_dont_use_global_vertex = set_variable; }

 protected:
  std::string m_mother_name_Tools;
  int m_num_intermediate_states {-1};
  std::vector<int> m_num_tracks_from_intermediate;
  std::vector<std::string> m_daughter_name;
  std::vector<int> m_daughter_charge;
  int m_num_tracks {-1};

  bool m_has_intermediates;
  std::vector<std::string> m_intermediate_name;
  std::vector<int> m_intermediate_charge;
  std::vector<std::pair<float, float>> m_intermediate_mass_range;
  std::vector<float> m_intermediate_min_pt;
  std::vector<float> m_intermediate_min_dira;
  std::vector<float> m_intermediate_min_fdchi2;
  std::vector<float> m_intermediate_min_ip_xy;
  std::vector<float> m_intermediate_max_ip_xy;
  std::vector<float> m_intermediate_min_ipchi2_xy;
  std::vector<float> m_intermediate_max_ipchi2_xy;
  std::vector<float> m_intermediate_min_ip;
  std::vector<float> m_intermediate_max_ip;
  std::vector<float> m_intermediate_min_ipchi2;
  std::vector<float> m_intermediate_max_ipchi2;
  std::vector<float> m_intermediate_vertex_volume;

  bool m_use_PID{false};
  float m_dEdx_band_width {0.2}; //Fraction of expected dE/dx
  
  TF1 *f_pion_plus{nullptr};
  TF1 *f_kaon_plus{nullptr};
  TF1 *f_proton_plus{nullptr};
  TF1 *f_pion_minus{nullptr};
  TF1 *f_kaon_minus{nullptr};
  TF1 *f_proton_minus{nullptr};

  std::map<int, TF1*> pidMap;

  float m_min_mass {-1};

  float m_max_mass {-1};

  float m_min_decayTime_xy {-1000};

  float m_max_decayTime_xy {std::numeric_limits<float>::max()};

  float m_min_decayLength_xy {-1000};

  float m_max_decayLength_xy {std::numeric_limits<float>::max()};

  float m_min_decayTime {-1000};

  float m_max_decayTime {std::numeric_limits<float>::max()};

  float m_min_decayLength {-1000};

  float m_max_decayLength {std::numeric_limits<float>::max()};

  float m_mother_min_decay_time_significance {-1};

  float m_mother_min_decay_length_significance {-1};

  float m_mother_min_decay_length_xy_significance {-1};

  float m_track_pt {-1};

  float m_track_ptchi2 {std::numeric_limits<float>::max()};

  float m_track_ip_xy {-1};

  float m_track_ipchi2_xy {-1};

  float m_track_ip {-1};

  float m_track_ipchi2 {-1};

  float m_track_chi2ndof {std::numeric_limits<float>::max()};

  int m_nMVTXStates {2};

  int m_nINTTStates {1};

  int m_nTPCStates {20};

  int m_nTPOTStates {0};

  float m_comb_DCA_xy {std::numeric_limits<float>::max()};

  float m_comb_DCA {std::numeric_limits<float>::max()};

  float m_vertex_chi2ndof {std::numeric_limits<float>::max()};

  float m_fdchi2 {-1};

  float m_dira_xy_min {-1};

  float m_dira_xy_max {1};

  float m_dira_min {-1};

  float m_dira_max {1};

  float m_mother_pt {-1};

  float m_mother_ip {std::numeric_limits<float>::max()};

  float m_mother_ipchi2 {std::numeric_limits<float>::max()};

  float m_mother_ip_xy {std::numeric_limits<float>::max()};

  float m_mother_ipchi2_xy {std::numeric_limits<float>::max()};

  float m_mother_vertex_volume {std::numeric_limits<float>::max()};

  float m_mva_cut_value {-1};

  bool m_get_charge_conjugate {false};

  bool m_extrapolateTracksToSV {true};

  bool m_allowZeroMassTracks {false};

  bool m_use_2D_matching_tools {false};

  float m_min_radial_SV = -1.;

  bool m_bunch_crossing_zero_only {false};  

  bool m_require_bunch_crossing_match {true};

  bool m_use_mbd_vertex {false};

  bool m_dont_use_global_vertex {false};

  bool m_require_track_and_vertex_match {false};

  std::string m_vtx_map_node_name;
  std::string m_trk_map_node_name;
  GlobalVertexMap *m_dst_globalvertexmap {nullptr};
  GlobalVertex *m_dst_globalvertex {nullptr};
  MbdVertexMap *m_dst_mbdvertexmap {nullptr};
  MbdVertex *m_dst_mbdvertex {nullptr};
  SvtxTrackMap *m_dst_trackmap {nullptr};
  SvtxTrack *m_dst_track {nullptr};
  SvtxVertexMap *m_dst_vertexmap {nullptr};
  SvtxVertex *m_dst_vertex {nullptr};
  TrkrClusterContainer *m_cluster_map {nullptr};
  PHG4TpcCylinderGeomContainer *m_geom_container {nullptr};

  void removeDuplicates(std::vector<double> &v);
  void removeDuplicates(std::vector<int> &v);
  void removeDuplicates(std::vector<std::vector<int>> &v);
  void removeDuplicates(std::vector<std::vector<std::string>> &v);
};

#endif  // KFPARTICLESPHENIX_KFPARTICLETOOLS_H
