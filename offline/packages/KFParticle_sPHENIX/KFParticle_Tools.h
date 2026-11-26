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

class KFParticle_eventView;
class Vertex;
class SvtxTrack;
class SvtxTrack_v5;

class KFParticle_Tools
{
 public:
  /// Used by KFParticle_sPHENIX. The config is owned there and is not yet constructed when this runs:
  /// store the reference only, never read through it here.
  explicit KFParticle_Tools(const KFParticle_Config &config);

  virtual ~KFParticle_Tools() = default;

  KFParticle makeParticle(const SvtxTrack *track) const;

  std::vector<KFParticle> makeAllDaughterParticles(const KFParticle_eventView &eventView) const;

  bool isGoodTrack(const KFParticle_eventView &eventView, const KFParticle &particle, const std::vector<KFParticle> &primaryVertices) const;

  bool isTPConlyAtZeroCrossing(const KFParticle_eventView &eventView, const KFParticle &track) const;

  int calcMinPV_DCA(const KFParticle_eventView &eventView, const KFParticle &track, const std::vector<KFParticle> &PVs, float &minimumPV_DCA, float &minimumPV_DCA_stddev, bool do3D = true) const;

  std::vector<std::vector<int>> findTwoProngs(const KFParticle_eventView &eventView, const std::vector<KFParticle> &daughterParticles, std::vector<int> goodTrackIndex, int nTracks, const std::vector<KFParticle> &primaryVertices) const;

  std::vector<std::vector<int>> findNProngs(const KFParticle_eventView &eventView, const std::vector<KFParticle> &daughterParticles,
                                            const std::vector<int> &goodTrackIndex,
                                            std::vector<std::vector<int>> goodTracksThatMeet,
                                            int nRequiredTracks, unsigned int nProngs, const std::vector<KFParticle> &primaryVertices) const;

  std::vector<std::vector<int>> appendTracksToIntermediates(const KFParticle_eventView &eventView, KFParticle intermediateResonances[], const std::vector<KFParticle> &daughterParticles, const std::vector<int> &goodTrackIndex, int num_remaining_tracks, const std::vector<KFParticle> &primaryVertices) const;

  std::tuple<KFParticle, bool> buildMother(KFParticle vDaughters[], int daughterOrder[], bool isIntermediate, int intermediateNumber, int nTracks, bool constrainMass, float required_vertexID, KFParticle_eventView &eventView) const;

  void constrainToVertex(KFParticle &particle, bool &goodCandidate, KFParticle &vertex) const;

  std::tuple<KFParticle, bool> getCombination(KFParticle vDaughters[], int daughterOrder[], KFParticle vertex,
                                              bool constrain_to_vertex, bool isIntermediate, int intermediateNumber, int nTracks, bool constrainMass, float required_vertexID, KFParticle_eventView &eventView) const;

  void init_dEdx_fits();

  double get_dEdx_fitValue(float momentum, int PID) const;

  bool checkTrackAndVertexMatch(const KFParticle_eventView &eventView, KFParticle vDaughters[], int nTracks, const KFParticle &vertex) const;

  std::vector<int> getParticleCrossings(const KFParticle_eventView &eventView, const KFParticle &particle) const;

 protected:
  TF1 *f_pion_plus{nullptr};
  TF1 *f_kaon_plus{nullptr};
  TF1 *f_proton_plus{nullptr};
  TF1 *f_pion_minus{nullptr};
  TF1 *f_kaon_minus{nullptr};
  TF1 *f_proton_minus{nullptr};

  std::map<int, TF1 *> pidMap;

 private:
  const KFParticle_Config &m_config;
};

#endif  // KFPARTICLESPHENIX_KFPARTICLETOOLS_H
