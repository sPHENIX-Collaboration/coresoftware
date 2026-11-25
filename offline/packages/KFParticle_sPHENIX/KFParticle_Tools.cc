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

#include "KFParticle_Tools.h"
#include "KFParticle_util.h"

#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrack_v5.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/TrackAnalysisUtils.h>

#include <g4detectors/PHG4TpcGeom.h>
#include <g4detectors/PHG4TpcGeomContainer.h>

#include <trackbase/ActsGeometry.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrDefs.h>

#include <trackbase_historic/SvtxTrackSeed_v2.h>
#include <trackbase_historic/TrackSeedContainer_v1.h>
#include <trackbase_historic/TrackSeed_v2.h>

#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMap.h>
#include <globalvertex/SvtxVertex.h>
#include <globalvertex/SvtxVertexMap.h>

#include <phool/getClass.h>

#include <ffamodules/CDBInterface.h>

// KFParticle stuff
#include <KFParticle.h>
#include <KFVertex.h>

#include <Rtypes.h>
#include <TFile.h>

#include <algorithm>  // for max, remove, minmax_el...
#include <cmath>      // for sqrt, pow, M_PI
#include <cstdlib>    // for abs, NULL
#include <iostream>   // for operator<<, basic_ostream
#include <iterator>   // for end
#include <limits>
#include <map>        // for _Rb_tree_iterator, map
#include <memory>     // for allocator_traits<>::va...
#include <numeric>    // for iota

namespace
{
  // Shared default configuration for standalone KFParticle_Tools objects. A function-local static, so it
  // exists before first use even when the object using it is itself a global in another file.
  const KFParticle_Config &defaultConfig()
  {
    static const KFParticle_Config config;
    return config;
  }
using BunchCrossingMode = KFParticle_Config::BunchCrossingMode;

  // True if the two positions agree to within 1e-6 cm in every coordinate
  bool vertexToleranceCheck(float vertexMap_pos[3], float kfp_vertex_pos[3])
  {
    constexpr float epsilon = 1e-6F;
    for (unsigned int i = 0; i < 3; ++i)
    {
      if (std::abs(vertexMap_pos[i] - kfp_vertex_pos[i]) > epsilon)
      {
        return false;
      }
    }
    return true;
  }

  // The MBD measures only z. Its vertices have no x or y (MbdVertexv1 returns NaN, v2 and v3 return 0), no covariance
  // apart from the z error, and no chi2 or NDF, so they are read differently from SVTX vertices.
  bool isMbdVertex(const Vertex *vertex)
  {
    return dynamic_cast<const MbdVertex *>(vertex) != nullptr;
  }
}  // namespace

/// KFParticle constructor
KFParticle_Tools::KFParticle_Tools()
  : KFParticle_Tools(defaultConfig())
{
}

KFParticle_Tools::KFParticle_Tools(const KFParticle_Config &config)
  : m_config(config)
{
}

KFParticle KFParticle_Tools::makeVertex(const Vertex *vertex) const
{
  const bool isMbd = isMbdVertex(vertex);

  float vtxX = isMbd ? 0 : vertex->get_x();
  float vtxY = isMbd ? 0 : vertex->get_y();
  float vtxZ = vertex->get_z();

  float f_vertexParameters[6] = {vtxX, vtxY, vtxZ, 0, 0, 0};

  float f_vertexCovariance[21] = {0};
  unsigned int iterate = 0;
  if (isMbd)
  {
    f_vertexCovariance[5] = vertex->get_z_err();
  }
  else
  {
    for (unsigned int i = 0; i < 3; ++i)
    {
      for (unsigned int j = 0; j <= i; ++j)
      {
        f_vertexCovariance[iterate] = vertex->get_error(i, j);
        ++iterate;
      }
    }
  }

  KFParticle kfp_vertex;
  kfp_vertex.Create(f_vertexParameters, f_vertexCovariance, 0, -1);
  kfp_vertex.NDF() = isMbd ? 0 : vertex->get_ndof();
  kfp_vertex.Chi2() = isMbd ? 0 : vertex->get_chisq();

  return kfp_vertex;
}

const Vertex *KFParticle_Tools::findDstVertex(unsigned int id) const
{
  if (m_config.useMbdVertex)
  {
    return m_dst_mbdvertexmap->get(id);
  }
  return m_dst_vertexmap->get(id);
}

std::vector<KFParticle> KFParticle_Tools::makeAllPrimaryVertices(PHCompositeNode *topNode, const std::string &vertexMapName)
{
  const std::string &vtxMN = vertexMapName.empty() ? m_config.vertexMapNodeName : vertexMapName;

  std::vector<KFParticle> primaryVertices;

  if (m_config.useMbdVertex)
  {
    m_dst_mbdvertexmap = findNode::getClass<MbdVertexMap>(topNode, "MbdVertexMap");
  }
  else
  {
    m_dst_vertexmap = findNode::getClass<SvtxVertexMap>(topNode, vtxMN);
  }

  if (m_config.dontUseGlobalVertex)
  {
    if (m_config.useMbdVertex)
    {
      for (MbdVertexMap::ConstIter iter = m_dst_mbdvertexmap->begin(); iter != m_dst_mbdvertexmap->end(); ++iter)
      {
        primaryVertices.push_back(makeVertex(iter->second));
        primaryVertices.back().SetId(iter->first);
      }
    }
    else
    {
      for (SvtxVertexMap::ConstIter iter = m_dst_vertexmap->begin(); iter != m_dst_vertexmap->end(); ++iter)
      {
        primaryVertices.push_back(makeVertex(iter->second));
        primaryVertices.back().SetId(iter->first);
      }
    }

    return primaryVertices;
  }

  m_dst_globalvertexmap = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");
  if (!m_dst_globalvertexmap)
  {
    std::cout << "Can't continue in KFParticle_Tools::makeAllPrimaryVertices" << std::endl;
    return primaryVertices;
  }

  for (GlobalVertexMap::ConstIter iter = m_dst_globalvertexmap->begin(); iter != m_dst_globalvertexmap->end(); ++iter)
  {
    const GlobalVertex *globalVertex = iter->second;

    GlobalVertex::VTXTYPE whichVtx = m_config.useMbdVertex ? GlobalVertex::MBD : GlobalVertex::SVTX;

    auto svtxiter = globalVertex->find_vertexes(whichVtx);
    // check that it contains a track vertex
    if (svtxiter == globalVertex->end_vertexes())
    {
      continue;
    }

    const auto &svtxvertexvector = svtxiter->second;

    for (const auto &vertex : svtxvertexvector)
    {
      primaryVertices.push_back(makeVertex(findDstVertex(vertex->get_id())));
      primaryVertices.back().SetId(globalVertex->get_id());
    }
  }

  return primaryVertices;
}

KFParticle KFParticle_Tools::makeParticle(const SvtxTrack *track) const  /// Return a KFPTrack from track vector and covariance matrix. No mass or vertex constraints
{
  KFParticle kfp_particle;

  float f_trackParameters[6] = {track->get_x(),
                                track->get_y(),
                                track->get_z(),
                                track->get_px(),
                                track->get_py(),
                                track->get_pz()};

  float f_trackCovariance[21];
  unsigned int iterate = 0;
  for (unsigned int i = 0; i < 6; ++i)
  {
    for (unsigned int j = 0; j <= i; ++j)
    {
      f_trackCovariance[iterate] = track->get_error(i, j);
      ++iterate;
    }
  }

  kfp_particle.Create(f_trackParameters, f_trackCovariance, (Int_t) track->get_charge(), -1);
  kfp_particle.NDF() = track->get_ndf();
  kfp_particle.Chi2() = track->get_chisq();
  kfp_particle.SetId(track->get_id());

  return kfp_particle;
}

std::vector<KFParticle> KFParticle_Tools::makeAllDaughterParticles(PHCompositeNode *topNode)
{
  std::vector<KFParticle> daughterParticles;
  m_dst_trackmap = findNode::getClass<SvtxTrackMap>(topNode, m_config.trackMapNodeName);

  m_intermediate_crossings.clear();
  m_next_intermediate_id = -2;

  for (auto &iter : *m_dst_trackmap)
  {
    const SvtxTrack *track = iter.second;

    if (m_config.bunchCrossingZeroOnly && (track->get_crossing() != 0))
    {
      continue;
    }

    // Now lets check if we have the right number of tracker states
    int MVTX_states = 0;
    int INTT_states = 0;
    int TPC_states = 0;
    int TPOT_states = 0;

    // The per-detector cluster count getters only exist in SvtxTrack_v5 (not in the SvtxTrack base class)
    const auto *v5_track = dynamic_cast<const SvtxTrack_v5 *>(track);

    if (!v5_track)  // v1-v4: count the per-cluster track states
    {
      for (auto state_iter = track->begin_states();
           state_iter != track->end_states();
           ++state_iter)
      {
        SvtxTrackState *tstate = state_iter->second;
        if (tstate->get_pathlength() != 0)  // The first track state is an extrapolation so has no cluster
        {
          auto stateckey = tstate->get_cluskey();
          uint8_t id = TrkrDefs::getTrkrId(stateckey);

          switch (id)
          {
          case TrkrDefs::mvtxId:
            ++MVTX_states;
            break;
          case TrkrDefs::inttId:
            ++INTT_states;
            break;
          case TrkrDefs::tpcId:
            ++TPC_states;
            break;
          case TrkrDefs::micromegasId:
            ++TPOT_states;
            break;
          default:
            // std::cout << "Cluster key doesnt match a tracking system, could be related with projected track state to calorimeter system" << std::endl;
            break;
          }
        }
      }
    }
    else
    {
      MVTX_states = v5_track->get_nmvtx_clusters();
      INTT_states = v5_track->get_nintt_clusters();
      TPC_states = v5_track->get_ntpc_clusters();
      TPOT_states = v5_track->get_ntpot_clusters();
    }

    if (m_config.verbosity > 100)
    {
      KFParticle_util::printSelectionCheck("MVTX states", m_config.minMVTXStates, MVTX_states, 5);
      KFParticle_util::printSelectionCheck("INTT states", m_config.minINTTStates, INTT_states, 5);
      KFParticle_util::printSelectionCheck("TPC states", m_config.minTPCStates, TPC_states, 100);
      KFParticle_util::printSelectionCheck("TPOT states", m_config.minTPOTStates, TPOT_states, 5);
    }

    if (MVTX_states < m_config.minMVTXStates && !m_config.useTPCOnlyTracks)
    {
      continue;
    }
    if (INTT_states < m_config.minINTTStates && !m_config.useTPCOnlyTracks)
    {
      continue;
    }
    if (TPC_states < m_config.minTPCStates)
    {
      continue;
    }
    if (TPOT_states < m_config.minTPOTStates && !m_config.useTPCOnlyTracks)
    {
      continue;
    }
    if (m_config.useTPCOnlyTracks && (MVTX_states > 0 || INTT_states > 0))
    {
      continue;
    }

    const int trackndof = track->get_ndf();
    const float trackchi2ndof = trackndof > 0 ? track->get_chisq() / trackndof : std::numeric_limits<float>::quiet_NaN();
    if (m_config.verbosity > 11)
    {
      KFParticle_util::printSelectionCheck("Track chi^2/nDoF", 0, trackchi2ndof, m_config.trackMaxChi2nDOF);
    }
    if (!(trackchi2ndof <= m_config.trackMaxChi2nDOF)) //rejects nan
    {
      continue;
    }

    daughterParticles.push_back(makeParticle(track));  /// Turn all dst tracks in KFP tracks
    daughterParticles.back().SetId(iter.first);
  }

  return daughterParticles;
}

/*const*/ bool KFParticle_Tools::isGoodTrack(const KFParticle &particle, const std::vector<KFParticle> &primaryVertices)
{
  float min_PV_dca = 0;
  float min_PV_dca_stddev = 0;
  float min_PV_dca_xy = 0;
  float min_PV_dca_stddev_xy = 0;

  float pt = 0;
  float pterr = 0;
  particle.GetPt(pt, pterr);  // Both pt and pterr are passed by reference. GetPt() unhelpfully returns 0 for simulated silicon-only tracks, so we are changing to this function GetPt(p_t, sigma p_T)

  float ptchi2 = pow(pterr / pt, 2);
  calcMinPV_DCA(particle, primaryVertices, min_PV_dca, min_PV_dca_stddev);
  calcMinPV_DCA(particle, primaryVertices, min_PV_dca_xy, min_PV_dca_stddev_xy, false);

  bool goodTrack = KFParticle_util::isInRange(m_config.trackMinPt, pt, m_config.trackMaxPt) &&
                   ptchi2 <= m_config.trackMaxPtChi2 &&
                   min_PV_dca >= m_config.trackMinPVDCA &&
                   min_PV_dca_stddev >= m_config.trackMinPVDCAStdDev &&
                   min_PV_dca_xy >= m_config.trackMinPVDCAXY &&
                   min_PV_dca_stddev_xy >= m_config.trackMinPVDCAStdDevXY;

  if (m_config.verbosity >= 10)
  {
    KFParticle_util::printSelectionCheck("This track", "passed", "failed", "the selection", goodTrack);
    if (m_config.verbosity >= 11)
    {
      KFParticle_util::printSelectionCheck("Track pT", m_config.trackMinPt, pt, m_config.trackMaxPt);
      KFParticle_util::printSelectionCheck("Track pT chi^2", 0, ptchi2, m_config.trackMaxPtChi2);
      KFParticle_util::printSelectionCheck("PV DCA", m_config.trackMinPVDCA, min_PV_dca, std::numeric_limits<float>::max());
      KFParticle_util::printSelectionCheck("PV DCA Std. Dev.", m_config.trackMinPVDCAStdDev, min_PV_dca_stddev, std::numeric_limits<float>::max());
      KFParticle_util::printSelectionCheck("PV DCA xy", m_config.trackMinPVDCAXY, min_PV_dca_xy, std::numeric_limits<float>::max());
      KFParticle_util::printSelectionCheck("PV DCA xy Std. Dev.", m_config.trackMinPVDCAStdDevXY, min_PV_dca_stddev_xy, std::numeric_limits<float>::max());
    }
  }

  return goodTrack;
}

bool KFParticle_Tools::isTPConlyAtZeroCrossing(const KFParticle& track)
{
  const auto* svtxTrack = KFParticle_util::getTrack(track.Id(), m_dst_trackmap);

  if (!svtxTrack || svtxTrack->get_crossing() != 0)
  {
    return false;
  }

  // SvtxTrack_v5 keeps no seed pointers: TPC-only means no silicon clusters on the track
  if (const auto* v5 = dynamic_cast<const SvtxTrack_v5*>(svtxTrack))
  {
    return v5->get_nmvtx_clusters() == 0 && v5->get_nintt_clusters() == 0;
  }

  // SvtxTrack_v4 and earlier carry the seed pointer
  return svtxTrack->get_silicon_seed() == nullptr;
}

int KFParticle_Tools::calcMinPV_DCA(const KFParticle &track, const std::vector<KFParticle> &PVs,
                                float &minimumPV_DCA, float &minimumPV_DCA_stddev, bool do3D)
{
  std::vector<float> ip;
  std::vector<float> ip_significance;

  for (const auto &PV : PVs)
  {
    KFParticle trackArray[1] = {track};
    bool fromSameCrossing = checkTrackAndVertexMatch(trackArray, 1, PV);
    if (!isTPConlyAtZeroCrossing(track))
    {
      if (!fromSameCrossing)
      {
        continue;
      }
    }

    float thisPV_DCA_stddev = 0;

    if (do3D)
    {
      ip.push_back(track.GetDistanceFromVertex(PV));
      thisPV_DCA_stddev = track.GetDeviationFromVertex(PV);
    }
    else
    {
      ip.push_back(std::abs(track.GetDistanceFromVertexXY(PV)));
      thisPV_DCA_stddev = track.GetDeviationFromVertexXY(PV);
    }
    
    thisPV_DCA_stddev = std::max(thisPV_DCA_stddev, 0.F); // There are times where the PV_DCA_stddev calc fails
    ip_significance.push_back(thisPV_DCA_stddev);
  }

  if (ip.empty()) //Need to account for instances where the track has no associated primary vertex in its crossing. Track should always be rejected
  {
    minimumPV_DCA = std::numeric_limits<float>::quiet_NaN();
    minimumPV_DCA_stddev = std::numeric_limits<float>::quiet_NaN();
    return 0;
  }

  auto minmax_PV_dca = minmax_element(ip.begin(), ip.end());  // Order the PV_DCA from small to large
  minimumPV_DCA = *minmax_PV_dca.first;
  auto minmax_PV_dca_stddev = minmax_element(ip_significance.begin(), ip_significance.end());  // Order the PV_DCA chi2 from small to large
  minimumPV_DCA_stddev = *minmax_PV_dca_stddev.first;

  return 0;
}

std::vector<int> KFParticle_Tools::getParticleCrossings(const KFParticle &particle)
{
  std::vector<int> crossings;

  SvtxTrack *thisTrack = KFParticle_util::getTrack(particle.Id(), m_dst_trackmap);
  if (thisTrack)
  {
    crossings.push_back(thisTrack->get_crossing());
    return crossings;
  }

  auto it = m_intermediate_crossings.find(particle.Id());
  if (it != m_intermediate_crossings.end())
  {
    crossings = it->second;
  }

  return crossings;
}

std::vector<std::vector<int>> KFParticle_Tools::findTwoProngs(const std::vector<KFParticle> &daughterParticles, std::vector<int> goodTrackIndex, int nTracks, const std::vector<KFParticle> &primaryVertices)
{
  std::vector<std::vector<int>> goodTracksThatMeet;

  // The crossings of each particle do not depend on the pair, so get them once instead of for every pair
  const bool checkCrossings = m_config.bunchCrossingMode != BunchCrossingMode::Any;
  std::vector<std::vector<int>> particleCrossings;
  if (checkCrossings)
  {
    particleCrossings.resize(daughterParticles.size());
    for (const auto &idx : goodTrackIndex)
    {
      particleCrossings[idx] = getParticleCrossings(daughterParticles[idx]);
    }
  }
  std::vector<int> crossings;  // reused for every pair
  crossings.reserve(16);

  for (std::vector<int>::iterator i_it = goodTrackIndex.begin(); i_it != goodTrackIndex.end(); ++i_it)
  {
    for (std::vector<int>::iterator j_it = goodTrackIndex.begin(); j_it != goodTrackIndex.end(); ++j_it)
    {
      if (i_it < j_it)
      {
        if (checkCrossings)
        {
          crossings.clear();
          crossings.insert(crossings.end(), particleCrossings[*i_it].begin(), particleCrossings[*i_it].end());
          crossings.insert(crossings.end(), particleCrossings[*j_it].begin(), particleCrossings[*j_it].end());

          KFParticle_util::removeDuplicates(crossings);

          const bool sameCrossing = crossings.size() == 1;
          if ((m_config.bunchCrossingMode == BunchCrossingMode::Same && !sameCrossing) || (m_config.bunchCrossingMode == BunchCrossingMode::Mixed && sameCrossing))
          {
            continue;
          }
        }

        std::vector<KFParticle> dummy_tracks = {daughterParticles[*i_it], daughterParticles[*j_it]};

        KFParticle dummy_mother;
        dummy_mother.SetConstructMethod(2);

        for (auto &track : dummy_tracks)
        {
          dummy_mother.AddDaughter(track);
        }
        for (auto &track : dummy_tracks)
        {
          track.SetProductionVertex(dummy_mother);
        }

        float dca = dummy_tracks[0].GetDistanceFromParticle(dummy_tracks[1]);
        float dca_xy = std::abs(dummy_tracks[0].GetDistanceFromParticleXY(dummy_tracks[1]));

        if (m_config.verbosity >= 10)
        {
          KFParticle_util::printSelectionCheck("This track pair", "passed", "failed", "the DCA selection", (dca <= m_config.maxDaughterDCA) && (dca_xy <= m_config.maxDaughterDCAXY));
          if (m_config.verbosity >= 11)
          {
            KFParticle_util::printSelectionCheck("Pair DCA", 0., dca, m_config.maxDaughterDCA);
            KFParticle_util::printSelectionCheck("Pair DCA xy", 0., dca_xy, m_config.maxDaughterDCAXY);
          }
        }

        if (dca <= m_config.maxDaughterDCA && dca_xy <= m_config.maxDaughterDCAXY)
        {
          KFVertex twoParticleVertex;
          twoParticleVertex += dummy_tracks[0];
          twoParticleVertex += dummy_tracks[1];
          float vertexchi2ndof = twoParticleVertex.GetChi2() / twoParticleVertex.GetNDF();
          float sv_radial_position = sqrt(pow(twoParticleVertex.GetX(), 2) + pow(twoParticleVertex.GetY(), 2));
          std::vector<int> combination = {*i_it, *j_it};

          if (nTracks == 2 && m_config.verbosity >= 10)
          {
            KFParticle_util::printSelectionCheck("This track pair", "passed", "failed", "the quality and radius selection", (vertexchi2ndof <= m_config.maxVertexChi2nDOF) && (sv_radial_position >= m_config.minRadialSV));
            if (m_config.verbosity >= 11)
            {
              KFParticle_util::printSelectionCheck("SV chi^2/nDoF", 0., vertexchi2ndof, m_config.maxVertexChi2nDOF);
              KFParticle_util::printSelectionCheck("SV radius", m_config.minRadialSV, sv_radial_position, std::numeric_limits<float>::max());
            }
          }

          //Now check if tracks are good as we need full reco to make DCA calc make sense
          if (nTracks == 2)
          {
            if (vertexchi2ndof > m_config.maxVertexChi2nDOF || sv_radial_position < m_config.minRadialSV)
            {
              continue;
            }
       
            const bool allTracksPass = std::all_of(dummy_tracks.begin(), dummy_tracks.end(),
                                                   [&](const KFParticle &track) { return isGoodTrack(track, primaryVertices); });
            if (!allTracksPass)
            {
              continue;
            }

          }

          goodTracksThatMeet.push_back(combination);
        }
      }
    }
  }
  return goodTracksThatMeet;
}

std::vector<std::vector<int>> KFParticle_Tools::findNProngs(const std::vector<KFParticle> &daughterParticles,
                                                            const std::vector<int> &goodTrackIndex,
                                                            std::vector<std::vector<int>> goodTracksThatMeet,
                                                            int nRequiredTracks, unsigned int nProngs, const std::vector<KFParticle> &primaryVertices)
{
  unsigned int nGoodProngs = goodTracksThatMeet.size();

  std::vector<int> maxProngIndex;
  if (m_config.useOrderedProngCombinations)
  {
    maxProngIndex.reserve(nGoodProngs);
    for (const auto &prong : goodTracksThatMeet)
    {
      maxProngIndex.push_back(*std::max_element(prong.begin(), prong.end()));
    }
  }

  for (const auto &i_it : goodTrackIndex)
  {
    for (unsigned int i_prongs = 0; i_prongs < nGoodProngs; ++i_prongs)
    {
      if (m_config.useOrderedProngCombinations && i_it <= maxProngIndex[i_prongs])
      {
        continue;
      }

      const auto &prong = goodTracksThatMeet[i_prongs];
      const bool trackNotUsedAlready = std::find(prong.begin(), prong.end(), i_it) == prong.end();
      if (trackNotUsedAlready)
      {
        bool dcaMet = true;

        //Need to propagate all tracks first
        std::vector<int> combination;
        combination.reserve(nProngs);
        combination.push_back(i_it);
        combination.insert(combination.end(), prong.begin(), prong.end());

        KFParticle dummy_mother;
        std::vector<KFParticle> dummy_tracks;
        dummy_tracks.reserve(combination.size());
        for (auto &id : combination)
        {
          dummy_tracks.push_back(daughterParticles[id]);
        }
        dummy_mother.SetConstructMethod(2);

        for (auto &track : dummy_tracks)
        {
          dummy_mother.AddDaughter(track);
        }
        for (auto &track : dummy_tracks)
        {
          track.SetProductionVertex(dummy_mother);
        }

        for (unsigned int i = 1; i < combination.size(); ++i)
        {
          float dca = dummy_tracks[0].GetDistanceFromParticle(dummy_tracks[i]);
          float dca_xy = dummy_tracks[0].GetDistanceFromParticleXY(dummy_tracks[i]);

         if (m_config.verbosity >= 10)
         {
           KFParticle_util::printSelectionCheck("This track", "combined", "did not combine", "with a SV set", (dca <= m_config.maxDaughterDCA) && (dca_xy <= m_config.maxDaughterDCAXY));
           if (m_config.verbosity >= 11)
           {
             KFParticle_util::printSelectionCheck("Pair DCA", 0., dca, m_config.maxDaughterDCA);
             KFParticle_util::printSelectionCheck("Pair DCA xy", 0., dca_xy, m_config.maxDaughterDCAXY);
           }
         }

          if (dca > m_config.maxDaughterDCA || dca_xy > m_config.maxDaughterDCAXY)
          {
            dcaMet = false;
            break;  // this combination is rejected, no need to test the remaining tracks
          }
        }

        if (dcaMet)
        {
          // The vertex of the raw tracks is only needed for the quality and radius selection, so build it here
          KFVertex particleVertex;
          for (const int id : combination)
          {
            particleVertex += daughterParticles[id];
          }

          float vertexchi2ndof = particleVertex.GetChi2() / particleVertex.GetNDF();
          float sv_radial_position = sqrt(pow(particleVertex.GetX(), 2) + pow(particleVertex.GetY(), 2));

          if ((unsigned int) nRequiredTracks == nProngs && m_config.verbosity >= 10)
          {
            KFParticle_util::printSelectionCheck("This SV combination", "passed", "failed", "the quality and radius selection", (vertexchi2ndof <= m_config.maxVertexChi2nDOF) && (sv_radial_position >= m_config.minRadialSV));
            if (m_config.verbosity >= 11)
            {
              KFParticle_util::printSelectionCheck("SV chi^2/nDoF", 0., vertexchi2ndof, m_config.maxVertexChi2nDOF);
              KFParticle_util::printSelectionCheck("SV radius", m_config.minRadialSV, sv_radial_position, std::numeric_limits<float>::max());
            }
          }

          if ((unsigned int) nRequiredTracks == nProngs)
          {
            if (vertexchi2ndof > m_config.maxVertexChi2nDOF || sv_radial_position < m_config.minRadialSV)
            {
              continue;
            }

            const bool allTracksPass = std::all_of(dummy_tracks.begin(), dummy_tracks.end(),
                                                   [&](const KFParticle &track) { return isGoodTrack(track, primaryVertices); });
            if (!allTracksPass)
            {
              continue;
            }
          }

          goodTracksThatMeet.push_back(combination);
        }
      }
    }
  }

  goodTracksThatMeet.erase(goodTracksThatMeet.begin(), goodTracksThatMeet.begin() + nGoodProngs);
  for (auto &i : goodTracksThatMeet)
  {
    sort(i.begin(), i.end());
  }
  KFParticle_util::removeDuplicates(goodTracksThatMeet);

  return goodTracksThatMeet;
}

std::vector<std::vector<int>> KFParticle_Tools::appendTracksToIntermediates(KFParticle intermediateResonances[], const std::vector<KFParticle> &daughterParticles, const std::vector<int> &goodTrackIndex, int num_remaining_tracks, const std::vector<KFParticle> &primaryVertices)
{
  std::vector<std::vector<int>> goodTracksThatMeetIntermediates;

  // Do these tracks form a vertex together with the intermediates?
  auto combinesWithIntermediates = [&](const std::vector<int> &trackIds)
  {
    std::vector<KFParticle> v_intermediateResonances(intermediateResonances, intermediateResonances + m_config.numIntermediateStates);
    for (const int j : trackIds)
    {
      v_intermediateResonances.push_back(daughterParticles[j]);
    }
    std::vector<int> dummyTrackID(v_intermediateResonances.size());
    std::iota(dummyTrackID.begin(), dummyTrackID.end(), 0);
    std::vector<std::vector<int>> dummyTrackList = findTwoProngs(v_intermediateResonances, dummyTrackID, (int) v_intermediateResonances.size(), primaryVertices);
    for (unsigned int p = 3; p <= v_intermediateResonances.size(); ++p)
    {
      dummyTrackList = findNProngs(v_intermediateResonances, dummyTrackID, dummyTrackList, (int) v_intermediateResonances.size(), (int) p, primaryVertices);
    }
    return !dummyTrackList.empty();
  };

  if (num_remaining_tracks == 1)
  {
    for (const auto &i_it : goodTrackIndex)
    {
      if (combinesWithIntermediates({i_it}))
      {
        goodTracksThatMeetIntermediates.push_back({i_it});
      }
    }
  }
  else
  {
    std::vector<std::vector<int>> goodTracksThatMeet = findTwoProngs(daughterParticles, goodTrackIndex, num_remaining_tracks, primaryVertices);
    for (int p = 3; p <= num_remaining_tracks; ++p)
    {
      goodTracksThatMeet = findNProngs(daughterParticles, goodTrackIndex, goodTracksThatMeet, num_remaining_tracks, p, primaryVertices);
    }

    for (const auto &i : goodTracksThatMeet)
    {
      if (combinesWithIntermediates(i))
      {
        goodTracksThatMeetIntermediates.push_back(i);
      }
    }
  }

  return goodTracksThatMeetIntermediates;
}

std::tuple<KFParticle, bool> KFParticle_Tools::buildMother(KFParticle vDaughters[], int daughterOrder[],
                                                           bool isIntermediate, int intermediateNumber, int nTracks,
                                                           bool constrainMass, float required_vertexID, PHCompositeNode *topNode)
{
  KFParticle mother;
  std::vector<KFParticle> inputTracks(nTracks);

  mother.SetConstructMethod(2);

  bool daughterMassCheck = true;
  int particlesWithPID[] = {11, 211, 321, 2212};
  float unique_vertexID = 0;

  // Figure out if the decay has reco. tracks mixed with resonances
  int num_tracks_used_by_intermediates = 0;
  for (int i = 0; i < m_config.numIntermediateStates; ++i)
  {
    num_tracks_used_by_intermediates += m_config.numTracksFromIntermediate[i];
  }
  int num_remaining_tracks = m_config.numTracks - num_tracks_used_by_intermediates;

  for (int i = 0; i < nTracks; ++i)
  {
    float daughterMass = 0;

    if ((Int_t) vDaughters[i].GetQ() != 0)
    {
      // For charged particle, like p+/-, pi+/-, Sigma+/-, etc...
      // different charged particle has different PDGID
      // just to protect if they have different mass for different charge
      // but in EvtGen, there is no C-violation...so this is just a protection
      daughterMass = constrainMass ? KFParticle_util::getParticleMass((Int_t) vDaughters[i].GetQ() * daughterOrder[i]) : vDaughters[i].GetMass();
    }
    else if ((Int_t) vDaughters[i].GetQ() == 0)
    {
      // For neutral particle, like pi0, eta, J/psi, etc... who do not have an anti-particle with anti-PDGID
      // and other neutral particle, like Lambda0/anti-Lambda0 ... who have an anti-particle with anti-PDGID
      // avoid charge*PDGID=0 case and getting wrong mass
      daughterMass = constrainMass ? KFParticle_util::getParticleMass(daughterOrder[i]) : vDaughters[i].GetMass();
    }

    if ((num_remaining_tracks > 0 && i >= m_config.numIntermediateStates) || isIntermediate)
    {
      if ((Int_t) vDaughters[i].GetQ() != 0)
      {
        daughterMass = KFParticle_util::getParticleMass((Int_t) vDaughters[i].GetQ() * daughterOrder[i]);
      }
      else if ((Int_t) vDaughters[i].GetQ() == 0)
      {
        daughterMass = KFParticle_util::getParticleMass(daughterOrder[i]);
      }
    }
    inputTracks[i].Create(vDaughters[i].Parameters(),
                          vDaughters[i].CovarianceMatrix(),
                          (Int_t) vDaughters[i].GetQ(),
                          daughterMass);

    // Run PID check
    if (m_config.usePID)
    {
      int track_PDG_ID = (Int_t) vDaughters[i].GetQ() * daughterOrder[i];
      if (std::find(std::begin(particlesWithPID), std::end(particlesWithPID), std::abs(track_PDG_ID)) != std::end(particlesWithPID))
      {
        float calculated_dEdx_value = get_dEdx(topNode, vDaughters[i]);
        double expected_dEdx_value = get_dEdx_fitValue((Int_t) vDaughters[i].GetQ() * vDaughters[i].GetP(), track_PDG_ID);
        bool accept_dEdx = KFParticle_util::isInRange((1 - m_config.dEdxBandWidth) * expected_dEdx_value, calculated_dEdx_value, (1 + m_config.dEdxBandWidth) * expected_dEdx_value);

        if (m_config.verbosity >= 11)
        {
          KFParticle_util::printSelectionCheck("dE/dx check", (1 - m_config.dEdxBandWidth) * expected_dEdx_value, calculated_dEdx_value, (1 + m_config.dEdxBandWidth) * expected_dEdx_value);
        }

        if (!accept_dEdx)
        {
          return std::make_tuple(mother, false);
        }
      }
    }

    mother.AddDaughter(inputTracks[i]);
    mother.AddDaughterId(vDaughters[i].Id());
    unique_vertexID += (Int_t) vDaughters[i].GetQ() * KFParticle_util::getParticleMass(daughterOrder[i]);
  }

  if (isIntermediate)
  {
    mother.SetPDG(KFParticle_util::getParticleID(m_config.intermediateNames[intermediateNumber]));
  }
  if (!isIntermediate && !m_config.motherName.empty())
  {
    mother.SetPDG(KFParticle_util::getParticleID(m_config.motherName));
  }

  const bool chargeCheck = m_config.chargeConjugate ? std::abs(unique_vertexID) == std::abs(required_vertexID)
                                                  : unique_vertexID == required_vertexID;

  for (int j = 0; j < nTracks; ++j)
  {
    if (m_config.extrapolateTracksToSV)
    {
      inputTracks[j].SetProductionVertex(mother);
    }
    if (!m_config.allowZeroMassTracks && inputTracks[j].GetMass() == 0)
    {
      daughterMassCheck = false;
    }
  }

  float calculated_mass;
  float calculated_mass_err;
  mother.GetMass(calculated_mass, calculated_mass_err);
  float calculated_pt;
  float calculated_pt_err;
  mother.GetPt(calculated_pt, calculated_pt_err);

  float min_mass = isIntermediate ? m_config.intermediateMassRange[intermediateNumber].first : m_config.motherMinMass;
  float max_mass = isIntermediate ? m_config.intermediateMassRange[intermediateNumber].second : m_config.motherMaxMass;
  float min_pt = isIntermediate ? m_config.intermediateMinPt[intermediateNumber] : m_config.motherMinPt;

  float max_vertex_volume = isIntermediate ? m_config.intermediateMaxVertexVolume[intermediateNumber] : m_config.motherMaxVertexVolume;

  bool goodCandidate = KFParticle_util::isInRange(min_mass, calculated_mass, max_mass) &&
                       calculated_pt >= min_pt &&
                       daughterMassCheck &&
                       chargeCheck &&
                       KFParticle_util::calculateEllipsoidVolume(mother) <= max_vertex_volume;

  std::vector<int> crossings;
  for (int i = 0; i < nTracks; ++i)
  {
    std::vector<int> daughterCrossings = getParticleCrossings(vDaughters[i]);
    crossings.insert(crossings.end(), daughterCrossings.begin(), daughterCrossings.end());
  }
  KFParticle_util::removeDuplicates(crossings);

  if (goodCandidate && m_config.bunchCrossingMode != BunchCrossingMode::Any)
  {
    const bool sameCrossing = crossings.size() == 1;
    if ((m_config.bunchCrossingMode == BunchCrossingMode::Same && !sameCrossing) || (m_config.bunchCrossingMode == BunchCrossingMode::Mixed && sameCrossing))
    {
      goodCandidate = false;
    }

    if (m_config.verbosity >= 11)
    {
      KFParticle_util::printSelectionCheck("", "All tracks are from the same BC", "Tracks are from different BC", "", sameCrossing);
    }
  }

  if (isIntermediate && goodCandidate)
  {
    mother.SetId(m_next_intermediate_id);
    m_intermediate_crossings[m_next_intermediate_id] = crossings;
    --m_next_intermediate_id;
  }

  // Check the requirements of an intermediate states against this mother and re-do goodCandidate
  if (goodCandidate && m_config.hasIntermediates && !isIntermediate)  // The decay has intermediate states and we are now looking at the mother
  {
    for (int k = 0; k < m_config.numIntermediateStates; ++k)
    {
      float intermediate_DIRA = KFParticle_util::eventDIRA(vDaughters[k], mother);
      float intermediate_FDchi2 = KFParticle_util::flightDistanceChi2(vDaughters[k], mother);
      if (intermediate_DIRA < m_config.intermediateMinDIRA[k] ||
          intermediate_FDchi2 < m_config.intermediateMinFDChi2[k])
      {
        goodCandidate = false;
      }

      // Decay length and time of the intermediate state from the decay vertex of the mother, calculated in the same
      // way as for the mother in constrainToVertex. The intermediate is refitted to the mother vertex on a temporary
      // copy (SetProductionVertex changes the particle), so neither the intermediate nor the mother is altered
      const bool cutOnIntermediateDecay = k < (int) m_config.intermediateDecayLengthRange.size() || k < (int) m_config.intermediateDecayLengthRangeXY.size() || k < (int) m_config.intermediateDecayTimeRange.size() || k < (int) m_config.intermediateDecayTimeRangeXY.size();
      if (goodCandidate && cutOnIntermediateDecay)
      {
        KFParticle intermediateCopy = vDaughters[k];
        intermediateCopy.SetProductionVertex(mother);
        intermediateCopy.TransportToDecayVertex();

        float intermediate_decayTime;
        float intermediate_decayTimeErr;
        float intermediate_decayLength;
        float intermediate_decayLengthErr;
        float intermediate_decayLength_xy;
        float intermediate_decayLengthErr_xy;

        intermediateCopy.GetLifeTime(intermediate_decayTime, intermediate_decayTimeErr);
        intermediateCopy.GetDecayLength(intermediate_decayLength, intermediate_decayLengthErr);
        intermediateCopy.GetDecayLengthXY(intermediate_decayLength_xy, intermediate_decayLengthErr_xy);
        float intermediate_decayTime_xy = intermediateCopy.GetPseudoProperDecayTime(mother, intermediateCopy.GetMass());

        intermediate_decayTime /= KFParticle_util::speedOfLight;
        intermediate_decayTime_xy /= KFParticle_util::speedOfLight;  // As for the mother, the pseudo-proper decay time is returned as a length (c*tau)

        const std::pair<float, float> noCut(-std::numeric_limits<float>::max(), std::numeric_limits<float>::max());
        auto rangeOf = [&noCut, k](const std::vector<std::pair<float, float>> &ranges)
        { return k < (int) ranges.size() ? ranges[k] : noCut; };
        const auto lengthRange = rangeOf(m_config.intermediateDecayLengthRange);
        const auto lengthRange_xy = rangeOf(m_config.intermediateDecayLengthRangeXY);
        const auto timeRange = rangeOf(m_config.intermediateDecayTimeRange);
        const auto timeRange_xy = rangeOf(m_config.intermediateDecayTimeRangeXY);

        goodCandidate = KFParticle_util::isInRange(lengthRange.first, intermediate_decayLength, lengthRange.second) &&
                        KFParticle_util::isInRange(lengthRange_xy.first, intermediate_decayLength_xy, lengthRange_xy.second) &&
                        KFParticle_util::isInRange(timeRange.first, intermediate_decayTime, timeRange.second) &&
                        KFParticle_util::isInRange(timeRange_xy.first, intermediate_decayTime_xy, timeRange_xy.second);

        if (m_config.verbosity >= 11)
        {
          KFParticle_util::printSelectionCheck("Intermediate Decay Length", lengthRange.first, intermediate_decayLength, lengthRange.second);
          KFParticle_util::printSelectionCheck("Intermediate Decay Length xy", lengthRange_xy.first, intermediate_decayLength_xy, lengthRange_xy.second);
          KFParticle_util::printSelectionCheck("Intermediate Decay Time", timeRange.first, intermediate_decayTime, timeRange.second);
          KFParticle_util::printSelectionCheck("Intermediate Decay Time xy", timeRange_xy.first, intermediate_decayTime_xy, timeRange_xy.second);
        }
      }

      if (m_config.verbosity >= 10)
      {
        KFParticle_util::printSelectionCheck("", "Accepted", "Rejected", "the intermediate selection", goodCandidate);
        if (m_config.verbosity >= 11)
        {
          KFParticle_util::printSelectionCheck("Intermediate DIRA", m_config.intermediateMinDIRA[k], intermediate_DIRA, std::numeric_limits<float>::max());
          KFParticle_util::printSelectionCheck("Intermediate FD chi^2", m_config.intermediateMinFDChi2[k], intermediate_FDchi2, std::numeric_limits<float>::max());
        }
      }
    }
  }

  if (m_config.verbosity >= 10)
  {
    KFParticle_util::printSelectionCheck("", "Accepted", "Rejected", "the mother selection", goodCandidate);
    if (m_config.verbosity >= 11)
    {
      KFParticle_util::printSelectionCheck("Vertex charge is", "right", "wrong", "", chargeCheck);
      KFParticle_util::printSelectionCheck("Invariant Mass", min_mass, calculated_mass, max_mass);
      KFParticle_util::printSelectionCheck("Mother pT", min_pt, calculated_pt, std::numeric_limits<float>::max());
      KFParticle_util::printSelectionCheck("Mother SV volume", 0., KFParticle_util::calculateEllipsoidVolume(mother), max_vertex_volume);
    }
  }
  return std::make_tuple(mother, goodCandidate);
}

void KFParticle_Tools::constrainToVertex(KFParticle &particle, bool &goodCandidate, KFParticle &vertex) const
{
  KFParticle particleCopy = particle;
  particleCopy.SetProductionVertex(vertex);
  particleCopy.TransportToDecayVertex();

  float calculated_decayTime;
  float calculated_decayTimeErr;
  float calculated_decayLength;
  float calculated_decayLengthErr;
  float calculated_decayLength_xy;
  float calculated_decayLengthErr_xy;

  particleCopy.GetLifeTime(calculated_decayTime, calculated_decayTimeErr);
  particleCopy.GetDecayLength(calculated_decayLength, calculated_decayLengthErr);
  particleCopy.GetDecayLengthXY(calculated_decayLength_xy, calculated_decayLengthErr_xy);

  float calculated_decayTime_xy = particleCopy.GetPseudoProperDecayTime(vertex, particleCopy.GetMass());

  float calculated_fdchi2 = KFParticle_util::flightDistanceChi2(particle, vertex);

  float calculated_PV_dca_xy = std::abs(particle.GetDistanceFromVertexXY(vertex));
  float calculated_PV_dca_stddev_xy = particle.GetDeviationFromVertexXY(vertex);
  float calculated_dira_xy = KFParticle_util::eventDIRA(particle, vertex, false);
  float calculated_PV_dca = particle.GetDistanceFromVertex(vertex);
  float calculated_PV_dca_stddev = particle.GetDeviationFromVertex(vertex);
  float calculated_dira = KFParticle_util::eventDIRA(particle, vertex);

  float calculated_decay_time_significance = calculated_decayTime / calculated_decayTimeErr;
  float calculated_decay_length_significance = calculated_decayLength / calculated_decayLengthErr;
  float calculated_decay_length_xy_significance = calculated_decayLength_xy / calculated_decayLengthErr_xy;

  calculated_decayTime /= KFParticle_util::speedOfLight;
  calculated_decayTime_xy /= KFParticle_util::speedOfLight;  // The pseudo-proper decay time is returned as a length (c*tau), convert it to a time like the 3D one

  goodCandidate = calculated_fdchi2 >= m_config.motherMinFDChi2 &&
                  calculated_PV_dca <= m_config.motherMaxPVDCA &&
                  calculated_PV_dca_stddev <= m_config.motherMaxPVDCAStdDev &&
                  calculated_PV_dca_xy <= m_config.motherMaxPVDCAXY &&
                  calculated_PV_dca_stddev_xy <= m_config.motherMaxPVDCAStdDevXY &&
                  calculated_decay_time_significance >= m_config.motherMinDecayTimeSignificance &&
                  calculated_decay_length_significance >= m_config.motherMinDecayLengthSignificance &&
                  calculated_decay_length_xy_significance >= m_config.motherMinDecayLengthSignificanceXY &&
                  KFParticle_util::isInRange(m_config.motherMinDIRA, calculated_dira, m_config.motherMaxDIRA) &&
                  KFParticle_util::isInRange(m_config.motherMinDIRAXY, calculated_dira_xy, m_config.motherMaxDIRAXY) &&
                  KFParticle_util::isInRange(m_config.motherMinDecayTime, calculated_decayTime, m_config.motherMaxDecayTime) &&
                  KFParticle_util::isInRange(m_config.motherMinDecayTimeXY, calculated_decayTime_xy, m_config.motherMaxDecayTimeXY) &&
                  KFParticle_util::isInRange(m_config.motherMinDecayLength, calculated_decayLength, m_config.motherMaxDecayLength) &&
                  KFParticle_util::isInRange(m_config.motherMinDecayLengthXY, calculated_decayLength_xy, m_config.motherMaxDecayLengthXY);

  if (m_config.verbosity >= 10)
  {
    KFParticle_util::printSelectionCheck("", "Passed", "Failed", "the PV constraint", goodCandidate);
    if (m_config.verbosity >= 11)
    {
      KFParticle_util::printSelectionCheck("Mother DIRA", m_config.motherMinDIRA, calculated_dira, m_config.motherMaxDIRA);
      KFParticle_util::printSelectionCheck("Mother DIRA xy", m_config.motherMinDIRAXY, calculated_dira_xy, m_config.motherMaxDIRAXY);
      KFParticle_util::printSelectionCheck("Mother FD chi^2", m_config.motherMinFDChi2, calculated_fdchi2, std::numeric_limits<float>::max());
      KFParticle_util::printSelectionCheck("Mother PV DCA", 0, calculated_PV_dca, m_config.motherMaxPVDCA);
      KFParticle_util::printSelectionCheck("Mother PV DCA Std. Dev.", 0., calculated_PV_dca_stddev, m_config.motherMaxPVDCAStdDev);
      KFParticle_util::printSelectionCheck("Mother PV DCA xy", 0., calculated_PV_dca_xy, m_config.motherMaxPVDCAXY);
      KFParticle_util::printSelectionCheck("Mother PV DCA xy Std. Dev.", 0., calculated_PV_dca_stddev_xy, m_config.motherMaxPVDCAStdDevXY);
      KFParticle_util::printSelectionCheck("Mother Decay Time", m_config.motherMinDecayTime, calculated_decayTime, m_config.motherMaxDecayTime);
      KFParticle_util::printSelectionCheck("Mother Decay Time Significance", m_config.motherMinDecayTimeSignificance, calculated_decay_time_significance, std::numeric_limits<float>::max());
      KFParticle_util::printSelectionCheck("Mother Decay Time xy", m_config.motherMinDecayTimeXY, calculated_decayTime_xy, m_config.motherMaxDecayTimeXY);
      KFParticle_util::printSelectionCheck("Mother Decay Length", m_config.motherMinDecayLength, calculated_decayLength, m_config.motherMaxDecayLength);
      KFParticle_util::printSelectionCheck("Mother Decay Length Significance", m_config.motherMinDecayLengthSignificance, calculated_decay_length_significance, std::numeric_limits<float>::max());
      KFParticle_util::printSelectionCheck("Mother Decay Length xy", m_config.motherMinDecayLengthXY, calculated_decayLength_xy, m_config.motherMaxDecayLengthXY);
      KFParticle_util::printSelectionCheck("Mother Decay Length xy Significance", m_config.motherMinDecayLengthSignificanceXY, calculated_decay_length_xy_significance, std::numeric_limits<float>::max());
    }
  }
}

std::tuple<KFParticle, bool> KFParticle_Tools::getCombination(KFParticle vDaughters[], int daughterOrder[], KFParticle vertex, bool constrain_to_vertex, bool isIntermediate, int intermediateNumber, int nTracks, bool constrainMass, float required_vertexID, PHCompositeNode *topNode)
{
  auto [candidate, isGoodCandidate] = buildMother(vDaughters, daughterOrder, isIntermediate, intermediateNumber, nTracks, constrainMass, required_vertexID, topNode);

  if (constrain_to_vertex && isGoodCandidate && !isIntermediate)
  {
    if (m_config.requireTrackVertexBunchCrossingMatch && m_config.bunchCrossingMode != BunchCrossingMode::Mixed)
    {
      isGoodCandidate = checkTrackAndVertexMatch(vDaughters, nTracks, vertex);
    }

    if (isGoodCandidate)
    {
      constrainToVertex(candidate, isGoodCandidate, vertex);
    }
  }
  return std::make_tuple(candidate, isGoodCandidate);
}

float KFParticle_Tools::get_dEdx(PHCompositeNode *topNode, const KFParticle &daughter)
{
  SvtxTrackMap *trackmap = findNode::getClass<SvtxTrackMap>(topNode, m_config.trackMapNodeName);
  TrkrClusterContainer *clusterMap = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  if (!clusterMap)
  {
    clusterMap = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER_SEED");
  }
  auto *geomContainer = findNode::getClass<PHG4TpcGeomContainer>(topNode, "TPCGEOMCONTAINER");
  auto *geometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  if (!clusterMap || !geomContainer || !geometry)
  {
    return -1.0;
  }

  SvtxTrack *daughter_track = KFParticle_util::getTrack(daughter.Id(), trackmap);
  if (!daughter_track)
  {
    // No backing SvtxTrack -- e.g. daughter is itself a built intermediate.
    // dE/dx PID is only meaningful for a real reconstructed track.
    return -1.0;
  }

  bool isTrack_v5 = daughter_track->InheritsFrom("SvtxTrack_v5");
  TrackSeedContainer *tpcSeedContainer = isTrack_v5 ? findNode::getClass<TrackSeedContainer>(topNode, "TpcTrackSeedContainer") : nullptr;

  TrackSeed *tpcseed = isTrack_v5 ? KFParticle_util::findSeed(tpcSeedContainer, daughter_track->get_tpc_seed_id()) : daughter_track->get_tpc_seed();
  if (!tpcseed)
  {
    return -1;
  }

  float layerThicknesses[4] = {0.0, 0.0, 0.0, 0.0};
  // These are randomly chosen layer thicknesses for the TPC, to get the
  // correct region thicknesses in an easy to pass way to the helper fxn
  layerThicknesses[0] = geomContainer->GetLayerCellGeom(7)->get_thickness();
  layerThicknesses[1] = geomContainer->GetLayerCellGeom(8)->get_thickness();
  layerThicknesses[2] = geomContainer->GetLayerCellGeom(27)->get_thickness();
  layerThicknesses[3] = geomContainer->GetLayerCellGeom(50)->get_thickness();

  return TrackAnalysisUtils::calc_dedx(tpcseed, clusterMap, geometry, layerThicknesses);
}

void KFParticle_Tools::init_dEdx_fits()
{
  const std::string dedx_fitparams = m_config.useLocalPIDFile ? m_config.localPIDFilename : CDBInterface::instance()->getUrl("TPC_DEDX_FITPARAM");

  std::cout << PHWHERE << " opening " << dedx_fitparams << std::endl;

  TFile *filefit = TFile::Open(dedx_fitparams.c_str());

  if (!filefit->IsOpen())
  {
    std::cerr << "Error opening filefit!" << std::endl;
    return;
  }

  if (m_config.useLocalPIDFile)
  {
    if (m_config.verbosity > 4)
    {
      std::cout << PHWHERE << " using local file " << m_config.localPIDFilename << std::endl;
    }
    // new method is independent of charge
    filefit->GetObject("pi_band",f_pion_plus);
    filefit->GetObject("K_band",f_kaon_plus);
    filefit->GetObject("p_band",f_proton_plus);
    filefit->GetObject("pi_band",f_pion_minus);
    filefit->GetObject("K_band",f_kaon_minus);
    filefit->GetObject("p_band",f_proton_minus);
  }
  else
  {
    filefit->GetObject("f_piband", f_pion_plus);
    filefit->GetObject("f_Kband", f_kaon_plus);
    filefit->GetObject("f_pband", f_proton_plus);
    filefit->GetObject("f_piminus_band", f_pion_minus);
    filefit->GetObject("f_Kminus_band", f_kaon_minus);
    filefit->GetObject("f_pbar_band", f_proton_minus);
  }

  pidMap.emplace(-11, f_pion_plus);
  pidMap.emplace(211, f_pion_plus);
  pidMap.emplace(321, f_kaon_plus);
  pidMap.emplace(2212, f_proton_plus);
  pidMap.emplace(11, f_pion_minus);
  pidMap.emplace(-211, f_pion_minus);
  pidMap.emplace(-321, f_kaon_minus);
  pidMap.emplace(-2212, f_proton_minus);
}

double KFParticle_Tools::get_dEdx_fitValue(float momentum, int PID)
{
  return pidMap[PID]->Eval(momentum);
}

bool KFParticle_Tools::checkTrackAndVertexMatch(KFParticle vDaughters[], int nTracks, const KFParticle &vertex)
{
  int vertexCrossing = std::numeric_limits<int>::min();

  if (m_config.useFakePrimaryVertex) //Fake PV has no crossing and should always match to a track 
  {
    return true;
  }

  float obtained_vertex[3] = {0, 0, std::numeric_limits<float>::lowest()};
  float kfp_vertex[3] = {vertex.GetX(), vertex.GetY(), vertex.GetZ()};

  // Takes the crossing of a DST vertex at the position of the KFParticle vertex (MBD: z only, x and y stay 0)
  auto matchCrossing = [&](const Vertex *dstVertex)
  {
    if (!isMbdVertex(dstVertex))
    {
      obtained_vertex[0] = dstVertex->get_x();
      obtained_vertex[1] = dstVertex->get_y();
    }
    obtained_vertex[2] = dstVertex->get_z();
    if (vertexToleranceCheck(obtained_vertex, kfp_vertex))
    {
      vertexCrossing = dstVertex->get_beam_crossing();
    }
  };

  if (m_config.dontUseGlobalVertex)
  {
    const Vertex *dstVertex = findDstVertex(vertex.Id());
    if (!dstVertex)
    {
      return false;
    }
    matchCrossing(dstVertex);
  }
  else
  {
    const GlobalVertex *globalVertex = m_dst_globalvertexmap->get(vertex.Id());
    if (!globalVertex)
    {
      return false;
    }

    GlobalVertex::VTXTYPE whichVtx = m_config.useMbdVertex ? GlobalVertex::MBD : GlobalVertex::SVTX;

    auto svtxiter = globalVertex->find_vertexes(whichVtx);
    if (svtxiter == globalVertex->end_vertexes())
    {
      return false;
    }
    const auto &svtxvertexvector = svtxiter->second;

    for (const auto &gvertex : svtxvertexvector)
    {
      const Vertex *dstVertex = findDstVertex(gvertex->get_id());
      if (!dstVertex)
      {
        continue;
      }
      matchCrossing(dstVertex);
    }
  }

  for (int i = 0; i < nTracks; ++i)
  {
    for (const int trackCrossing : getParticleCrossings(vDaughters[i]))
    {
      if (trackCrossing != vertexCrossing)
      {
        return false;  // no point checking remaining tracks
      }
    }
  }

  return true;
}
