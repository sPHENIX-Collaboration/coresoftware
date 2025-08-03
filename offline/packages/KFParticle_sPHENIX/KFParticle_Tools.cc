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

#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>

#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrDefs.h> 
#include <g4detectors/PHG4TpcCylinderGeom.h>
#include <g4detectors/PHG4TpcCylinderGeomContainer.h>

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
#include <TDatabasePDG.h>
#include <TMatrixD.h>
#include "KFParticle_truthAndDetTools.h"

#include <TFile.h>
#include <TMatrixDfwd.h>  // for TMatrixD
#include <TMatrixT.h>     // for TMatrixT, operator*

#include <Eigen/Dense>

#include <algorithm>  // for max, remove, minmax_el...
#include <cmath>      // for sqrt, pow, M_PI
#include <cstdlib>    // for abs, NULL
#include <iostream>   // for operator<<, basic_ostream
#include <iterator>   // for end
#include <map>        // for _Rb_tree_iterator, map
#include <memory>     // for allocator_traits<>::va...

KFParticle_truthAndDetTools toolSet;

/// KFParticle constructor
KFParticle_Tools::KFParticle_Tools()
  : m_has_intermediates(false)
  , m_min_mass(0)
  , m_max_mass(0)
  , m_min_decayTime(-1 * FLT_MAX)
  , m_max_decayTime(FLT_MAX)
  , m_min_decayLength(-1 * FLT_MAX)
  , m_max_decayLength(FLT_MAX)
  , m_track_pt(0.2)
  , m_track_ptchi2(FLT_MAX)
  , m_track_ip_xy(-1.)
  , m_track_ipchi2_xy(-1)
  , m_track_ip(-1.)
  , m_track_ipchi2(-1)
  , m_track_chi2ndof(4.)
  , m_nMVTXStates(3)
  , m_nINTTStates(1)
  , m_nTPCStates(20)
  , m_nTPOTStates(0)
  , m_comb_DCA_xy(0.05)
  , m_comb_DCA(0.05)
  , m_vertex_chi2ndof(15.)
  , m_fdchi2(0.)
  , m_dira_min(0.90)
  , m_dira_max(1.01)
  , m_mother_pt(0.)
  , m_mother_ipchi2(FLT_MAX)
  , m_get_charge_conjugate(false)
  , m_extrapolateTracksToSV(true)
  , m_vtx_map_node_name("SvtxVertexMap")
  , m_trk_map_node_name("SvtxTrackMap")
  , m_dst_mbdvertexmap()
  , m_dst_mbdvertex()
  , m_dst_trackmap()
  , m_dst_track()
  , m_dst_vertexmap()
  , m_dst_vertex()
  , m_cluster_map()
  , m_geom_container()
{
}

KFParticle KFParticle_Tools::makeVertex(PHCompositeNode * /*topNode*/)
{
  float vtxX = m_use_mbd_vertex ? 0 : m_dst_vertex->get_x();
  float vtxY = m_use_mbd_vertex ? 0 : m_dst_vertex->get_y();
  float vtxZ = m_use_mbd_vertex ? m_dst_mbdvertex->get_z() : m_dst_vertex->get_z();

  float f_vertexParameters[6] = {vtxX, vtxY, vtxZ, 0, 0, 0};

  float f_vertexCovariance[21] = {0};
  unsigned int iterate = 0;
  if (m_use_mbd_vertex)
  {
    f_vertexCovariance[5] = m_dst_mbdvertex->get_z_err();
  }
  else
  {
    for (unsigned int i = 0; i < 3; ++i)
    {
      for (unsigned int j = 0; j <= i; ++j)
      {
        f_vertexCovariance[iterate] = m_dst_vertex->get_error(i, j);
        ++iterate;
      }
    }
  }

  KFParticle kfp_vertex;
  kfp_vertex.Create(f_vertexParameters, f_vertexCovariance, 0, -1);
  kfp_vertex.NDF() = m_use_mbd_vertex ? 0 : m_dst_vertex->get_ndof();
  kfp_vertex.Chi2() = m_use_mbd_vertex ? 0 : m_dst_vertex->get_chisq();

  return kfp_vertex;
}

std::vector<KFParticle> KFParticle_Tools::makeAllPrimaryVertices(PHCompositeNode *topNode, const std::string &vertexMapName)
{
  std::string vtxMN;

  unsigned int vertexID = 0;

  if (vertexMapName.empty())
  {
    vtxMN = m_vtx_map_node_name;
  }
  else
  {
    vtxMN = vertexMapName;
  }

  std::vector<KFParticle> primaryVertices;

  if (m_use_mbd_vertex)
  {
    m_dst_mbdvertexmap = findNode::getClass<MbdVertexMap>(topNode, "MbdVertexMap");
  }
  else
  {
    m_dst_vertexmap = findNode::getClass<SvtxVertexMap>(topNode, vtxMN);
  }

  if (m_dont_use_global_vertex)
  {
    if (m_use_mbd_vertex)
    {
      for (MbdVertexMap::ConstIter iter = m_dst_mbdvertexmap->begin(); iter != m_dst_mbdvertexmap->end(); ++iter)
      {
        m_dst_mbdvertex = iter->second;
        primaryVertices.push_back(makeVertex(topNode));
        primaryVertices[vertexID].SetId(iter->first);
        ++vertexID;
      }
    }
    else
    {
       for (SvtxVertexMap::ConstIter iter = m_dst_vertexmap->begin(); iter != m_dst_vertexmap->end(); ++iter)
      {
        m_dst_vertex = iter->second;
        primaryVertices.push_back(makeVertex(topNode));
        primaryVertices[vertexID].SetId(iter->first);
        ++vertexID;
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
    m_dst_globalvertex = iter->second;
    
    GlobalVertex::VTXTYPE whichVtx = m_use_mbd_vertex ? GlobalVertex::MBD : GlobalVertex::SVTX; 

    auto svtxiter = m_dst_globalvertex->find_vertexes(whichVtx);
    // check that it contains a track vertex
    if (svtxiter == m_dst_globalvertex->end_vertexes())
    {
      continue;
    }

    auto svtxvertexvector = svtxiter->second;

    for (auto &vertex : svtxvertexvector)
    {
      if (m_use_mbd_vertex)
      {
        m_dst_mbdvertex = m_dst_mbdvertexmap->find(vertex->get_id())->second;
      }
      else
      {
        m_dst_vertex = m_dst_vertexmap->find(vertex->get_id())->second;
      }

      primaryVertices.push_back(makeVertex(topNode));
      primaryVertices[vertexID].SetId(m_dst_globalvertex->get_id());
      ++vertexID;
    }
  }

  return primaryVertices;
}

KFParticle KFParticle_Tools::makeParticle(PHCompositeNode * /*topNode*/)  /// Return a KFPTrack from track vector and covariance matrix. No mass or vertex constraints
{
  KFParticle kfp_particle;

  float f_trackParameters[6] = {m_dst_track->get_x(),
                                m_dst_track->get_y(),
                                m_dst_track->get_z(),
                                m_dst_track->get_px(),
                                m_dst_track->get_py(),
                                m_dst_track->get_pz()};

  float f_trackCovariance[21];
  unsigned int iterate = 0;
  for (unsigned int i = 0; i < 6; ++i)
  {
    for (unsigned int j = 0; j <= i; ++j)
    {
      f_trackCovariance[iterate] = m_dst_track->get_error(i, j);
      ++iterate;
    }
  }

  kfp_particle.Create(f_trackParameters, f_trackCovariance, (Int_t) m_dst_track->get_charge(), -1);
  kfp_particle.NDF() = m_dst_track->get_ndf();
  kfp_particle.Chi2() = m_dst_track->get_chisq();
  kfp_particle.SetId(m_dst_track->get_id());

  return kfp_particle;
}

std::vector<KFParticle> KFParticle_Tools::makeAllDaughterParticles(PHCompositeNode *topNode)
{
  std::vector<KFParticle> daughterParticles;
  m_dst_trackmap = findNode::getClass<SvtxTrackMap>(topNode, m_trk_map_node_name.c_str());
  unsigned int trackID = 0;

  for (auto &iter : *m_dst_trackmap)
  {
    m_dst_track = iter.second;

    if (m_bunch_crossing_zero_only && (m_dst_track->get_crossing() != 0))
    {
      continue;
    }

    //Now lets check if we have the right number of tracker states
    int MVTX_states = 0;
    int INTT_states = 0;
    int TPC_states = 0;
    int TPOT_states = 0;

    for (auto state_iter = m_dst_track->begin_states();
         state_iter != m_dst_track->end_states();
         ++state_iter)
    {
      SvtxTrackState* tstate = state_iter->second;
      if (tstate->get_pathlength() != 0) //The first track state is an extrapolation so has no cluster
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
            //std::cout << "Cluster key doesnt match a tracking system, could be related with projected track state to calorimeter system" << std::endl;
            break; 
        }
      }
    }

    if (MVTX_states < m_nMVTXStates)
    {
      continue;
    }
    if (INTT_states < m_nINTTStates)
    {
      continue;
    }
    if (TPC_states < m_nTPCStates)
    {
      continue;
    }
    if (TPOT_states < m_nTPOTStates)
    {
      continue;
    }
 
    daughterParticles.push_back(makeParticle(topNode));  /// Turn all dst tracks in KFP tracks
    daughterParticles[trackID].SetId(iter.first);
    ++trackID;
  }

  return daughterParticles;
}

void KFParticle_Tools::getTracksFromBC(PHCompositeNode *topNode, const int &bunch_crossing, const std::string &vertexMapName, int &nTracks, int &nPVs)
{
  if (m_use_mbd_vertex) //If you're using the MBD vertex then there is no way to know which tracks are associated to it
  {
    return;
  }

  std::string vtxMN;
  if (vertexMapName.empty())
  {
    vtxMN = m_vtx_map_node_name;
  }
  else
  {
    vtxMN = vertexMapName;
  }

  nTracks = 0;
  nPVs = 0;
  m_dst_vertexmap = findNode::getClass<SvtxVertexMap>(topNode, vtxMN);
  for (SvtxVertexMap::ConstIter iter = m_dst_vertexmap->begin(); iter != m_dst_vertexmap->end(); ++iter)
  {
    m_dst_vertex = iter->second;
    if ((int) m_dst_vertex->get_beam_crossing() == bunch_crossing)
    {
      nTracks += m_dst_vertex->size_tracks();
      ++nPVs;
    }
  }
}

int KFParticle_Tools::getTracksFromVertex(PHCompositeNode *topNode, const KFParticle &vertex, const std::string &vertexMapName)
{
  if (m_use_mbd_vertex) //If you're using the MBD vertex then there is no way to know which tracks are associated to it
  {
    return 0;
  }

  std::string vtxMN;
  if (vertexMapName.empty())
  {
    vtxMN = m_vtx_map_node_name;
  }
  else
  {
    vtxMN = vertexMapName;
  }

  SvtxVertex* associated_vertex = nullptr;
  if (m_dont_use_global_vertex)
  {
    m_dst_vertexmap = findNode::getClass<SvtxVertexMap>(topNode, vtxMN);
    associated_vertex = m_dst_vertexmap->get(vertex.Id());
  }
  else
  {
    m_dst_globalvertexmap = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");
    if (!m_dst_globalvertexmap)
    {
      std::cout << "Can't continue in KFParticle_Tools::makeAllPrimaryVertices" << std::endl;
      return 0;
    }
    auto associated_gvertex = m_dst_globalvertexmap->get(vertex.Id());

    auto svtxiter = associated_gvertex->find_vertexes(GlobalVertex::SVTX);
    auto svtxvertexvector = svtxiter->second;

    for (auto &gvertex : svtxvertexvector)
    {
      associated_vertex = m_dst_vertexmap->find(gvertex->get_id())->second;
    }
  }
  if (associated_vertex)
  {
    return associated_vertex->size_tracks();
  }
  else
  {
    return 0;
  }
}

/*const*/ bool KFParticle_Tools::isGoodTrack(const KFParticle &particle, const std::vector<KFParticle> &primaryVertices)
{
  bool goodTrack = false;

  float min_ip = 0;
  float min_ipchi2 = 0;
  float min_ip_xy = 0;
  float min_ipchi2_xy = 0;

  float pt = particle.GetPt();
  float pterr = particle.GetErrPt();
  float ptchi2 = pow(pterr / pt, 2);
  float trackchi2ndof = particle.GetChi2() / particle.GetNDF();
  calcMinIP(particle, primaryVertices, min_ip, min_ipchi2);
  calcMinIP(particle, primaryVertices, min_ip_xy, min_ipchi2_xy, false);

  if (pt >= m_track_pt
   && ptchi2 <= m_track_ptchi2
   && min_ip >= m_track_ip
   && min_ipchi2 >= m_track_ipchi2 
   && min_ip_xy >= m_track_ip_xy
   && min_ipchi2_xy >= m_track_ipchi2_xy 
   && trackchi2ndof <= m_track_chi2ndof)
  {
    goodTrack = true;
  }
  return goodTrack;
}

int KFParticle_Tools::calcMinIP(const KFParticle &track, const std::vector<KFParticle> &PVs,
                                float &minimumIP, float &minimumIPchi2, bool do3D)
{
  std::vector<float> ip, ipchi2;

  for (const auto &PV : PVs)
  {
    float thisIPchi2 = 0;

    if (do3D) 
    {
      ip.push_back(track.GetDistanceFromVertex(PV));
      track.GetDeviationFromVertex(PV);
    }
    else
    {
      ip.push_back(abs(track.GetDistanceFromVertexXY(PV)));
      track.GetDeviationFromVertexXY(PV);
    }

    if (thisIPchi2 < 0)
    {
      thisIPchi2 = 0;
    }
    ipchi2.push_back(thisIPchi2);  //Τhere are times where the IPchi2 calc fails
  }

  auto minmax_ip = minmax_element(ip.begin(), ip.end());  // Order the IP from small to large
  minimumIP = *minmax_ip.first;
  auto minmax_ipchi2 = minmax_element(ipchi2.begin(), ipchi2.end());  // Order the IP chi2 from small to large
  minimumIPchi2 = *minmax_ipchi2.first;

  return 0;
}

std::vector<int> KFParticle_Tools::findAllGoodTracks(const std::vector<KFParticle> &daughterParticles, const std::vector<KFParticle> &primaryVertices)
{
  std::vector<int> goodTrackIndex;

  for (unsigned int i_parts = 0; i_parts < daughterParticles.size(); ++i_parts)
  {
    if (isGoodTrack(daughterParticles[i_parts], primaryVertices))
    {
      goodTrackIndex.push_back(i_parts);
    }
  }

  removeDuplicates(goodTrackIndex);

  return goodTrackIndex;
}

std::vector<std::vector<int>> KFParticle_Tools::findTwoProngs(std::vector<KFParticle> daughterParticles, std::vector<int> goodTrackIndex, int nTracks)
{
  std::vector<std::vector<int>> goodTracksThatMeet;

  for (std::vector<int>::iterator i_it = goodTrackIndex.begin(); i_it != goodTrackIndex.end(); ++i_it)
  {
    for (std::vector<int>::iterator j_it = goodTrackIndex.begin(); j_it != goodTrackIndex.end(); ++j_it)
    {
      if (i_it < j_it)
      {
        float dca = daughterParticles[*i_it].GetDistanceFromParticle(daughterParticles[*j_it]);
        float dca_xy = abs(daughterParticles[*i_it].GetDistanceFromParticleXY(daughterParticles[*j_it]));

        if (dca <= m_comb_DCA && dca_xy <= m_comb_DCA_xy)
        {
          KFVertex twoParticleVertex;
          twoParticleVertex += daughterParticles[*i_it];
          twoParticleVertex += daughterParticles[*j_it];
          float vertexchi2ndof = twoParticleVertex.GetChi2() / twoParticleVertex.GetNDF();
          float sv_radial_position = sqrt(pow(twoParticleVertex.GetX(), 2) + pow(twoParticleVertex.GetY(), 2));
          std::vector<int> combination = {*i_it, *j_it};

          if (nTracks == 2 && vertexchi2ndof > m_vertex_chi2ndof)
          {
            continue;
          }
          else
          {
            if (nTracks == 2 && sv_radial_position < m_min_radial_SV)
            {
              continue;
            }
            else
            {
              goodTracksThatMeet.push_back(combination);
            }
          }
        }
      }
    }
  }

  return goodTracksThatMeet;
}

std::vector<std::vector<int>> KFParticle_Tools::findNProngs(std::vector<KFParticle> daughterParticles,
                                                            const std::vector<int> &goodTrackIndex,
                                                            std::vector<std::vector<int>> goodTracksThatMeet,
                                                            int nRequiredTracks, unsigned int nProngs)
{
  unsigned int nGoodProngs = goodTracksThatMeet.size();

  for (auto &i_it : goodTrackIndex)
  {
    for (unsigned int i_prongs = 0; i_prongs < nGoodProngs; ++i_prongs)
    {
      bool trackNotUsedAlready = true;
      for (unsigned int i_trackCheck = 0; i_trackCheck < nProngs - 1; ++i_trackCheck)
      {
        if (i_it == goodTracksThatMeet[i_prongs][i_trackCheck])
        {
          trackNotUsedAlready = false;
        }
      }
      if (trackNotUsedAlready)
      {
        bool dcaMet = true;
        for (unsigned int i = 0; i < nProngs - 1; ++i)
        {
          float dca = daughterParticles[i_it].GetDistanceFromParticle(daughterParticles[goodTracksThatMeet[i_prongs][i]]);
          float dca_xy = abs(daughterParticles[i_it].GetDistanceFromParticleXY(daughterParticles[goodTracksThatMeet[i_prongs][i]]));

          if (dca > m_comb_DCA || dca_xy > m_comb_DCA_xy)
          {
            dcaMet = false;
          }
        }

        if (dcaMet)
        {
          KFVertex particleVertex;
          particleVertex += daughterParticles[i_it];
          std::vector<int> combination;
          combination.push_back(i_it);
          for (unsigned int i = 0; i < nProngs - 1; ++i)
          {
            particleVertex += daughterParticles[goodTracksThatMeet[i_prongs][i]];
            combination.push_back(goodTracksThatMeet[i_prongs][i]);
          }
          float vertexchi2ndof = particleVertex.GetChi2() / particleVertex.GetNDF();
          float sv_radial_position = sqrt(pow(particleVertex.GetX(), 2) + pow(particleVertex.GetY(), 2));

          if ((unsigned int) nRequiredTracks == nProngs && vertexchi2ndof > m_vertex_chi2ndof)
          {
            continue;
          }
          else
          {
            if ((unsigned int) nRequiredTracks == nProngs && sv_radial_position < m_min_radial_SV)
            {
              continue;
            }
            else
            {
              goodTracksThatMeet.push_back(combination);
            }
          }
        }
      }
    }
  }

  goodTracksThatMeet.erase(goodTracksThatMeet.begin(), goodTracksThatMeet.begin() + nGoodProngs);
  for (auto &i : goodTracksThatMeet)
  {
    sort(i.begin(), i.end());
  }
  removeDuplicates(goodTracksThatMeet);

  return goodTracksThatMeet;
}

std::vector<std::vector<int>> KFParticle_Tools::appendTracksToIntermediates(KFParticle intermediateResonances[], const std::vector<KFParticle> &daughterParticles, const std::vector<int> &goodTrackIndex, int num_remaining_tracks)
{
  std::vector<std::vector<int>> goodTracksThatMeet, goodTracksThatMeetIntermediates;  //, vectorOfGoodTracks;
  if (num_remaining_tracks == 1)
  {
    for (auto &i_it : goodTrackIndex)
    {
      std::vector<KFParticle> v_intermediateResonances(intermediateResonances, intermediateResonances + m_num_intermediate_states);
      std::vector<std::vector<int>> dummyTrackList;
      std::vector<int> dummyTrackID;  // I already have the track ids stored in goodTracksThatMeet[i]
      v_intermediateResonances.insert(end(v_intermediateResonances), daughterParticles[i_it]);
      for (unsigned int k = 0; k < v_intermediateResonances.size(); ++k)
      {
        dummyTrackID.push_back(k);
      }
      dummyTrackList = findTwoProngs(v_intermediateResonances, dummyTrackID, (int) v_intermediateResonances.size());
      if (v_intermediateResonances.size() > 2)
      {
        for (unsigned int p = 3; p <= v_intermediateResonances.size(); ++p)
        {
          dummyTrackList = findNProngs(v_intermediateResonances,
                                       dummyTrackID, dummyTrackList,
                                       (int) v_intermediateResonances.size(), (int) p);
        }
      }

      if (dummyTrackList.size() != 0)
      {
        std::vector<int> goodTrack{i_it};
        goodTracksThatMeetIntermediates.push_back(goodTrack);
      }
    }
  }
  else
  {
    goodTracksThatMeet = findTwoProngs(daughterParticles, goodTrackIndex, num_remaining_tracks);

    for (int p = 3; p <= num_remaining_tracks; ++p)
    {
      goodTracksThatMeet = findNProngs(daughterParticles, goodTrackIndex, goodTracksThatMeet, num_remaining_tracks, p);
    }

    for (auto &i : goodTracksThatMeet)
    {
      std::vector<KFParticle> v_intermediateResonances(intermediateResonances, intermediateResonances + m_num_intermediate_states);
      std::vector<std::vector<int>> dummyTrackList;
      std::vector<int> dummyTrackID;  // I already have the track ids stored in goodTracksThatMeet[i]
      for (int j : i)
      {
        v_intermediateResonances.push_back(daughterParticles[i[j]]);
      }
      for (unsigned int k = 0; k < v_intermediateResonances.size(); ++k)
      {
        dummyTrackID.push_back(k);
      }
      dummyTrackList = findTwoProngs(v_intermediateResonances, dummyTrackID, (int) v_intermediateResonances.size());
      for (unsigned int p = 3; p <= v_intermediateResonances.size(); ++p)
      {
        dummyTrackList = findNProngs(v_intermediateResonances, dummyTrackID, dummyTrackList, (int) v_intermediateResonances.size(), (int) p);
      }

      if (dummyTrackList.size() != 0)
      {
        goodTracksThatMeetIntermediates.push_back(i);
      }
    }
  }

  return goodTracksThatMeetIntermediates;
}

float KFParticle_Tools::eventDIRA(const KFParticle &particle, const KFParticle &vertex, bool do3D)
{
  const int nDimensions = do3D ? 3 : 2;
  TMatrixD flightVector(nDimensions, 1);
  TMatrixD momVector(nDimensions, 1);
  flightVector(0, 0) = particle.GetX() - vertex.GetX();
  flightVector(1, 0) = particle.GetY() - vertex.GetY();

  momVector(0, 0) = particle.GetPx();
  momVector(1, 0) = particle.GetPy();
  
  if (do3D)
  {
    flightVector(2, 0) = particle.GetZ() - vertex.GetZ();
    momVector(2, 0) = particle.GetPz();
  }

  TMatrixD momDotFD(1, 1);  // Calculate momentum dot flight distance
  momDotFD = TMatrixD(momVector, TMatrixD::kTransposeMult, flightVector);
  float f_momDotFD = momDotFD(0, 0);

  TMatrixD sizeOfMom(1, 1);  // Calculates the size of the momentum vector
  sizeOfMom = TMatrixD(momVector, TMatrixD::kTransposeMult, momVector);
  float f_sizeOfMom = sqrt(sizeOfMom(0, 0));

  TMatrixD sizeOfFD(1, 1);  // Calculates the size of the flight distance vector
  sizeOfFD = TMatrixD(flightVector, TMatrixD::kTransposeMult, flightVector);
  float f_sizeOfFD = sqrt(sizeOfFD(0, 0));

  return f_momDotFD / (f_sizeOfMom * f_sizeOfFD);
}

float KFParticle_Tools::flightDistanceChi2(const KFParticle &particle, const KFParticle &vertex)
{
  TMatrixD flightVector(3, 1);
  TMatrixD flightDistanceCovariance(3, 3);

  const KFParticle &kfp_vertex(vertex);

  flightVector(0, 0) = particle.GetX() - kfp_vertex.GetX();
  flightVector(1, 0) = particle.GetY() - kfp_vertex.GetY();
  flightVector(2, 0) = particle.GetZ() - kfp_vertex.GetZ();

  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      flightDistanceCovariance(i, j) = particle.GetCovariance(i, j) + kfp_vertex.GetCovariance(i, j);
    }
  }

  TMatrixD anInverseMatrix(3, 3);
  anInverseMatrix = flightDistanceCovariance.Invert();
  TMatrixD m_chi2Value(1, 1);
  m_chi2Value = TMatrixD(flightVector, TMatrixD::kTransposeMult, anInverseMatrix * flightVector);
  return m_chi2Value(0, 0);
}

std::tuple<KFParticle, bool> KFParticle_Tools::buildMother(KFParticle vDaughters[], int daughterOrder[],
                                                           bool isIntermediate, int intermediateNumber, int nTracks,
                                                           bool constrainMass, float required_vertexID, PHCompositeNode* topNode)
{
  KFParticle mother;
  KFParticle *inputTracks = new KFParticle[nTracks];

  mother.SetConstructMethod(2);

  bool daughterMassCheck = true;
  int particlesWithPID[] = {11, 211, 321, 2212};
  float unique_vertexID = 0;

  // Figure out if the decay has reco. tracks mixed with resonances
  int num_tracks_used_by_intermediates = 0;
  for (int i = 0; i < m_num_intermediate_states; ++i)
  {
    num_tracks_used_by_intermediates += m_num_tracks_from_intermediate[i];
  }
  int num_remaining_tracks = m_num_tracks - num_tracks_used_by_intermediates;

  for (int i = 0; i < nTracks; ++i)
  {
    float daughterMass = 0;

    if ((Int_t) vDaughters[i].GetQ() != 0)
    {
      // For charged particle, like p+/-, pi+/-, Sigma+/-, etc...
      // different charged particle has different PDGID
      // just to protect if they have different mass for different charge
      // but in EvtGen, there is no C-violation...so this is just a protection
      daughterMass = constrainMass ? getParticleMass((Int_t) vDaughters[i].GetQ() * daughterOrder[i]) : vDaughters[i].GetMass();
    }
    else if ((Int_t) vDaughters[i].GetQ() == 0)
    {
      // For neutral particle, like pi0, eta, J/psi, etc... who do not have an anti-particle with anti-PDGID
      // and other neutral particle, like Lambda0/anti-Lambda0 ... who have an anti-particle with anti-PDGID
      // avoid charge*PDGID=0 case and getting wrong mass
      daughterMass = constrainMass ? getParticleMass(daughterOrder[i]) : vDaughters[i].GetMass();
    }

    if ((num_remaining_tracks > 0 && i >= m_num_intermediate_states) || isIntermediate)
    {
      if ((Int_t) vDaughters[i].GetQ() != 0)
      {
        daughterMass = getParticleMass((Int_t) vDaughters[i].GetQ() * daughterOrder[i]);
      }
      else if ((Int_t) vDaughters[i].GetQ() == 0)
      {
        daughterMass = getParticleMass(daughterOrder[i]);
      }

    }
    inputTracks[i].Create(vDaughters[i].Parameters(),
                          vDaughters[i].CovarianceMatrix(),
                          (Int_t) vDaughters[i].GetQ(),
                          daughterMass);

    //Run PID check
    if (m_use_PID)
    {
      int track_PDG_ID = (Int_t) vDaughters[i].GetQ()*daughterOrder[i];
      if (std::find(std::begin(particlesWithPID), std::end(particlesWithPID), std::abs(track_PDG_ID)) != std::end(particlesWithPID))
      {
        float calculated_dEdx_value = get_dEdx(topNode, vDaughters[i]);
        double expected_dEdx_value = get_dEdx_fitValue((Int_t) vDaughters[i].GetQ() * vDaughters[i].GetP(), track_PDG_ID);
        bool accept_dEdx = isInRange((1-m_dEdx_band_width)*expected_dEdx_value, calculated_dEdx_value, (1+m_dEdx_band_width)*expected_dEdx_value);
        if(!accept_dEdx)
        {
         delete [] inputTracks;
         return std::make_tuple(mother, false);
        }
      }
    }

    mother.AddDaughter(inputTracks[i]);
    mother.AddDaughterId(vDaughters[i].Id());
    unique_vertexID += (Int_t) vDaughters[i].GetQ() * getParticleMass(daughterOrder[i]);
  }

  if (isIntermediate)
  {
    mother.SetPDG(getParticleID(m_intermediate_name[intermediateNumber].c_str()));
  }
  if (!isIntermediate && !m_mother_name_Tools.empty())
  {
    mother.SetPDG(getParticleID(m_mother_name_Tools));
  }

  bool chargeCheck;
  if (m_get_charge_conjugate)
  {
    chargeCheck = std::abs(unique_vertexID) == std::abs(required_vertexID) ? true : false;
  }
  else
  {
    chargeCheck = unique_vertexID == required_vertexID ? true : false;
  }

  for (int j = 0; j < nTracks; ++j)
  {
    if (m_extrapolateTracksToSV)
    {
      inputTracks[j].SetProductionVertex(mother);
    }
    if (!m_allowZeroMassTracks)
    {
      if (inputTracks[j].GetMass() == 0)
      {
        daughterMassCheck = false;
      }
    }
  }


  float calculated_mass, calculated_mass_err;
  mother.GetMass(calculated_mass, calculated_mass_err);
  float calculated_pt = mother.GetPt();

  float min_mass = isIntermediate ? m_intermediate_mass_range[intermediateNumber].first : m_min_mass;
  float max_mass = isIntermediate ? m_intermediate_mass_range[intermediateNumber].second : m_max_mass;
  float min_pt = isIntermediate ? m_intermediate_min_pt[intermediateNumber] : m_mother_pt;

  float max_vertex_volume = isIntermediate ? m_intermediate_vertex_volume[intermediateNumber] : m_mother_vertex_volume;

  bool goodCandidate = false;

  if (calculated_mass >= min_mass && calculated_mass <= max_mass &&
      calculated_pt >= min_pt && daughterMassCheck && chargeCheck && calculateEllipsoidVolume(mother) <= max_vertex_volume)
  {
    goodCandidate = true;
  }

  if (goodCandidate && m_require_bunch_crossing_match)
  {
    std::vector<int> crossings;
    for (int i = 0; i < nTracks; ++i)
    {
      SvtxTrack *thisTrack = toolSet.getTrack(vDaughters[i].Id(), m_dst_trackmap);
      if (thisTrack)//This protects against intermediates which have no track but I need a way to assign the bunch crossing to an interemdiate as this was already checked when it was actually built
      {
        crossings.push_back(thisTrack->get_crossing());
      }
    }

    removeDuplicates(crossings);

    if (crossings.size() != 1)
    {
      goodCandidate = false;
    }
  }

  // Check the requirements of an intermediate states against this mother and re-do goodCandidate
  if (goodCandidate && m_has_intermediates && !isIntermediate)  // The decay has intermediate states and we are now looking at the mother
  {
    for (int k = 0; k < m_num_intermediate_states; ++k)
    {
      float intermediate_DIRA = eventDIRA(vDaughters[k], mother);
      float intermediate_FDchi2 = flightDistanceChi2(vDaughters[k], mother);
      if (intermediate_DIRA < m_intermediate_min_dira[k] ||
          intermediate_FDchi2 < m_intermediate_min_fdchi2[k])
      {
        goodCandidate = false;
      }
    }
  }
  delete [] inputTracks;
  return std::make_tuple(mother, goodCandidate);
}

void KFParticle_Tools::constrainToVertex(KFParticle &particle, bool &goodCandidate, KFParticle &vertex)
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

  float calculated_fdchi2 = flightDistanceChi2(particle, vertex);

  float calculated_ip_xy = abs(particle.GetDistanceFromVertexXY(vertex));
  float calculated_ipchi2_xy = particle.GetDeviationFromVertexXY(vertex);
  float calculated_dira_xy = eventDIRA(particle, vertex, false);
  float calculated_ip = particle.GetDistanceFromVertex(vertex);
  float calculated_ipchi2 = particle.GetDeviationFromVertex(vertex);
  float calculated_dira = eventDIRA(particle, vertex);

  float calculated_decay_time_significance = calculated_decayTime/calculated_decayTimeErr;
  float calculated_decay_length_significance = calculated_decayLength/calculated_decayLengthErr;
  float calculated_decay_length_xy_significance = calculated_decayLength_xy/calculated_decayLengthErr_xy;

  goodCandidate = false;

  const float speed = 2.99792458e-2;
  calculated_decayTime /= speed;

  if (calculated_fdchi2 >= m_fdchi2 
   && calculated_ip <= m_mother_ip 
   && calculated_ipchi2 <= m_mother_ipchi2 
   && calculated_ip_xy <= m_mother_ip_xy
   && calculated_ipchi2_xy <= m_mother_ipchi2_xy
   && calculated_decay_time_significance >= m_mother_min_decay_time_significance
   && calculated_decay_length_significance >= m_mother_min_decay_length_significance
   && calculated_decay_length_xy_significance >= m_mother_min_decay_length_xy_significance
   && isInRange(m_dira_min, calculated_dira, m_dira_max) 
   && isInRange(m_dira_xy_min, calculated_dira_xy, m_dira_xy_max) 
   && isInRange(m_min_decayTime, calculated_decayTime, m_max_decayTime) 
   && isInRange(m_min_decayTime_xy, calculated_decayTime_xy, m_max_decayTime_xy) 
   && isInRange(m_min_decayLength, calculated_decayLength, m_max_decayLength)
   && isInRange(m_min_decayLength_xy, calculated_decayLength_xy, m_max_decayLength_xy))
  {
    goodCandidate = true;
  }
}

std::tuple<KFParticle, bool> KFParticle_Tools::getCombination(KFParticle vDaughters[], int daughterOrder[], KFParticle vertex, bool constrain_to_vertex, bool isIntermediate, int intermediateNumber, int nTracks, bool constrainMass, float required_vertexID, PHCompositeNode* topNode)
{
  KFParticle candidate;
  bool isGoodCandidate;

  std::tie(candidate, isGoodCandidate) = buildMother(vDaughters, daughterOrder, isIntermediate, intermediateNumber, nTracks, constrainMass, required_vertexID, topNode);

  if (constrain_to_vertex && isGoodCandidate && !isIntermediate)
  {
    if (m_require_track_and_vertex_match)
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

std::vector<std::vector<int>> KFParticle_Tools::findUniqueDaughterCombinations(int start, int end)
{
  std::vector<int> vect_permutations;
  std::vector<std::vector<int>> uniqueCombinations;
  std::map<int, int> daughterMap;
  for (int i = start; i < end; i++)
  {
    daughterMap.insert(std::pair<int, int>(i, abs(getParticleID(m_daughter_name[i].c_str()))));
    vect_permutations.push_back(i);
  }
  int *permutations = &vect_permutations[0];

  do
  {
    std::vector<int> combination;
    combination.reserve((end - start));
    for (int i = 0; i < (end - start); i++)
    {
      combination.push_back(daughterMap.find(permutations[i])->second);
    }
    uniqueCombinations.push_back(combination);
  } while (std::next_permutation(permutations, permutations + vect_permutations.size()));

  removeDuplicates(uniqueCombinations);

  return uniqueCombinations;
}

double KFParticle_Tools::calculateEllipsoidRadius(int posOrNeg, double sigma_ii, double sigma_jj, double sigma_ij)
{  // Note - Only works for a 2D ellipsoid OR rotated nD ellipsoid to avoid projections
  if (std::abs(posOrNeg) != 1)
  {
    std::cout << "You have set posOrNeg to " << posOrNeg << ". This value must be  +/- 1! Skipping" << std::endl;
    return 0;
  }

  double r_ij = sqrt((sigma_ii + sigma_jj) / 2 + posOrNeg * (sqrt(pow((sigma_ii - sigma_jj) / 2, 2) + pow(sigma_ij, 2))));

  return r_ij;
}

float KFParticle_Tools::calculateEllipsoidVolume(const KFParticle &particle)
{
  TMatrixD cov_matrix(3, 3);

  for (int i = 0; i < 3; ++i)
  {
    for (int j = 0; j < 3; ++j)
    {
      cov_matrix(i, j) = particle.GetCovariance(i, j);
    }
  }

  float volume;
  if (cov_matrix(0, 0) * cov_matrix(1, 1) * cov_matrix(2, 2) == 0)
  {
    volume = 0;
  }
  else
  {
    volume = (4. / 3.) * M_PI * sqrt((std::abs(cov_matrix.Determinant())));  // The covariance matrix is error-squared
  }

  return volume;
}

float KFParticle_Tools::calculateJT(const KFParticle &mother, const KFParticle &daughter)
{
  Eigen::Vector3f motherP = Eigen::Vector3f(mother.GetPx(), mother.GetPy(), mother.GetPz());
  Eigen::Vector3f daughterP = Eigen::Vector3f(daughter.GetPx(), daughter.GetPy(), daughter.GetPz());

  Eigen::Vector3f motherP_X_daughterP = motherP.cross(daughterP);

  float jT = (motherP_X_daughterP.norm()) / motherP.norm();

  return jT;
}

bool KFParticle_Tools::isInRange(float min, float value, float max)
{
  return min <= value && value <= max;
}

void KFParticle_Tools::removeDuplicates(std::vector<double> &v)
{
  auto end = v.end();
  for (auto it = v.begin(); it != end; ++it)
  {
    end = remove(it + 1, end, *it);
  }
  v.erase(end, v.end());
}

void KFParticle_Tools::removeDuplicates(std::vector<int> &v)
{
  auto end = v.end();
  for (auto it = v.begin(); it != end; ++it)
  {
    end = remove(it + 1, end, *it);
  }
  v.erase(end, v.end());
}

void KFParticle_Tools::removeDuplicates(std::vector<std::vector<int>> &v)
{
  auto end = v.end();
  for (auto it = v.begin(); it != end; ++it)
  {
    end = remove(it + 1, end, *it);
  }
  v.erase(end, v.end());
}

void KFParticle_Tools::removeDuplicates(std::vector<std::vector<std::string>> &v)
{
  auto end = v.end();
  for (auto it = v.begin(); it != end; ++it)
  {
    end = remove(it + 1, end, *it);
  }
  v.erase(end, v.end());
}

bool KFParticle_Tools::findParticle(const std::string &particle)
{
  bool particleFound = true;
  if (!TDatabasePDG::Instance()->GetParticle(particle.c_str()))
  {
    std::cout << "The particle, " << particle << " is not recognized" << std::endl;
    std::cout << "Check TDatabasePDG for a list of available particles" << std::endl;
    particleFound = false;
  }

  return particleFound;
}

int KFParticle_Tools::getParticleID(const std::string &particle)
{
  return TDatabasePDG::Instance()->GetParticle(particle.c_str())->PdgCode();
}

float KFParticle_Tools::getParticleMass(const std::string &particle)
{
  return TDatabasePDG::Instance()->GetParticle(particle.c_str())->Mass();
}

float KFParticle_Tools::getParticleMass(const int PDGID)
{
  return TDatabasePDG::Instance()->GetParticle(PDGID)->Mass();
}

void KFParticle_Tools::identify(const KFParticle &particle)
{
  std::cout << "Track ID: " << particle.Id() << std::endl;
  std::cout << "PDG ID: " << particle.GetPDG() << ", charge: " << (int) particle.GetQ() << ", mass: " << particle.GetMass() << " GeV" << std::endl;
  std::cout << "(px,py,pz) = (" << particle.GetPx() << " +/- " << std::sqrt(particle.GetCovariance(3, 3)) << ", ";
  std::cout << particle.GetPy() << " +/- " << std::sqrt(particle.GetCovariance(4, 4)) << ", ";
  std::cout << particle.GetPz() << " +/- " << std::sqrt(particle.GetCovariance(5, 5)) << ") GeV" << std::endl;
  std::cout << "(x,y,z) = (" << particle.GetX() << " +/- " << std::sqrt(particle.GetCovariance(0, 0)) << ", ";
  std::cout << particle.GetY() << " +/- " << std::sqrt(particle.GetCovariance(1, 1)) << ", ";
  std::cout << particle.GetZ() << " +/- " << std::sqrt(particle.GetCovariance(2, 2)) << ") cm\n"
            << std::endl;
}

float KFParticle_Tools::get_dEdx(PHCompositeNode *topNode, const KFParticle &daughter)
{
  m_dst_trackmap = findNode::getClass<SvtxTrackMap>(topNode, m_trk_map_node_name.c_str());
  m_cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  m_geom_container = findNode::getClass<PHG4TpcCylinderGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");
  if(!m_cluster_map || !m_geom_container)
  {
    std::cout << "Can't continue in KFParticle_Tools::get_dEdx, returning -1" << std::endl;
    return -1.0;
  }

  SvtxTrack *daughter_track = toolSet.getTrack(daughter.Id(), m_dst_trackmap);
  TrackSeed *tpcseed = daughter_track->get_tpc_seed();

  std::vector<TrkrDefs::cluskey> clusterKeys;
  clusterKeys.insert(clusterKeys.end(), tpcseed->begin_cluster_keys(), tpcseed->end_cluster_keys());
  std::vector<float> dedxlist;
  for (unsigned long cluster_key : clusterKeys)
  {
    unsigned int layer_local = TrkrDefs::getLayer(cluster_key);
    if(TrkrDefs::getTrkrId(cluster_key) != TrkrDefs::TrkrId::tpcId)
    {
        continue;
    }
    TrkrCluster* cluster = m_cluster_map->findCluster(cluster_key);

    float adc = cluster->getAdc();
    PHG4TpcCylinderGeom* GeoLayer_local = m_geom_container->GetLayerCellGeom(layer_local);
    float thick = GeoLayer_local->get_thickness();
    
    float r = GeoLayer_local->get_radius();
    float alpha = (r * r) / (2 * r * TMath::Abs(1.0 / tpcseed->get_qOverR()));
    float beta = atan(tpcseed->get_slope());

    float alphacorr = cos(alpha);
    if(alphacorr<0||alphacorr>4)
    {
      alphacorr=4;
    }

    float betacorr = cos(beta);
    if(betacorr<0||betacorr>4)
    {
      betacorr=4;
    }

    adc/=thick;
    adc*=alphacorr;
    adc*=betacorr;
    dedxlist.push_back(adc);
    sort(dedxlist.begin(), dedxlist.end());
  }

  int trunc_min = 0;
  int trunc_max = (int)dedxlist.size()*0.7;
  float sumdedx = 0;
  int ndedx = 0;
  for(int j = trunc_min; j<=trunc_max;j++)
  {
    sumdedx+=dedxlist.at(j);
    ndedx++;
  }
  sumdedx/=ndedx;
  return sumdedx;
}

void KFParticle_Tools::init_dEdx_fits()
{
  std::string dedx_fitparams = CDBInterface::instance()->getUrl("TPC_DEDX_FITPARAM");
  TFile *filefit = TFile::Open(dedx_fitparams.c_str());

  if (!filefit->IsOpen())
  {
      std::cerr << "Error opening filefit!" << std::endl;
      return;
  }

  filefit->GetObject("f_piband", f_pion_plus);
  filefit->GetObject("f_Kband", f_kaon_plus);
  filefit->GetObject("f_pband", f_proton_plus);
  filefit->GetObject("f_piminus_band", f_pion_minus);
  filefit->GetObject("f_Kminus_band", f_kaon_minus);
  filefit->GetObject("f_pbar_band", f_proton_minus);

  pidMap.insert(std::pair<int, TF1*>( -11,  f_pion_plus));
  pidMap.insert(std::pair<int, TF1*>( 211,  f_pion_plus));
  pidMap.insert(std::pair<int, TF1*>( 321,  f_kaon_plus));
  pidMap.insert(std::pair<int, TF1*>( 2212, f_proton_plus));
  pidMap.insert(std::pair<int, TF1*>(  11,  f_pion_minus));
  pidMap.insert(std::pair<int, TF1*>(-211,  f_pion_minus));
  pidMap.insert(std::pair<int, TF1*>(-321,  f_kaon_minus));
  pidMap.insert(std::pair<int, TF1*>(-2212, f_proton_minus));
}

double KFParticle_Tools::get_dEdx_fitValue(float momentum, int PID)
{
  return pidMap[PID]->Eval(momentum);
}

bool KFParticle_Tools::checkTrackAndVertexMatch(KFParticle vDaughters[], int nTracks, KFParticle vertex)
{
  bool vertexAndTrackMatch = true;

  int vertexCrossing = 1e5;

  if (m_dont_use_global_vertex)
  {
    if (m_use_mbd_vertex)
    {
      m_dst_mbdvertex = m_dst_mbdvertexmap->get(vertex.Id());
      vertexCrossing = m_dst_mbdvertex->get_beam_crossing();
    }
    else
    {
      m_dst_vertex = m_dst_vertexmap->get(vertex.Id());
      vertexCrossing = m_dst_vertex->get_beam_crossing();
    }
  }
  else
  {
    m_dst_globalvertex = m_dst_globalvertexmap->get(vertex.Id());

    GlobalVertex::VTXTYPE whichVtx = m_use_mbd_vertex ? GlobalVertex::MBD : GlobalVertex::SVTX; 

    auto svtxiter = m_dst_globalvertex->find_vertexes(whichVtx);
    auto svtxvertexvector = svtxiter->second;

    for (auto &gvertex : svtxvertexvector)
    {
      if (m_use_mbd_vertex)
      {
        m_dst_mbdvertex = m_dst_mbdvertexmap->find(gvertex->get_id())->second;
        vertexCrossing = m_dst_mbdvertex->get_beam_crossing();
      }
      else
      {
        m_dst_vertex = m_dst_vertexmap->find(gvertex->get_id())->second;
        vertexCrossing = m_dst_vertex->get_beam_crossing();
      }
    }
  }

  for (int i = 0; i < nTracks; ++i)
  {
    SvtxTrack *thisTrack = toolSet.getTrack(vDaughters[i].Id(), m_dst_trackmap);
    if (thisTrack)//This protects against intermediates which have no track
    {
      int trackCrossing = thisTrack->get_crossing();
      vertexAndTrackMatch = trackCrossing != vertexCrossing ? false : vertexAndTrackMatch;
    }
  }
   
  return vertexAndTrackMatch;
}
