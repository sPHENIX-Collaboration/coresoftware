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

#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMap.h>
#include <globalvertex/SvtxVertex.h>
#include <globalvertex/SvtxVertexMap.h>

#include <phool/getClass.h>

// KFParticle stuff
#include <KFParticle.h>
#include <KFVertex.h>

#include <Rtypes.h>
#include <TDatabasePDG.h>
#include <TMatrixD.h>
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
  , m_track_ip(-1.)
  , m_track_ipchi2(0.)
  , m_track_chi2ndof(4.)
  , m_nMVTXHits(3)
  , m_nTPCHits(20)
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
  , m_dst_vertexmap()
  , m_dst_trackmap()
  , m_dst_vertex()
  , m_dst_track()
{
}

KFParticle KFParticle_Tools::makeVertex(PHCompositeNode * /*topNode*/)
{
  float f_vertexParameters[6] = {m_dst_vertex->get_x(),
                                 m_dst_vertex->get_y(),
                                 m_dst_vertex->get_z(), 0, 0, 0};

  float f_vertexCovariance[21];
  unsigned int iterate = 0;
  for (unsigned int i = 0; i < 3; ++i)
  {
    for (unsigned int j = 0; j <= i; ++j)
    {
      f_vertexCovariance[iterate] = m_dst_vertex->get_error(i, j);
      ++iterate;
    }
  }

  KFParticle kfp_vertex;
  kfp_vertex.Create(f_vertexParameters, f_vertexCovariance, 0, -1);
  kfp_vertex.NDF() = m_dst_vertex->get_ndof();
  kfp_vertex.Chi2() = m_dst_vertex->get_chisq();

  return kfp_vertex;
}

std::vector<KFParticle> KFParticle_Tools::makeAllPrimaryVertices(PHCompositeNode *topNode, const std::string &vertexMapName)
{
  std::string vtxMN;
  if (vertexMapName.empty())
  {
    vtxMN = m_vtx_map_node_name;
  }
  else
  {
    vtxMN = vertexMapName;
  }

  std::vector<KFParticle> primaryVertices;
  m_dst_vertexmap = findNode::getClass<SvtxVertexMap>(topNode, vtxMN);
  auto globalvertexmap = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");
  if (!globalvertexmap)
  {
    std::cout << "Can't continue in KFParticle_Tools::makeAllPrimaryVertices" << std::endl;
  }

  unsigned int vertexID = 0;

  for (GlobalVertexMap::ConstIter iter = globalvertexmap->begin(); iter != globalvertexmap->end(); ++iter)
  {
    GlobalVertex *gvertex = iter->second;
    auto svtxiter = gvertex->find_vertexes(GlobalVertex::SVTX);
    // check that it contains a track vertex
    if (svtxiter == gvertex->end_vertexes())
    {
      continue;
    }

    auto svtxvertexvector = svtxiter->second;

    for (auto &vertex : svtxvertexvector)
    {
      m_dst_vertex = m_dst_vertexmap->find(vertex->get_id())->second;

      primaryVertices.push_back(makeVertex(topNode));
      primaryVertices[vertexID].SetId(gvertex->get_id());
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

    // First check if we have the required number of MVTX and TPC hits
    TrackSeed *tpcseed = m_dst_track->get_tpc_seed();
    TrackSeed *silseed = m_dst_track->get_silicon_seed();
    int MVTX_hits = 0;
    int TPC_hits = 0;

    if (silseed)
    {
      for (auto cluster_iter = silseed->begin_cluster_keys();
           cluster_iter != silseed->end_cluster_keys(); ++cluster_iter)
      {
        const auto &cluster_key = *cluster_iter;
        const auto trackerID = TrkrDefs::getTrkrId(cluster_key);

        if (trackerID == TrkrDefs::mvtxId)
        {
          ++MVTX_hits;
        }
      }
      if (MVTX_hits < m_nMVTXHits)
      {
        continue;
      }
    }
    if (tpcseed)
    {
      for (auto cluster_iter = tpcseed->begin_cluster_keys(); cluster_iter != tpcseed->end_cluster_keys(); ++cluster_iter)
      {
        const auto &cluster_key = *cluster_iter;
        const auto trackerID = TrkrDefs::getTrkrId(cluster_key);

        if (trackerID == TrkrDefs::tpcId)
        {
          ++TPC_hits;
        }
      }
      if (TPC_hits < m_nTPCHits)
      {
        continue;
      }
    }

    daughterParticles.push_back(makeParticle(topNode));  /// Turn all dst tracks in KFP tracks
    daughterParticles[trackID].SetId(iter.first);
    ++trackID;
  }

  return daughterParticles;
}

int KFParticle_Tools::getTracksFromVertex(PHCompositeNode *topNode, const KFParticle &vertex, const std::string &vertexMapName)
{
  std::string vtxMN;
  if (vertexMapName.empty())
  {
    vtxMN = m_vtx_map_node_name;
  }
  else
  {
    vtxMN = vertexMapName;
  }

  SvtxVertex *associatedVertex = nullptr;
  m_dst_vertexmap = findNode::getClass<SvtxVertexMap>(topNode, vtxMN);
  auto globalvertexmap = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");
  GlobalVertex *associatedgvertex = globalvertexmap->find(vertex.Id())->second;
  auto svtxvtx_id = associatedgvertex->find_vtxids(GlobalVertex::SVTX)->second;
  associatedVertex = m_dst_vertexmap->find(svtxvtx_id)->second;

  return associatedVertex->size_tracks();
}

/*const*/ bool KFParticle_Tools::isGoodTrack(const KFParticle &particle, const std::vector<KFParticle> &primaryVertices)
{
  bool goodTrack = false;

  float min_ip = 0;
  float min_ipchi2 = 0;

  float pt = particle.GetPt();
  float pterr = particle.GetErrPt();
  float ptchi2 = pow(pterr / pt, 2);
  float trackchi2ndof = particle.GetChi2() / particle.GetNDF();
  calcMinIP(particle, primaryVertices, min_ip, min_ipchi2);

  if (pt >= m_track_pt && ptchi2 <= m_track_ptchi2 && min_ip >= m_track_ip && min_ipchi2 >= m_track_ipchi2 && trackchi2ndof <= m_track_chi2ndof)
  {
    goodTrack = true;
  }
  return goodTrack;
}

int KFParticle_Tools::calcMinIP(const KFParticle &track, const std::vector<KFParticle> &PVs,
                                float &minimumIP, float &minimumIPchi2)
{
  std::vector<float> ip, ipchi2;

  for (const auto &PV : PVs)
  {
    ip.push_back(track.GetDistanceFromVertex(PV));
    float thisIPchi2 = track.GetDeviationFromVertex(PV);
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

std::vector<int> KFParticle_Tools::findAllGoodTracks(std::vector<KFParticle> daughterParticles, const std::vector<KFParticle> &primaryVertices)
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
        if (daughterParticles[*i_it].GetDistanceFromParticle(daughterParticles[*j_it]) <= m_comb_DCA)
        {
          KFVertex twoParticleVertex;
          twoParticleVertex += daughterParticles[*i_it];
          twoParticleVertex += daughterParticles[*j_it];
          float vertexchi2ndof = twoParticleVertex.GetChi2() / twoParticleVertex.GetNDF();
          std::vector<int> combination = {*i_it, *j_it};

          if (nTracks == 2 && vertexchi2ndof > m_vertex_chi2ndof)
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
          if (daughterParticles[i_it].GetDistanceFromParticle(daughterParticles[goodTracksThatMeet[i_prongs][i]]) > m_comb_DCA)
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

          if ((unsigned int) nRequiredTracks == nProngs && vertexchi2ndof > m_vertex_chi2ndof)
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

  goodTracksThatMeet.erase(goodTracksThatMeet.begin(), goodTracksThatMeet.begin() + nGoodProngs);
  for (auto &i : goodTracksThatMeet)
  {
    sort(i.begin(), i.end());
  }
  removeDuplicates(goodTracksThatMeet);

  return goodTracksThatMeet;
}

std::vector<std::vector<int>> KFParticle_Tools::appendTracksToIntermediates(KFParticle intermediateResonances[], std::vector<KFParticle> daughterParticles, const std::vector<int> &goodTrackIndex, int num_remaining_tracks)
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

float KFParticle_Tools::eventDIRA(const KFParticle &particle, const KFParticle &vertex)
{
  TMatrixD flightVector(3, 1);
  TMatrixD momVector(3, 1);
  flightVector(0, 0) = particle.GetX() - vertex.GetX();
  flightVector(1, 0) = particle.GetY() - vertex.GetY();
  flightVector(2, 0) = particle.GetZ() - vertex.GetZ();

  momVector(0, 0) = particle.GetPx();
  momVector(1, 0) = particle.GetPy();
  momVector(2, 0) = particle.GetPz();

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
                                                           bool constrainMass, float required_vertexID)
{
  KFParticle mother;
  KFParticle inputTracks[nTracks];

  mother.SetConstructMethod(2);

  bool daughterMassCheck = true;
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
    mother.AddDaughter(inputTracks[i]);
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

  particleCopy.GetLifeTime(calculated_decayTime, calculated_decayTimeErr);
  particleCopy.GetDecayLength(calculated_decayLength, calculated_decayLengthErr);

  float calculated_fdchi2 = flightDistanceChi2(particle, vertex);
  float calculated_dira = eventDIRA(particle, vertex);
  float calculated_ipchi2 = particle.GetDeviationFromVertex(vertex);

  goodCandidate = false;

  const float speed = 2.99792458e-2;
  calculated_decayTime /= speed;

  if (calculated_fdchi2 >= m_fdchi2 && calculated_ipchi2 <= m_mother_ipchi2 && isInRange(m_dira_min, calculated_dira, m_dira_max) && isInRange(m_min_decayTime, calculated_decayTime, m_max_decayTime) && isInRange(m_min_decayLength, calculated_decayLength, m_max_decayLength))
  {
    goodCandidate = true;
  }
}

std::tuple<KFParticle, bool> KFParticle_Tools::getCombination(KFParticle vDaughters[], int daughterOrder[], KFParticle vertex, bool constrain_to_vertex, bool isIntermediate, int intermediateNumber, int nTracks, bool constrainMass, float required_vertexID)
{
  KFParticle candidate;
  bool isGoodCandidate;

  std::tie(candidate, isGoodCandidate) = buildMother(vDaughters, daughterOrder, isIntermediate, intermediateNumber, nTracks, constrainMass, required_vertexID);

  if (constrain_to_vertex && isGoodCandidate && !isIntermediate)
  {
    constrainToVertex(candidate, isGoodCandidate, vertex);
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
