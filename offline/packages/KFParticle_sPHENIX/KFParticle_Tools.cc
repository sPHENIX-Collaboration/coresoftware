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
#include <trackbase_historic/SvtxVertex.h>
#include <trackbase_historic/SvtxVertexMap.h>

//sPHENIX stuff
#include <g4eval/SvtxClusterEval.h>
#include <g4eval/SvtxEvalStack.h>
#include <g4eval/SvtxTrackEval.h>

#include <g4main/PHG4Particle.h>

#include <phool/getClass.h>

//KFParticle stuff
#include <KFPTrack.h>
#include <KFParticle.h>
#include <KFParticleDatabase.h>
#include <KFVertex.h>

#include <TMatrixD.h>

/// Create necessary objects
typedef pair<int, float> particle_pair;
KFParticle_particleList kfp_particleList;

//Particle masses are in GeV
map<string, particle_pair> particleMasses = kfp_particleList.getParticleList();

/// KFParticle constructor
KFParticle_Tools::KFParticle_Tools()
  : m_daughter_name{"pion", "pion", "pion", "pion"}
  , m_daughter_charge{1, -1, 1, -1}
  , m_min_mass(0)
  , m_max_mass(1e1)
  , m_track_pt(0.25)
  , m_track_ptchi2(FLT_MAX)
  , m_track_ipchi2(10.)
  , m_track_chi2ndof(4.)
  , m_comb_DCA(0.05)
  , m_vertex_chi2ndof(15.)
  , m_fdchi2(50.)
  , m_dira_min(0.95)
  , m_dira_max(1.01)
  , m_mother_pt(0.)
  , m_mother_ipchi2(FLT_MAX)
  , m_get_charge_conjugate(true)
  , m_vtx_map_node_name("SvtxVertexMap")
  , m_trk_map_node_name("SvtxTrackMap")
  , m_dst_vertexmap()
  , m_dst_trackmap()
  , m_dst_vertex()
  , m_dst_track()
{
}

KFParticle KFParticle_Tools::makeVertex(PHCompositeNode *topNode)
{
  KFParticle kfp_vertex;

  float f_vertexParameters[6] = {m_dst_vertex->get_x(),
                                 m_dst_vertex->get_y(),
                                 m_dst_vertex->get_z(), 0, 0, 0};

  float f_vertexCovariance[21];
  unsigned int iterate = 0;
  for (unsigned int i = 0; i < 3; ++i)
    for (unsigned int j = 0; j <= i; ++j)
    {
      f_vertexCovariance[iterate] = m_dst_vertex->get_error(i, j);
      ++iterate;
    }

  kfp_vertex.Create(f_vertexParameters, f_vertexCovariance, 0, -1);
  kfp_vertex.NDF() = m_dst_vertex->get_ndof();
  kfp_vertex.Chi2() = m_dst_vertex->get_chisq();

  return kfp_vertex;
}

vector<KFParticle> KFParticle_Tools::makeAllPrimaryVertices(PHCompositeNode *topNode)
{
  vector<KFParticle> primaryVertices;
  m_dst_vertexmap = findNode::getClass<SvtxVertexMap>(topNode, m_vtx_map_node_name.c_str());
  unsigned int vertexID = 0;

  for (SvtxVertexMap::ConstIter iter = m_dst_vertexmap->begin(); iter != m_dst_vertexmap->end(); ++iter)
  {
    m_dst_vertex = iter->second;
    primaryVertices.push_back(makeVertex(topNode));
    primaryVertices[vertexID].SetId(iter->first);
    ++vertexID;
  }

  return primaryVertices;
}

KFParticle KFParticle_Tools::makeParticle(PHCompositeNode *topNode)  ///Return a KFPTrack from track vector and covariance matrix. No mass or vertex constraints
{
  float f_trackParameters[6] = {m_dst_track->get_x(),
                                m_dst_track->get_y(),
                                m_dst_track->get_z(),
                                m_dst_track->get_px(),
                                m_dst_track->get_py(),
                                m_dst_track->get_pz()};

  float f_trackCovariance[21];
  unsigned int iterate = 0;
  for (unsigned int i = 0; i < 6; ++i)
    for (unsigned int j = 0; j <= i; ++j)
    {
      f_trackCovariance[iterate] = m_dst_track->get_error(i, j);
      ++iterate;
    }

  KFParticle kfp_particle;
  kfp_particle.Create(f_trackParameters, f_trackCovariance, (Int_t) m_dst_track->get_charge(), -1);
  kfp_particle.NDF() = m_dst_track->get_ndf();
  kfp_particle.Chi2() = m_dst_track->get_chisq();
  kfp_particle.SetId(m_dst_track->get_id());

  return kfp_particle;
}

vector<KFParticle> KFParticle_Tools::makeAllDaughterParticles(PHCompositeNode *topNode)
{
  vector<KFParticle> daughterParticles;
  m_dst_trackmap = findNode::getClass<SvtxTrackMap>(topNode, m_trk_map_node_name.c_str());
  unsigned int trackID = 0;

  for (SvtxTrackMap::Iter iter = m_dst_trackmap->begin(); iter != m_dst_trackmap->end(); ++iter)
  {
    m_dst_track = iter->second;
    daughterParticles.push_back(makeParticle(topNode));  ///Turn all dst tracks in KFP tracks
    daughterParticles[trackID].SetId(iter->first);
    ++trackID;
  }

  return daughterParticles;
}

int KFParticle_Tools::getTracksFromVertex(PHCompositeNode *topNode, KFParticle vertex)
{
  SvtxVertex *associatedVertex = NULL;
  m_dst_vertexmap = findNode::getClass<SvtxVertexMap>(topNode, m_vtx_map_node_name);

  if (m_dst_vertexmap->size() == 0) m_dst_vertexmap = findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMap");

  associatedVertex = m_dst_vertexmap->find(vertex.Id())->second;

  return associatedVertex->size_tracks();
}

const bool KFParticle_Tools::isGoodTrack(KFParticle particle, vector<KFParticle> primaryVertices)
{
  bool goodTrack = false;

  float pt = particle.GetPt();
  float pterr = particle.GetErrPt();
  float ptchi2 = pow(pterr / pt, 2);
  float trackchi2ndof = particle.GetChi2() / particle.GetNDF();
  vector<float> ipchi2;

  for (unsigned int i_verts = 0; i_verts < primaryVertices.size(); ++i_verts)
    ipchi2.push_back(particle.GetDeviationFromVertex(primaryVertices[i_verts]));

  auto minmax_ipchi2 = minmax_element(ipchi2.begin(), ipchi2.end());  //Order the IP chi2 from small to large
  float min_ipchi2 = *minmax_ipchi2.first;

  if (pt >= m_track_pt && ptchi2 <= m_track_ptchi2 && min_ipchi2 >= m_track_ipchi2 && trackchi2ndof <= m_track_chi2ndof) goodTrack = true;

  return goodTrack;
}

vector<int> KFParticle_Tools::findAllGoodTracks(vector<KFParticle> daughterParticles, const vector<KFParticle> primaryVertices)
{
  vector<int> goodTrackIndex;

  for (unsigned int i_parts = 0; i_parts < daughterParticles.size(); ++i_parts)
  {
    if (isGoodTrack(daughterParticles[i_parts], primaryVertices)) goodTrackIndex.push_back(i_parts);
  }

  removeDuplicates(goodTrackIndex);

  return goodTrackIndex;
}

vector<vector<int>> KFParticle_Tools::findTwoProngs(vector<KFParticle> daughterParticles, vector<int> goodTrackIndex, int nTracks)
{
  vector<vector<int>> goodTracksThatMeet;

  for (vector<int>::iterator i_it = goodTrackIndex.begin(); i_it != goodTrackIndex.end(); ++i_it)
  {
    for (vector<int>::iterator j_it = goodTrackIndex.begin(); j_it != goodTrackIndex.end(); ++j_it)
    {
      if (i_it < j_it)
      {
        if (daughterParticles[*i_it].GetDistanceFromParticle(daughterParticles[*j_it]) < m_comb_DCA)
        {
          KFVertex twoParticleVertex;
          twoParticleVertex += daughterParticles[*i_it];
          twoParticleVertex += daughterParticles[*j_it];
          float vertexchi2ndof = twoParticleVertex.GetChi2() / twoParticleVertex.GetNDF();
          vector<int> combination = {*i_it, *j_it};

          if (nTracks == 2 && vertexchi2ndof <= m_vertex_chi2ndof)
          {
            goodTracksThatMeet.push_back(combination);
          }
          else if (nTracks == 2 && vertexchi2ndof > m_vertex_chi2ndof)
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

vector<vector<int>> KFParticle_Tools::findNProngs(vector<KFParticle> daughterParticles,
                                                  vector<int> goodTrackIndex,
                                                  vector<vector<int>> goodTracksThatMeet,
                                                  int nRequiredTracks, unsigned int nProngs)
{
  unsigned int nGoodProngs = goodTracksThatMeet.size();

  for (vector<int>::iterator i_it = goodTrackIndex.begin(); i_it != goodTrackIndex.end(); ++i_it)
  {
    for (unsigned int i_prongs = 0; i_prongs < nGoodProngs; ++i_prongs)
    {
      bool trackNotUsedAlready = true;
      for (unsigned int i_trackCheck = 0; i_trackCheck < nProngs - 1; ++i_trackCheck)
        if (*i_it == goodTracksThatMeet[i_prongs][i_trackCheck]) trackNotUsedAlready = false;
      if (trackNotUsedAlready)
      {
        bool dcaMet = 1;
        for (unsigned int i = 0; i < nProngs - 1; ++i)
        {
          if (daughterParticles[*i_it].GetDistanceFromParticle(daughterParticles[goodTracksThatMeet[i_prongs][i]]) > m_comb_DCA)
          {
            dcaMet = 0;
          }
        }

        if (dcaMet)
        {
          KFVertex particleVertex;
          particleVertex += daughterParticles[*i_it];
          vector<int> combination;
          combination.push_back(*i_it);
          for (unsigned int i = 0; i < nProngs - 1; ++i)
          {
            particleVertex += daughterParticles[goodTracksThatMeet[i_prongs][i]];
            combination.push_back(goodTracksThatMeet[i_prongs][i]);
          }
          float vertexchi2ndof = particleVertex.GetChi2() / particleVertex.GetNDF();

          if ((unsigned int) nRequiredTracks == nProngs && vertexchi2ndof <= m_vertex_chi2ndof)
          {
            goodTracksThatMeet.push_back(combination);
          }
          else if ((unsigned int) nRequiredTracks == nProngs && vertexchi2ndof > m_vertex_chi2ndof)
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
  for (unsigned int i = 0; i < goodTracksThatMeet.size(); ++i) sort(goodTracksThatMeet[i].begin(), goodTracksThatMeet[i].end());
  removeDuplicates(goodTracksThatMeet);

  return goodTracksThatMeet;
}

vector<vector<int>> KFParticle_Tools::appendTracksToIntermediates(KFParticle intermediateResonances[], vector<KFParticle> daughterParticles, vector<int> goodTrackIndex, int num_remaining_tracks)
{
  vector<vector<int>> goodTracksThatMeet, goodTracksThatMeetIntermediates;  //, vectorOfGoodTracks;
  if (num_remaining_tracks == 1)
  {
    for (vector<int>::iterator i_it = goodTrackIndex.begin(); i_it != goodTrackIndex.end(); ++i_it)
    {
      vector<KFParticle> v_intermediateResonances(intermediateResonances, intermediateResonances + m_num_intermediate_states);
      vector<vector<int>> dummyTrackList;
      vector<int> dummyTrackID;  //I already have the track ids stored in goodTracksThatMeet[i]
      v_intermediateResonances.insert(end(v_intermediateResonances), daughterParticles[*i_it]);
      for (unsigned int k = 0; k < v_intermediateResonances.size(); ++k) 
      {
        dummyTrackID.push_back(k);
      }
      dummyTrackList = findTwoProngs(v_intermediateResonances, dummyTrackID, (int) v_intermediateResonances.size());
      if (v_intermediateResonances.size() > 2)
      {
        for (unsigned int p = 3; p <= v_intermediateResonances.size(); ++p) dummyTrackList = findNProngs(v_intermediateResonances,
                                                                                                         dummyTrackID, dummyTrackList,
                                                                                                         (int) v_intermediateResonances.size(), (int) p);
      }

      if (dummyTrackList.size() != 0)
      {
        vector<int> goodTrack{*i_it};
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

    for (unsigned int i = 0; i < goodTracksThatMeet.size(); ++i)
    {
      vector<KFParticle> v_intermediateResonances(intermediateResonances, intermediateResonances + m_num_intermediate_states);
      vector<vector<int>> dummyTrackList;
      vector<int> dummyTrackID;  //I already have the track ids stored in goodTracksThatMeet[i]
      for (unsigned int j = 0; j < goodTracksThatMeet[i].size(); ++j) 
      {
        v_intermediateResonances.push_back(daughterParticles[goodTracksThatMeet[i][j]]);
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

      if (dummyTrackList.size() != 0) goodTracksThatMeetIntermediates.push_back(goodTracksThatMeet[i]);
    }
  }

  return goodTracksThatMeetIntermediates;
}

float KFParticle_Tools::eventDIRA(KFParticle particle, KFParticle vertex)
{
  TMatrixD flightVector(3, 1);
  TMatrixD momVector(3, 1);
  flightVector(0, 0) = particle.GetX() - vertex.GetX();
  flightVector(1, 0) = particle.GetY() - vertex.GetY();
  flightVector(2, 0) = particle.GetZ() - vertex.GetZ();

  momVector(0, 0) = particle.GetPx();
  momVector(1, 0) = particle.GetPy();
  momVector(2, 0) = particle.GetPz();

  TMatrixD momDotFD(1, 1);  //Calculate momentum dot flight distance
  momDotFD = TMatrixD(momVector, TMatrixD::kTransposeMult, flightVector);
  float f_momDotFD = momDotFD(0, 0);

  TMatrixD sizeOfMom(1, 1);  //Calculates the size of the momentum vector
  sizeOfMom = TMatrixD(momVector, TMatrixD::kTransposeMult, momVector);
  float f_sizeOfMom = sqrt(sizeOfMom(0, 0));

  TMatrixD sizeOfFD(1, 1);  //Calculates the size of the flight distance vector
  sizeOfFD = TMatrixD(flightVector, TMatrixD::kTransposeMult, flightVector);
  float f_sizeOfFD = sqrt(sizeOfFD(0, 0));

  return f_momDotFD / (f_sizeOfMom * f_sizeOfFD);  //returns the DIRA
}

float KFParticle_Tools::flightDistanceChi2(KFParticle particle, KFParticle vertex)
{
  TMatrixD flightVector(3, 1);
  TMatrixD flightDistanceCovariance(3, 3);

  KFParticle kfp_vertex(vertex);

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

tuple<KFParticle, bool> KFParticle_Tools::buildMother(KFParticle vDaughters[], string daughterOrder[], bool isIntermediate, int intermediateNumber, int nTracks, bool constrainMass, float required_vertexID)
{
  KFParticle mother;
  KFParticle inputTracks[nTracks];

  mother.SetConstructMethod(2);

  bool daughterMassCheck = true;
  float unique_vertexID = 0;

  for (int i = 0; i < nTracks; ++i)
  {
    float daughterMass = constrainMass ? particleMasses.find(daughterOrder[i].c_str())->second.second : vDaughters[i].GetMass();
    inputTracks[i].Create(vDaughters[i].Parameters(),
                          vDaughters[i].CovarianceMatrix(),
                          (Int_t) vDaughters[i].GetQ(),
                          daughterMass);
    mother.AddDaughter(inputTracks[i]);
    if (inputTracks[i].GetMass() == 0) daughterMassCheck = false;
    unique_vertexID += vDaughters[i].GetQ() * particleMasses.find(daughterOrder[i].c_str())->second.second;
  }

  if (isIntermediate) mother.SetPDG(particleMasses.find(m_intermediate_name[intermediateNumber].c_str())->second.first);
  if (!isIntermediate && !m_mother_name_Tools.empty()) mother.SetPDG(particleMasses.find(m_mother_name_Tools.c_str())->second.first);

  bool chargeCheck;
  if (m_get_charge_conjugate)
    chargeCheck = abs(unique_vertexID) == abs(required_vertexID) ? 1 : 0;
  else
    chargeCheck = unique_vertexID == required_vertexID ? 1 : 0;

  for (int j = 0; j < nTracks; ++j) inputTracks[j].SetProductionVertex(mother);

  float calculated_mass, calculated_mass_err;
  mother.GetMass(calculated_mass, calculated_mass_err);
  float calculated_pt = mother.GetPt();

  float min_mass = isIntermediate ? m_intermediate_mass_range[intermediateNumber].first : m_min_mass;
  float max_mass = isIntermediate ? m_intermediate_mass_range[intermediateNumber].second : m_max_mass;
  float min_pt = isIntermediate ? m_intermediate_min_pt[intermediateNumber] : m_mother_pt;

  bool goodCandidate = false;
  if (calculated_mass >= min_mass && calculated_mass <= max_mass &&
      calculated_pt >= min_pt &&
      daughterMassCheck && chargeCheck)
    goodCandidate = true;

  return make_tuple(mother, goodCandidate);
}

void KFParticle_Tools::constrainToVertex(KFParticle &particle, bool &goodCandidate, KFParticle &vertex)
{
  //KFParticle prod_vertex( vertex );
  //prod_vertex.AddDaughter( particle );
  //particle.SetProductionVertex( prod_vertex );
  //particle.SetAtProductionVertex( true );

  float calculated_fdchi2 = flightDistanceChi2(particle, vertex);
  float calculated_dira = eventDIRA(particle, vertex);
  //float calculated_ipchi2  = particle.GetDeviationFromVertex( prod_vertex );
  float calculated_ipchi2 = particle.GetDeviationFromVertex(vertex);
  //float calculated_lifetime, calculated_lifetime_error;
  //mother.GetLifeTime( calculated_lifetime, calculated_lifetime_error );
  goodCandidate = false;
  if (calculated_fdchi2 >= m_fdchi2 && calculated_ipchi2 <= m_mother_ipchi2 &&
      calculated_dira >= m_dira_min && calculated_dira <= m_dira_max) goodCandidate = true;
}

tuple<KFParticle, bool> KFParticle_Tools::getCombination(KFParticle vDaughters[], string daughterOrder[], KFParticle vertex, bool constrain_to_vertex, bool isIntermediate, int intermediateNumber, int nTracks, bool constrainMass, float required_vertexID)
{
  KFParticle candidate;
  bool isGoodCandidate;

  tie(candidate, isGoodCandidate) = buildMother(vDaughters, daughterOrder, isIntermediate, intermediateNumber, nTracks, constrainMass, required_vertexID);
  if (constrain_to_vertex && isGoodCandidate && !isIntermediate) constrainToVertex(candidate, isGoodCandidate, vertex);

  return make_tuple(candidate, isGoodCandidate);
}

vector<vector<string>> KFParticle_Tools::findUniqueDaughterCombinations(int start, int end)
{
  vector<int> vect_permutations;
  vector<vector<string>> uniqueCombinations;
  map<int, string> daughterMap;
  for (int i = start; i < end; i++)
  {
    daughterMap.insert(pair<int, string>(i, m_daughter_name[i].c_str()));
    vect_permutations.push_back(i);
  }
  int *permutations = &vect_permutations[0];

  do
  {
    vector<string> combination;
    for (int i = 0; i < (end - start); i++) combination.push_back(daughterMap.find(permutations[i])->second);
    uniqueCombinations.push_back(combination);
  } while (next_permutation(permutations, permutations + vect_permutations.size()));

  removeDuplicates(uniqueCombinations);

  return uniqueCombinations;
}

double KFParticle_Tools::calculateEllipsoidRadius(int posOrNeg, double sigma_ii, double sigma_jj, double sigma_ij)
{  //Note - Only works for a 2D ellipsoid OR rotated nD ellipsoid to avoid projections
  if (abs(posOrNeg) != 1)
  {
    cout << "You have set posOrNeg to " << posOrNeg << ". This value must be  +/- 1! Skipping" << endl;
    return 0;
  }

  double r_ij = sqrt((sigma_ii + sigma_jj) / 2 + posOrNeg * (sqrt(pow((sigma_ii - sigma_jj) / 2, 2) + pow(sigma_ij, 2))));

  return r_ij;
}

float KFParticle_Tools::calculateEllipsoidVolume(KFParticle particle)
{
  TMatrixD cov_matrix(3, 3);

  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      cov_matrix(i, j) = particle.GetCovariance(i, j);

  float volume;
  if (cov_matrix(0, 0) * cov_matrix(1, 1) * cov_matrix(2, 2) == 0)
    volume = 0;
  else
    volume = (4 / 3) * M_PI * sqrt((abs(cov_matrix.Determinant())));  //The covariance matrix is error-squared

  return volume;
}

void KFParticle_Tools::removeDuplicates(vector<double> &v)
{
  auto end = v.end();
  for (auto it = v.begin(); it != end; ++it)
  {
    end = remove(it + 1, end, *it);
  }
  v.erase(end, v.end());
}

void KFParticle_Tools::removeDuplicates(vector<int> &v)
{
  auto end = v.end();
  for (auto it = v.begin(); it != end; ++it)
  {
    end = remove(it + 1, end, *it);
  }
  v.erase(end, v.end());
}

void KFParticle_Tools::removeDuplicates(vector<vector<int>> &v)
{
  auto end = v.end();
  for (auto it = v.begin(); it != end; ++it)
  {
    end = remove(it + 1, end, *it);
  }
  v.erase(end, v.end());
}

void KFParticle_Tools::removeDuplicates(vector<vector<string>> &v)
{
  auto end = v.end();
  for (auto it = v.begin(); it != end; ++it)
  {
    end = remove(it + 1, end, *it);
  }
  v.erase(end, v.end());
}

void KFParticle_Tools::identify(KFParticle particle)
{
  cout << "Track ID: " << particle.Id() << endl;
  cout << "PDG ID: " << particle.GetPDG() << ", charge: " << (int) particle.GetQ() << ", mass: " << particle.GetMass() << " GeV" << endl;
  cout << "(px,py,pz) = (" << particle.GetPx() << " +/- " << sqrt(particle.GetCovariance(3, 3)) << ", ";
  cout << particle.GetPy() << " +/- " << sqrt(particle.GetCovariance(4, 4)) << ", ";
  cout << particle.GetPz() << " +/- " << sqrt(particle.GetCovariance(5, 5)) << ") GeV" << endl;
  cout << "(x,y,z) = (" << particle.GetX() << " +/- " << sqrt(particle.GetCovariance(0, 0)) << ", ";
  cout << particle.GetY() << " +/- " << sqrt(particle.GetCovariance(1, 1)) << ", ";
  cout << particle.GetZ() << " +/- " << sqrt(particle.GetCovariance(2, 2)) << ") cm\n"
       << endl;
}
