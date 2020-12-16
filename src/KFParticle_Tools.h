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
#include <fun4all/Fun4AllReturnCodes.h>
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

using namespace std;

class PHCompositeNode;
class SvtxVertexMap;
class SvtxTrackMap;
class SvtxVertex;
class SvtxTrack;

class KFParticle_Tools : public KFParticle_particleList, protected KFParticle_MVA 
{
 public:
  
  KFParticle_Tools();

  ~KFParticle_Tools();

  KFParticle makeVertex(PHCompositeNode *topNode);

  vector<KFParticle> makeAllPrimaryVertices(PHCompositeNode *topNode);

  KFParticle makeParticle(PHCompositeNode *topNode);

  vector<KFParticle> makeAllDaughterParticles(PHCompositeNode *topNode);

  int getTracksFromVertex( PHCompositeNode *topNode, KFParticle vertex );

  const bool isGoodTrack(KFParticle particle, vector<KFParticle> primaryVertices);

  vector<int> findAllGoodTracks(vector<KFParticle> daughterParticles, vector<KFParticle> primaryVertices);

  vector<vector<int>> findTwoProngs(vector<KFParticle> daughterParticles, vector<int> goodTrackIndex, int nTracks);

  vector<vector<int>> findNProngs( vector<KFParticle> daughterParticles, 
                                             vector<int> goodTrackIndex, 
                                             vector<vector<int>> goodTracksThatMeet, 
                                             int nRequiredTracks, unsigned int nProngs );

  vector<vector<int>> appendTracksToIntermediates( KFParticle intermediateResonances[], vector<KFParticle> daughterParticles, vector<int> goodTrackIndex, int num_remaining_tracks);

  float eventDIRA(KFParticle particle, KFParticle vertex);

  float flightDistanceChi2(KFParticle particle, KFParticle vertex);

  tuple<KFParticle, bool> buildMother( KFParticle vDaughters[], string daughterOrder[], bool isIntermediate, int intermediateNumber, int nTracks, bool constrainMass, float required_vertexID );

  void constrainToVertex(KFParticle& particle, bool& goodCandidate, KFParticle& vertex );

  tuple<KFParticle, bool> getCombination(KFParticle vDaughters[], string daughterOrder[], KFParticle vertex, 
                                              bool constrain_to_vertex, bool isIntermediate, int intermediateNumber, int nTracks, bool constrainMass, float required_vertexID );

  vector<vector<string>> findUniqueDaughterCombinations( int start, int end );

  double calculateEllipsoidRadius( int posOrNeg, double sigma_ii, double sigma_jj, double sigma_ij );
  
  float calculateEllipsoidVolume( KFParticle particle );

  void identify( KFParticle particle );
 
 protected:
   
  string m_mother_name_Tools;
  int m_num_intermediate_states;
  int m_num_tracks_from_intermediate[99];
  string m_daughter_name[99];
  int m_daughter_charge[99];
  string m_intermediate_name[99];
  float m_min_mass;
  float m_max_mass;
  pair<float, float> m_intermediate_mass_range[99];
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
  bool m_get_charge_conjugate;
  string m_vtx_map_node_name;
  string m_trk_map_node_name;
  SvtxVertexMap *m_dst_vertexmap;
  SvtxTrackMap *m_dst_trackmap;
  SvtxVertex *m_dst_vertex;
  SvtxTrack *m_dst_track;

 private:
 
  void removeDuplicates(vector<double> &v);
  void removeDuplicates(vector<int> &v);
  void removeDuplicates(vector<vector<int>> &v);
  void removeDuplicates(vector<vector<string>> &v);
};

#endif //KFParticle_Tools_H

