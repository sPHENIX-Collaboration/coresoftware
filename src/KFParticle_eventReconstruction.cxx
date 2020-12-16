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

#include "KFParticle_eventReconstruction.h"

//KFParticle stuff
#include "KFPTrack.h"
#include "KFParticle.h"
#include "KFVertex.h"
#include "KFParticleDatabase.h"

//#include <g4main/PHG4Particle.h>
#include <g4eval/SvtxTrackEval.h>
#include <g4eval/SvtxClusterEval.h>
#include <g4eval/SvtxEvalStack.h>

using namespace std;

/// Create necessary objects 
typedef pair<int, float> particle_pair;
KFParticle_particleList kfp_particleList_evtReco;

//Particle masses are in GeV
map<string, particle_pair> particleMasses_evtReco = kfp_particleList_evtReco.getParticleList(); 

/// KFParticle constructor
KFParticle_eventReconstruction::KFParticle_eventReconstruction():
    m_has_intermediates( false ),
    m_num_tracks( 2 ),
    m_daughter_name_evt{"pion", "pion", "pion", "pion"},
    m_daughter_charge_evt{ 1, -1, 1, -1}, 
    m_constrain_to_vertex( true ),
    m_constrain_int_mass( false ) 
{}

KFParticle_eventReconstruction::~KFParticle_eventReconstruction(){} /// KFParticle destructor


void KFParticle_eventReconstruction::createDecay( PHCompositeNode *topNode, vector<KFParticle>& selectedMother, vector<KFParticle>& selectedVertex,
                                                              vector<vector<KFParticle>>& selectedDaughters,
                                                              vector<vector<KFParticle>>& selectedIntermediates,
                                                              int& nPVs, int& multiplicity)
{
  KFParticle::SetField( -1.5e0 );

  vector<KFParticle> primaryVertices = makeAllPrimaryVertices( topNode );
  vector<KFParticle> daughterParticles = makeAllDaughterParticles( topNode );

  nPVs = primaryVertices.size();
  multiplicity = daughterParticles.size();

  vector<int> goodTrackIndex = findAllGoodTracks( daughterParticles, primaryVertices );

  if ( !m_has_intermediates ) buildBasicChain(selectedMother, selectedVertex,  selectedDaughters, daughterParticles, goodTrackIndex, primaryVertices);
  else buildChain(selectedMother, selectedVertex,  selectedDaughters, selectedIntermediates, daughterParticles, goodTrackIndex, primaryVertices);
 }

/*
 *  This function is used to build a basic n-body decay without any intermediate particles such as D's or J/psi's
 */
void KFParticle_eventReconstruction::buildBasicChain(vector<KFParticle>& selectedMother, 
	                               vector<KFParticle>& selectedVertex, 
	                               vector<vector<KFParticle>>& selectedDaughters, 
                                       const vector<KFParticle> daughterParticles,
                                       const vector<int> goodTrackIndex,
	                               const vector<KFParticle> primaryVertices)
{
  vector<vector<int>> goodTracksThatMeet = findTwoProngs( daughterParticles, goodTrackIndex, m_num_tracks );
  for ( int p = 3; p < m_num_tracks + 1; ++p) goodTracksThatMeet = findNProngs( daughterParticles, goodTrackIndex, goodTracksThatMeet, m_num_tracks, p);

  getCandidateDecay(selectedMother, selectedVertex, selectedDaughters, daughterParticles,
                    goodTracksThatMeet, primaryVertices, 0, m_num_tracks, false, 0, true);
}

/*
 *  This function is used to build a more complicated decay with intermediate particles such as D's or J/psi's
 */
void KFParticle_eventReconstruction::buildChain(vector<KFParticle>& selectedMother,
                                  vector<KFParticle>& selectedVertex,
                                  vector<vector<KFParticle>>& selectedDaughters,
                                  vector<vector<KFParticle>>& selectedIntermediates,
                                  vector<KFParticle> daughterParticles,
                                  const vector<int> goodTrackIndex,
                                  vector<KFParticle> primaryVertices)
{
    int track_start = 0;
    int track_stop = m_num_tracks_from_intermediate[0];

    vector<KFParticle> goodCandidates, goodVertex, goodDaughters[ m_num_tracks ], goodIntermediates[ m_num_intermediate_states ];
    vector<KFParticle> potentialIntermediates[ m_num_intermediate_states ];
    vector<vector<KFParticle>> potentialDaughters[ m_num_intermediate_states ];

    for (int i = 0; i < m_num_intermediate_states; ++i) 
    {
        vector<KFParticle> vertices;

    	vector<vector<int>> goodTracksThatMeet = findTwoProngs( daughterParticles, goodTrackIndex, m_num_tracks_from_intermediate[i] );
        for ( int p = 3; p <= m_num_tracks_from_intermediate[i]; ++p) goodTracksThatMeet = findNProngs( daughterParticles, 
                                                                                                           goodTrackIndex, 
                                                                                                           goodTracksThatMeet, 
                                                                                                           m_num_tracks_from_intermediate[i], p);

        getCandidateDecay(potentialIntermediates[i], vertices, potentialDaughters[i], daughterParticles, goodTracksThatMeet, primaryVertices, track_start, track_stop, true, i, m_constrain_int_mass);

 	track_start += track_stop;
 	track_stop += m_num_tracks_from_intermediate[ i + 1 ];
    }
    int num_tracks_used_by_intermediates = 0; 
    for (int i = 0; i < m_num_intermediate_states; ++i) num_tracks_used_by_intermediates += m_num_tracks_from_intermediate[i];
    int num_remaining_tracks = m_num_tracks - num_tracks_used_by_intermediates;   
    unsigned int num_pot_inter_a, num_pot_inter_b, num_pot_inter_c, num_pot_inter_d; //Number of potential intermediates found
    num_pot_inter_a = potentialIntermediates[0].size();
    num_pot_inter_b = m_num_intermediate_states < 2 ? 1 : potentialIntermediates[1].size(); //Ensure the code inside the loop below is executed
    num_pot_inter_c = m_num_intermediate_states < 3 ? 1 : potentialIntermediates[2].size();
    num_pot_inter_d = m_num_intermediate_states < 4 ? 1 : potentialIntermediates[3].size();

    for (unsigned int a = 0; a < num_pot_inter_a; ++a)
    {
    for (unsigned int b = 0; b < num_pot_inter_b; ++b)
    {  
    for (unsigned int c = 0; c < num_pot_inter_c; ++c)
    {
    for (unsigned int d = 0; d < num_pot_inter_d; ++d)
    {  
    for (unsigned int i_pv = 0; i_pv < primaryVertices.size(); ++i_pv)
    {  
        KFParticle candidate;
        bool isGood = false;
        unsigned int matchIterators[4] = { a, b, c, d };

        int num_mother_decay_products = m_num_intermediate_states + num_remaining_tracks;
        KFParticle motherDecayProducts[ num_mother_decay_products ];
        vector<KFParticle> finalTracks = potentialDaughters[0][a]; 
        for (int i = 0; i < m_num_intermediate_states; ++i) motherDecayProducts[i] = potentialIntermediates[i][matchIterators[i]]; 
        for (int j = 1; j < m_num_intermediate_states; ++j) finalTracks.insert(finalTracks.end(), potentialDaughters[j][matchIterators[j]].begin(), potentialDaughters[j][matchIterators[j]].end());

        // If there are daughter tracks coming from the mother not an intermediate, need to ensure that the intermeditate decay tracks aren't used again
        vector<int> goodTrackIndex_withoutIntermediates = goodTrackIndex;
        for (int m = 0; m < m_num_intermediate_states; ++m) 
        { 
          int trackID_to_remove =  finalTracks[m].Id();
          goodTrackIndex_withoutIntermediates.erase(remove(goodTrackIndex_withoutIntermediates.begin(),
                                                                goodTrackIndex_withoutIntermediates.end(), trackID_to_remove),
                                                                goodTrackIndex_withoutIntermediates.end()); 
        }

        float required_unique_vertexID = 0;
        for ( int n = 0; n < m_num_intermediate_states; ++n ) 
           required_unique_vertexID += m_intermediate_charge[n]*particleMasses_evtReco.find( m_intermediate_name[n].c_str() )->second.second ;

        if (num_remaining_tracks == 0) tie( candidate, isGood ) = getCombination( motherDecayProducts, m_intermediate_name, primaryVertices[i_pv], 
                                                                       m_constrain_to_vertex, false, 0, num_mother_decay_products, m_constrain_int_mass, required_unique_vertexID );
        else//Build n-prong from remaining tracks if needed
        {
          for ( int i = num_tracks_used_by_intermediates; i < m_num_tracks; ++i ) 
             required_unique_vertexID += m_daughter_charge[i]*particleMasses_evtReco.find( m_daughter_name[i].c_str() )->second.second ;

          vector<vector<int>> goodTracksThatMeet_withoutIntermediates;
      	  if (num_remaining_tracks > 1) goodTracksThatMeet_withoutIntermediates = findTwoProngs( daughterParticles, goodTrackIndex_withoutIntermediates, num_remaining_tracks );
          for ( int p = 3; p <= num_remaining_tracks; ++p) goodTracksThatMeet_withoutIntermediates = findNProngs( daughterParticles, goodTrackIndex_withoutIntermediates, 
                                                                                                                  goodTracksThatMeet_withoutIntermediates, num_remaining_tracks, p);
      
         vector<vector<string>> uniqueCombinations = findUniqueDaughterCombinations( num_tracks_used_by_intermediates , m_num_tracks ); //Unique comb of remaining trackIDs
         vector<string> v_intermediate_name(m_intermediate_name, m_intermediate_name + m_num_intermediate_states);

         vector<vector<int>> listOfTracksToAppend = appendTracksToIntermediates( motherDecayProducts, daughterParticles, goodTrackIndex_withoutIntermediates, num_remaining_tracks);
         for (unsigned int n_names = 0; n_names < uniqueCombinations.size(); ++n_names)
             uniqueCombinations[n_names].insert(begin(uniqueCombinations[n_names]), begin(v_intermediate_name), end(v_intermediate_name));

         for (unsigned int n_tracks = 0; n_tracks < listOfTracksToAppend.size(); ++n_tracks) 
         {
           for (int n_trackID = 0; n_trackID < num_remaining_tracks; ++n_trackID) 
           {
             int mDP_trackElem = m_num_intermediate_states + n_trackID;
             int dP_trackElem = listOfTracksToAppend[n_tracks][n_trackID];
             motherDecayProducts[mDP_trackElem] = daughterParticles[dP_trackElem];
           }
           for (unsigned int n_names = 0; n_names < uniqueCombinations.size(); ++n_names) 
           {
             tie( candidate, isGood ) = getCombination( motherDecayProducts, &uniqueCombinations[n_names][0], primaryVertices[i_pv], m_constrain_to_vertex, false, 0, num_mother_decay_products, m_constrain_int_mass, required_unique_vertexID );
             if (isGood)
             {
               goodCandidates.push_back( candidate );
               if (m_constrain_to_vertex) goodVertex.push_back( primaryVertices[i_pv] );
               for (int k = 0; k < m_num_intermediate_states; ++k) goodIntermediates[ k ].push_back(motherDecayProducts[ k ] );
               for ( int k = 0; k < num_tracks_used_by_intermediates; ++k ) goodDaughters[ k ].push_back( finalTracks[ k ] );
               for ( int k = 0; k < num_remaining_tracks; ++k ) goodDaughters[ k + num_tracks_used_by_intermediates ].push_back( motherDecayProducts[ k + m_num_intermediate_states ] ); 
             }
           }
         }
       }

       if (isGood && num_remaining_tracks == 0)
       {
         goodCandidates.push_back( candidate );
         if (m_constrain_to_vertex) goodVertex.push_back( primaryVertices[i_pv] );
         for (int k = 0; k < m_num_intermediate_states; ++k) goodIntermediates[ k ].push_back(motherDecayProducts[ k ] );
         for ( int k = 0; k < m_num_tracks; ++k ) goodDaughters[ k ].push_back( finalTracks[ k ] );

       }    
   } //Close PVs
   } //Close forth intermediate
   } //Close third intermediate
   } //Close second intermediate
   } //Close first intermediate

     if ( goodCandidates.size() != 0 )
     {
       KFParticle smallestMassError = goodCandidates[0];
       int bestCombinationIndex = 0;
       for ( unsigned int i = 0; i < goodCandidates.size(); ++i)
       {
         if ( goodCandidates[i].GetErrMass() < smallestMassError.GetErrMass() ) { smallestMassError = goodCandidates[i]; bestCombinationIndex = i; }
       }
       selectedMother.push_back( goodCandidates[ bestCombinationIndex ] );
       if ( m_constrain_to_vertex ) selectedVertex.push_back( goodVertex[ bestCombinationIndex ] );
       vector<KFParticle> intermediates;
       for ( int i = 0; i < m_num_intermediate_states; ++i ) intermediates.push_back( goodIntermediates[i][ bestCombinationIndex ] );
       selectedIntermediates.push_back( intermediates );
       vector<KFParticle> particles;
       for ( int i = 0; i < m_num_tracks; ++i ) particles.push_back( goodDaughters[i][ bestCombinationIndex ] );
       selectedDaughters.push_back( particles );
     }
     goodCandidates.clear();
     goodVertex.clear();
     for (int j = 0; j < m_num_intermediate_states; ++j ) goodIntermediates[j].clear();
     for (int j = 0; j < m_num_tracks; ++j ) goodDaughters[j].clear();
	
}


void KFParticle_eventReconstruction::getCandidateDecay(vector<KFParticle>& selectedMother,
                                         vector<KFParticle>& selectedVertex,
                                         vector<vector<KFParticle>>& selectedDaughters,
                                         vector<KFParticle> daughterParticles,
                                         vector<vector<int>> goodTracksThatMeet,
                                         vector<KFParticle> primaryVertices,
                                         int n_track_start, int n_track_stop, 
                                         bool isIntermediate, int intermediateNumber, bool constrainMass)
{
  int nTracks =  n_track_stop - n_track_start;
  vector<vector<string>> uniqueCombinations = findUniqueDaughterCombinations( n_track_start, n_track_stop );
  vector<KFParticle> goodCandidates, goodVertex, goodDaughters[ nTracks ];
  KFParticle candidate;
  bool isGood;
  bool fixToPV = m_constrain_to_vertex && !isIntermediate;

  float required_unique_vertexID = 0;
  for ( int i = n_track_start; i < n_track_stop; ++i ) required_unique_vertexID += m_daughter_charge[i]*particleMasses_evtReco.find( m_daughter_name[i].c_str() )->second.second ;

  for ( unsigned int i_comb = 0; i_comb < goodTracksThatMeet.size(); ++i_comb ) //Loop over all good track combinations
  {
    KFParticle daughterTracks[ nTracks ];

    for ( int i_track = 0; i_track < nTracks; ++i_track ) { daughterTracks[ i_track ] = daughterParticles[ goodTracksThatMeet[ i_comb ][ i_track ] ]; } //Build array of the good tracks in that combination

      for ( unsigned int i_uc = 0; i_uc < uniqueCombinations.size(); ++i_uc ) //Loop over unique track PID assignments
      {
        for ( unsigned int i_pv = 0; i_pv < primaryVertices.size(); ++i_pv ) //Loop over all PVs in the event
        {
          string *names = &uniqueCombinations[ i_uc ][0];
          tie( candidate, isGood ) = getCombination( daughterTracks, names, primaryVertices[ i_pv ], m_constrain_to_vertex, 
                                                          isIntermediate, intermediateNumber, nTracks, constrainMass, required_unique_vertexID );

          if (isGood)
          {
            goodCandidates.push_back( candidate );
            goodVertex.push_back( primaryVertices[ i_pv ] );
            for ( int i = 0; i < nTracks; ++i )
            {
                KFParticle intParticle;
                intParticle.Create( daughterTracks[ i ].Parameters(),
                                    daughterTracks[ i ].CovarianceMatrix(),
                            (Int_t) daughterTracks[ i ].GetQ(),
                                particleMasses_evtReco.find( names[ i ].c_str() )->second.second );
                intParticle.SetId( daughterTracks[ i ].Id() );
                intParticle.SetPDG( daughterTracks[ i ].GetQ()*particleMasses_evtReco.find( names[ i ].c_str() )->second.first );
                goodDaughters[ i ].push_back(intParticle);
            }
          }

        }
      }

     if ( goodCandidates.size() != 0 )
     {
       KFParticle smallestMassError = goodCandidates[0];
       int bestCombinationIndex = 0;
       for ( unsigned int i = 0; i < goodCandidates.size(); ++i)
       {
         if ( goodCandidates[i].GetErrMass() < smallestMassError.GetErrMass() ) { smallestMassError = goodCandidates[i]; bestCombinationIndex = i; }
       }
       selectedMother.push_back( goodCandidates[ bestCombinationIndex ] );
       if ( fixToPV ) selectedVertex.push_back( goodVertex[ bestCombinationIndex ] );
       vector<KFParticle> particles;
       for ( int i = 0; i < nTracks; ++i ) particles.push_back( goodDaughters[i][ bestCombinationIndex ] );
       selectedDaughters.push_back( particles );
     }
     goodCandidates.clear();
     goodVertex.clear();
     for (int j = 0; j < nTracks; ++j ) goodDaughters[j].clear();
  }
}
