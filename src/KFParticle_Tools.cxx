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

//KFParticle stuff
#include "KFParticle.h"
#include "KFPTrack.h"
#include "KFPVertex.h"
#include "KFVertex.h"
#include "KFParticleDatabase.h"

/// Create necessary objects 
//KFParticleDatabase kfpDatabase;
KFParticle_particleList kfp_particleList;

//Particle masses are in GeV
std::map<std::string, float> particleMasses = kfp_particleList.getParticleList(); 

/// KFParticle constructor
KFParticle_Tools::KFParticle_Tools():
    m_has_intermediates( false ),
    m_num_tracks( 2 ),
    m_daughter_name{"pion", "pion", "pion", "pion"},
    m_daughter_charge{ 1, -1, 1, -1}, 
    m_min_mass( 0 ),
    m_max_mass( 1e1 ),
    m_min_lifetime( 0 ),
    m_max_lifetime( 1e1 ),
    m_track_pt( 0.25 ),
    m_track_ptchi2( FLT_MAX ),
    m_track_ipchi2( 10. ),
    m_track_chi2ndof( 4. ),
    m_comb_DCA( 0.05 ),
    m_vertex_chi2ndof( 15. ),
    m_fdchi2( 50. ),
    m_dira_min( 0.95 ),
    m_dira_max( 1.01 ),
    m_mother_pt( 0. ),
    m_mother_ipchi2( FLT_MAX ),
    m_constrain_to_vertex( true ),
    m_constrain_int_mass( false ), 
    m_get_charge_conjugate( true ), 
    m_dst_vertexmap(),
    m_dst_trackmap(),
    m_dst_vertex(),
    m_dst_track()
{}

KFParticle_Tools::~KFParticle_Tools(){} /// KFParticle destructor


void KFParticle_Tools::createDecay( PHCompositeNode *topNode, std::vector<KFParticle>& selectedMother, std::vector<KFParticle>& selectedVertex,
                                                              std::vector<std::vector<KFParticle>>& selectedDaughters,
                                                              std::vector<std::vector<KFParticle>>& selectedIntermediates,
                                                              int& nPVs, int& multiplicity)
{
  KFParticle::SetField( -1.5e0 );

  std::vector<KFPVertex> primaryVertices = makeAllPrimaryVertices( topNode );
  std::vector<KFParticle> daughterParticles = makeAllDaughterParticles( topNode );

  nPVs = primaryVertices.size();
  multiplicity = daughterParticles.size();

  std::vector<int> goodTrackIndex = findAllGoodTracks( daughterParticles, primaryVertices );

  if ( !m_has_intermediates ) buildBasicChain(selectedMother, selectedVertex,  selectedDaughters, daughterParticles, goodTrackIndex, primaryVertices);
  else buildChain(selectedMother, selectedVertex,  selectedDaughters, selectedIntermediates, daughterParticles, goodTrackIndex, primaryVertices);
 }

/*
 *  This function is used to build a basic n-body decay without any intermediate particles such as D's or J/psi's
 */
void KFParticle_Tools::buildBasicChain(std::vector<KFParticle>& selectedMother, 
	                               std::vector<KFParticle>& selectedVertex, 
	                               std::vector<std::vector<KFParticle>>& selectedDaughters, 
                                       std::vector<KFParticle> daughterParticles,
                                       std::vector<int> goodTrackIndex,
	                               std::vector<KFPVertex> primaryVertices)
{
  std::vector<std::vector<int>> goodTracksThatMeet = findTwoProngs( daughterParticles, goodTrackIndex, m_num_tracks );
  for ( int p = 3; p < m_num_tracks + 1; ++p) goodTracksThatMeet = findNProngs( daughterParticles, goodTrackIndex, goodTracksThatMeet, m_num_tracks, p);

  getCandidateDecay(selectedMother, selectedVertex, selectedDaughters, daughterParticles,
                    goodTracksThatMeet, primaryVertices, 0, m_num_tracks, false, 0, true);

}

/*
 *  This function is used to build a more complicated decay with intermediate particles such as D's or J/psi's
 */
void KFParticle_Tools::buildChain(std::vector<KFParticle>& selectedMother,
                                  std::vector<KFParticle>& selectedVertex,
                                  std::vector<std::vector<KFParticle>>& selectedDaughters,
                                  std::vector<std::vector<KFParticle>>& selectedIntermediates,
                                  std::vector<KFParticle> daughterParticles,
                                  std::vector<int> goodTrackIndex,
                                  std::vector<KFPVertex> primaryVertices)
{
    int track_start = 0;
    int track_stop = m_num_tracks_from_intermediate[0];

    std::vector<KFParticle> goodCandidates, goodVertex, goodDaughters[ m_num_tracks ], goodIntermediates[ m_num_intermediate_states ];
    std::vector<KFParticle> potentialIntermediates[ m_num_intermediate_states ];
    std::vector<std::vector<KFParticle>> potentialDaughters[ m_num_intermediate_states ];

    for (int i = 0; i < m_num_intermediate_states; ++i) 
    {
        std::vector<KFParticle> vertices;

    	std::vector<std::vector<int>> goodTracksThatMeet = findTwoProngs( daughterParticles, goodTrackIndex, m_num_tracks_from_intermediate[i] );
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
//Need to ensure all loops are run at least once or else internal code wont run
        unsigned int c, d; c = d = 0;
        unsigned int matchIterators[4] = { a, b, c, d };

        int num_mother_decay_products = m_num_intermediate_states + num_remaining_tracks;
        KFParticle motherDecayProducts[ num_mother_decay_products ];
        std::vector<KFParticle> finalTracks = potentialDaughters[0][a]; 
        for (int i = 0; i < m_num_intermediate_states; ++i) motherDecayProducts[i] = potentialIntermediates[i][matchIterators[i]]; 
        for (int j = 1; j < m_num_intermediate_states; ++j) finalTracks.insert(finalTracks.end(), potentialDaughters[j][matchIterators[j]].begin(), potentialDaughters[j][matchIterators[j]].end());

        // If there are daughter tracks coming from the mother not an intermediate, need to ensure that the intermeditate decay tracks aren't used again
        std::vector<int> goodTrackIndex_withoutIntermediates = goodTrackIndex;
        for (int m = 0; m < m_num_intermediate_states; ++m) 
        { 
          int trackID_to_remove =  finalTracks[m].Id();
          goodTrackIndex_withoutIntermediates.erase(std::remove(goodTrackIndex_withoutIntermediates.begin(),
                                                                goodTrackIndex_withoutIntermediates.end(), trackID_to_remove),
                                                                goodTrackIndex_withoutIntermediates.end()); 
        }

        float required_unique_vertexID = 0;
        for ( int n = 0; n < m_num_intermediate_states; ++n ) 
           required_unique_vertexID += m_intermediate_charge[n]*particleMasses.find( m_intermediate_name[n].c_str() )->second ;

        if (num_remaining_tracks == 0) std::tie( candidate, isGood ) = getCombination( motherDecayProducts, m_intermediate_name, primaryVertices[i_pv], 
                                                                       true, false, 0, num_mother_decay_products, m_constrain_int_mass, required_unique_vertexID );
        else//Build n-prong from remaining tracks if needed
        {
          for ( int i = num_tracks_used_by_intermediates; i < m_num_tracks; ++i ) 
             required_unique_vertexID += m_daughter_charge[i]*particleMasses.find( m_daughter_name[i].c_str() )->second ;

          std::vector<std::vector<int>> goodTracksThatMeet_withoutIntermediates;
      	  if (num_remaining_tracks > 1) goodTracksThatMeet_withoutIntermediates = findTwoProngs( daughterParticles, goodTrackIndex_withoutIntermediates, num_remaining_tracks );
          for ( int p = 3; p <= num_remaining_tracks; ++p) goodTracksThatMeet_withoutIntermediates = findNProngs( daughterParticles, goodTrackIndex_withoutIntermediates, 
                                                                                                                  goodTracksThatMeet_withoutIntermediates, num_remaining_tracks, p);
      
         std::vector<std::vector<std::string>> uniqueCombinations = findUniqueDaughterCombinations( num_tracks_used_by_intermediates , m_num_tracks ); //Unique comb of remaining trackIDs
         std::vector<std::string> v_intermediate_name(m_intermediate_name, m_intermediate_name + m_num_intermediate_states);

         std::vector<std::vector<int>> listOfTracksToAppend = appendTracksToIntermediates( motherDecayProducts, daughterParticles, goodTrackIndex_withoutIntermediates, num_remaining_tracks);
         for (unsigned int n_names = 0; n_names < uniqueCombinations.size(); ++n_names)
             uniqueCombinations[n_names].insert(std::begin(uniqueCombinations[n_names]), std::begin(v_intermediate_name), std::end(v_intermediate_name));

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
             std::tie( candidate, isGood ) = getCombination( motherDecayProducts, &uniqueCombinations[n_names][0], primaryVertices[i_pv], true, false, 0, num_mother_decay_products, m_constrain_int_mass, required_unique_vertexID );
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
       std::vector<KFParticle> intermediates;
       for ( int i = 0; i < m_num_intermediate_states; ++i ) intermediates.push_back( goodIntermediates[i][ bestCombinationIndex ] );
       selectedIntermediates.push_back( intermediates );
       std::vector<KFParticle> particles;
       for ( int i = 0; i < m_num_tracks; ++i ) particles.push_back( goodDaughters[i][ bestCombinationIndex ] );
       selectedDaughters.push_back( particles );
     }
     goodCandidates.clear();
     goodVertex.clear();
     for (int j = 0; j < m_num_intermediate_states; ++j ) goodIntermediates[j].clear();
     for (int j = 0; j < m_num_tracks; ++j ) goodDaughters[j].clear();
	
}


void KFParticle_Tools::getCandidateDecay(std::vector<KFParticle>& selectedMother,
                                         std::vector<KFParticle>& selectedVertex,
                                         std::vector<std::vector<KFParticle>>& selectedDaughters,
                                         std::vector<KFParticle> daughterParticles,
                                         std::vector<std::vector<int>> goodTracksThatMeet,
                                         std::vector<KFPVertex> primaryVertices,
                                         int n_track_start, int n_track_stop, 
                                         bool isIntermediate, int intermediateNumber, bool constrainMass)
{
  int nTracks =  n_track_stop - n_track_start;
  std::vector<std::vector<std::string>> uniqueCombinations = findUniqueDaughterCombinations( n_track_start, n_track_stop );
  std::vector<KFParticle> goodCandidates, goodVertex, goodDaughters[ nTracks ];
  std::vector<KFParticle> goodDaughter_vector;
  KFParticle candidate;
  bool isGood;
  bool fixToPV = m_constrain_to_vertex && !isIntermediate;

  float required_unique_vertexID = 0;
  for ( int i = n_track_start; i < n_track_stop; ++i ) required_unique_vertexID += m_daughter_charge[i]*particleMasses.find( m_daughter_name[i].c_str() )->second ;

  for ( unsigned int i_comb = 0; i_comb < goodTracksThatMeet.size(); ++i_comb ) //Loop over all good track combinations
  {
    KFParticle daughterTracks[ nTracks ];

    for ( int i_track = 0; i_track < nTracks; ++i_track ) { daughterTracks[ i_track ] = daughterParticles[ goodTracksThatMeet[ i_comb ][ i_track ] ]; } //Build array of the good tracks in that combination

      for ( unsigned int i_uc = 0; i_uc < uniqueCombinations.size(); ++i_uc ) //Loop over unique track PID assignments
      {
        for ( unsigned int i_pv = 0; i_pv < primaryVertices.size(); ++i_pv ) //Loop over all PVs in the event
        {
          std::string *names = &uniqueCombinations[ i_uc ][0];
          std::tie( candidate, isGood ) = getCombination( daughterTracks, names, primaryVertices[ i_pv ], m_constrain_to_vertex, 
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
                                    particleMasses.find( names[ i ].c_str() )->second );
                intParticle.SetId(  daughterTracks[ i ].Id() );
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
       std::vector<KFParticle> particles;
       for ( int i = 0; i < nTracks; ++i ) particles.push_back( goodDaughters[i][ bestCombinationIndex ] );
       selectedDaughters.push_back( particles );
     }
     goodCandidates.clear();
     goodVertex.clear();
     for (int j = 0; j < nTracks; ++j ) goodDaughters[j].clear();
  }
}

KFPVertex KFParticle_Tools::makeVertex( PHCompositeNode *topNode )
{ 
   KFPVertex kfp_vertex;

   kfp_vertex.SetXYZ( m_dst_vertex->get_x(),
                      m_dst_vertex->get_y(),
                      m_dst_vertex->get_z() );

   kfp_vertex.SetCovarianceMatrix( m_dst_vertex->get_error( 0,0 ),
                                   0,//m_dst_vertex->get_error( 1,0 ),
                                   m_dst_vertex->get_error( 1,1 ),
                                   0,//m_dst_vertex->get_error( 2,0 ),
                                   0,//m_dst_vertex->get_error( 2,1 ),
                                   m_dst_vertex->get_error( 2,2 ) );

   kfp_vertex.SetNDF( m_dst_vertex->get_ndof() );

   kfp_vertex.SetChi2( m_dst_vertex->get_chisq() );

   return kfp_vertex;

}


std::vector<KFPVertex> KFParticle_Tools::makeAllPrimaryVertices( PHCompositeNode *topNode )
{ 
  std::vector<KFPVertex> primaryVertices;
  m_dst_vertexmap = findNode::getClass<SvtxVertexMap>( topNode, "SvtxVertexMap" );
    
  for ( SvtxVertexMap::ConstIter iter = m_dst_vertexmap->begin(); iter != m_dst_vertexmap->end(); ++iter )
  { 
    m_dst_vertex = iter->second;
    primaryVertices.push_back( makeVertex( topNode ) );
  }
  
  return primaryVertices;
}


KFParticle KFParticle_Tools::makeParticle( PHCompositeNode *topNode ) ///Return a KFPTrack from track vector and covariance matrix. No mass or vertex constraints
{ 
  float f_trackParameters[6] = { m_dst_track->get_x(),
                                 m_dst_track->get_y(),
                                 m_dst_track->get_z(),
                                 m_dst_track->get_px(),
                                 m_dst_track->get_py(),
                                 m_dst_track->get_pz() };

  float f_trackCovariance[21];
  unsigned int iterate = 0;
  for (unsigned int i = 0; i < 6; ++i) 
    for (unsigned int j = 0; j <= i; ++j) 
      { f_trackCovariance[iterate] = m_dst_track->get_error( i,j ); if (i != j) f_trackCovariance[iterate] = 0; ++iterate;}

  KFParticle kfp_particle;
  kfp_particle.Create( f_trackParameters, f_trackCovariance, (Int_t) m_dst_track->get_charge(), -1);
  kfp_particle.NDF() = m_dst_track->get_ndf();
  kfp_particle.Chi2() = m_dst_track->get_chisq();
  kfp_particle.SetId( m_dst_track->get_id() );

  return kfp_particle;
}


std::vector<KFParticle> KFParticle_Tools::makeAllDaughterParticles( PHCompositeNode *topNode )
{ 
  std::vector<KFParticle> daughterParticles;
  unsigned int trackID = 0;
  m_dst_trackmap = findNode::getClass<SvtxTrackMap>( topNode, "SvtxTrackMap" );

  for ( SvtxTrackMap::Iter iter = m_dst_trackmap->begin(); iter != m_dst_trackmap->end(); ++iter )
  { 
     m_dst_track = iter->second;
     daughterParticles.push_back( makeParticle( topNode ) ); ///Turn all dst tracks in KFP tracks
     daughterParticles[trackID].SetId( iter->first );
     ++trackID;
  }

  return daughterParticles;
}


bool KFParticle_Tools::isGoodTrack( KFParticle particle, std::vector<KFPVertex> primaryVertices )
{ 
  bool goodTrack = false;

  float pt = particle.GetPt();
  float pterr = particle.GetErrPt();
  float ptchi2 = std::pow( pterr/pt, 2);
  float trackchi2ndof = particle.GetChi2()/particle.GetNDF();
  std::vector<float> ipchi2;

  for ( unsigned int i_verts = 0; i_verts < primaryVertices.size(); ++i_verts )
  { 
    KFParticle kfparticlevertex( primaryVertices[ i_verts ] );
    ipchi2.push_back( particle.GetDeviationFromVertex( kfparticlevertex ) );
  }

  auto minmax_ipchi2 = std::minmax_element( ipchi2.begin(), ipchi2.end() ); //Order the IP chi2 from small to large
  float min_ipchi2 = *minmax_ipchi2.first;

  if ( pt > m_track_pt && ptchi2 < m_track_ptchi2 && min_ipchi2 > m_track_ipchi2 && trackchi2ndof < m_track_chi2ndof ) goodTrack = true;

  return goodTrack;
}


std::vector<int> KFParticle_Tools::findAllGoodTracks( std::vector<KFParticle> daughterParticles, std::vector<KFPVertex> primaryVertices ) 
{ 
  std::vector<int> goodTrackIndex;

  for ( unsigned int i_parts = 0; i_parts < daughterParticles.size(); ++i_parts )
  { 
    if ( isGoodTrack( daughterParticles[ i_parts ], primaryVertices ) ) goodTrackIndex.push_back( i_parts );
  }
  
  removeDuplicates( goodTrackIndex );

  return goodTrackIndex;
}


std::vector<std::vector<int>> KFParticle_Tools::findTwoProngs( std::vector<KFParticle> daughterParticles, std::vector<int> goodTrackIndex, int nTracks )
{ 
  std::vector<std::vector<int>> goodTracksThatMeet;
  KFVertex twoParticleVertex;  

  for ( std::vector <int>::iterator i_it = goodTrackIndex.begin(); i_it != goodTrackIndex.end(); ++i_it )
  { 
    for ( std::vector <int>::iterator j_it = goodTrackIndex.begin(); j_it != goodTrackIndex.end(); ++j_it )
    { 
      if( i_it < j_it )
      {
        if( daughterParticles[*i_it ].GetDistanceFromParticle( daughterParticles[*j_it ] ) < m_comb_DCA ) 
        { 
          twoParticleVertex += daughterParticles[*i_it ];
          twoParticleVertex += daughterParticles[*j_it ];
          float vertexchi2ndof = twoParticleVertex.GetChi2()/twoParticleVertex.GetNDF();
          std::vector<int> combination = { *i_it, *j_it };
          if ( nTracks == 2 && vertexchi2ndof < m_vertex_chi2ndof ) goodTracksThatMeet.push_back( combination );
          else goodTracksThatMeet.push_back( combination ); 
        }
      }
    }
  }
  
  return goodTracksThatMeet;
}

std::vector<std::vector<int>>  KFParticle_Tools::findNProngs( std::vector<KFParticle> daughterParticles, 
                                                              std::vector<int> goodTrackIndex, 
                                                              std::vector<std::vector<int>> goodTracksThatMeet, 
                                                              int nRequiredTracks, unsigned int nProngs )
{
  unsigned int nGoodProngs = goodTracksThatMeet.size();
  KFVertex particleVertex;

  for ( std::vector <int>::iterator i_it = goodTrackIndex.begin(); i_it != goodTrackIndex.end(); ++i_it )
  {
    for ( unsigned int i_prongs = 0; i_prongs < nGoodProngs; ++i_prongs )
    {
      bool trackNotUsedAlready = true;
      for (unsigned int i_trackCheck = 0; i_trackCheck < nProngs - 1; ++i_trackCheck)
  	if( *i_it == goodTracksThatMeet[ i_prongs ][ i_trackCheck ] ) trackNotUsedAlready = false;
      if( trackNotUsedAlready )
      {
        bool dcaMet = 1;
        for ( unsigned int i = 0; i < nProngs - 1; ++i ) 
        {
          if( daughterParticles[ *i_it ].GetDistanceFromParticle( daughterParticles[ goodTracksThatMeet[ i_prongs ][ i ] ] ) > m_comb_DCA ) 
            { dcaMet = 0; }
            //{ dcaMet *= 0; }
        }

        if( dcaMet )
        {
            particleVertex += daughterParticles[*i_it ];
            std::vector<int> combination; combination.push_back( *i_it );
            for (unsigned int i = 0; i < nProngs - 1; ++i) 
            {
              particleVertex += daughterParticles[ goodTracksThatMeet[ i_prongs ][ i ] ];
              combination.push_back( goodTracksThatMeet[ i_prongs ][ i ] );
            }
            float vertexchi2ndof = particleVertex.GetChi2()/particleVertex.GetNDF();
            if ( (unsigned int) nRequiredTracks == nProngs &&  vertexchi2ndof < m_vertex_chi2ndof ) goodTracksThatMeet.push_back( combination );
            else goodTracksThatMeet.push_back( combination );
        }
      }
    }
  }

  goodTracksThatMeet.erase(  goodTracksThatMeet.begin(), goodTracksThatMeet.begin() + nGoodProngs );
  for ( unsigned int i = 0; i < goodTracksThatMeet.size(); ++i ) std::sort( goodTracksThatMeet[i].begin(), goodTracksThatMeet[i].end() ); 
 removeDuplicates( goodTracksThatMeet );

  return goodTracksThatMeet;
}


std::vector<std::vector<int>> KFParticle_Tools::appendTracksToIntermediates( KFParticle intermediateResonances[], std::vector<KFParticle> daughterParticles, std::vector<int> goodTrackIndex, int num_remaining_tracks)
{
  std::vector<std::vector<int>> goodTracksThatMeet, goodTracksThatMeetIntermediates;//, vectorOfGoodTracks;
  if (num_remaining_tracks == 1)
  {
    for ( std::vector<int>::iterator i_it = goodTrackIndex.begin(); i_it != goodTrackIndex.end(); ++i_it )
    {
      std::vector<KFParticle> v_intermediateResonances(intermediateResonances, intermediateResonances + m_num_intermediate_states);
      std::vector<std::vector<int>> dummyTrackList;
      std::vector<int> dummyTrackID; //I already have the track ids stored in goodTracksThatMeet[i] 
      v_intermediateResonances.insert(std::end(v_intermediateResonances), daughterParticles[ *i_it ]);
      for ( unsigned int k = 0; k < v_intermediateResonances.size(); ++k) dummyTrackID.push_back(k);
      dummyTrackList = findTwoProngs( v_intermediateResonances, dummyTrackID, (int) v_intermediateResonances.size() );
      if (v_intermediateResonances.size() > 2)
      {
      for ( unsigned int p = 3; p <= v_intermediateResonances.size(); ++p) dummyTrackList = findNProngs( v_intermediateResonances,
                                                                                                         dummyTrackID, dummyTrackList,
                                                                                                         (int) v_intermediateResonances.size(), (int) p);
      }

      if (dummyTrackList.size() != 0) { std::vector<int> goodTrack{ *i_it };  goodTracksThatMeetIntermediates.push_back( goodTrack ); }

    }
  }
  else
  { 
    goodTracksThatMeet = findTwoProngs( daughterParticles, goodTrackIndex, num_remaining_tracks );
    for ( int p = 3; p <= num_remaining_tracks; ++p) goodTracksThatMeet = findNProngs( daughterParticles,
                                                                                       goodTrackIndex,
                                                                                       goodTracksThatMeet,
                                                                                       num_remaining_tracks, p);
    for (unsigned int i = 0; i < goodTracksThatMeet.size(); ++i)
    {
      std::vector<KFParticle> v_intermediateResonances(intermediateResonances, intermediateResonances + m_num_intermediate_states);
      std::vector<std::vector<int>> dummyTrackList;
      std::vector<int> dummyTrackID; //I already have the track ids stored in goodTracksThatMeet[i] 
      for ( unsigned int j = 0; j < goodTracksThatMeet[i].size(); ++j)    v_intermediateResonances.push_back(daughterParticles[goodTracksThatMeet[i][j]]);
      for ( unsigned int k = 0; k < v_intermediateResonances.size(); ++k) dummyTrackID.push_back(k);
      dummyTrackList = findTwoProngs( v_intermediateResonances, dummyTrackID, (int) v_intermediateResonances.size() );
      for ( unsigned int p = 3; p <= v_intermediateResonances.size(); ++p) dummyTrackList = findNProngs( v_intermediateResonances,
                                                                                                         dummyTrackID, dummyTrackList,
                                                                                                         (int) v_intermediateResonances.size(), (int) p);

      if (dummyTrackList.size() != 0) goodTracksThatMeetIntermediates.push_back(goodTracksThatMeet[i]); 
    }
  }

  return goodTracksThatMeetIntermediates;
}


float KFParticle_Tools::eventDIRA( KFParticle particle, KFPVertex vertex )
{ 

  TMatrixD flightVector( 3, 1 ), momVector( 3, 1 );
  flightVector( 0,0 ) = particle.GetX() - vertex.GetX();
  flightVector( 1,0 ) = particle.GetY() - vertex.GetY();
  flightVector( 2,0 ) = particle.GetZ() - vertex.GetZ();

  momVector( 0, 0 ) = particle.GetPx();
  momVector( 1, 0 ) = particle.GetPy();
  momVector( 2, 0 ) = particle.GetPz();

  TMatrixD momDotFD( 1,1 ); //Calculate momentum dot flight distance
  momDotFD = TMatrixD( momVector, TMatrixD::kTransposeMult, flightVector );
  float f_momDotFD = momDotFD( 0,0 );
  
  TMatrixD sizeOfMom( 1,1 ); //Calculates the size of the momentum vector
  sizeOfMom = TMatrixD( momVector, TMatrixD::kTransposeMult, momVector );
  float f_sizeOfMom = std::sqrt( sizeOfMom( 0,0 ) );
  
  TMatrixD sizeOfFD( 1,1 ); //Calculates the size of the flight distance vector
  sizeOfFD = TMatrixD( flightVector, TMatrixD::kTransposeMult, flightVector );
  float f_sizeOfFD = std::sqrt( sizeOfFD( 0,0 ) );
  
  return f_momDotFD/( f_sizeOfMom*f_sizeOfFD ); //returns the DIRA
}


float KFParticle_Tools::flightDistanceChi2( KFParticle particle, KFPVertex vertex )
{ 
  TMatrixD flightVector( 3, 1 ), flightDistanceCovariance( 3, 3 );
  
  KFParticle kfp_vertex( vertex );

  flightVector( 0,0 ) = particle.GetX() - kfp_vertex.GetX();
  flightVector( 1,0 ) = particle.GetY() - kfp_vertex.GetY();
  flightVector( 2,0 ) = particle.GetZ() - kfp_vertex.GetZ();  

  for ( int i = 0; i < 3; i++ ) { for ( int j = 0; j < 3; j++ ) { flightDistanceCovariance( i, j ) = particle.GetCovariance( i, j ) + kfp_vertex.GetCovariance( i, j ); } }

  TMatrixD anInverseMatrix( 3,3 );
  anInverseMatrix = flightDistanceCovariance.Invert();
  TMatrixD m_chi2Value( 1,1 );
  m_chi2Value = TMatrixD( flightVector, TMatrixD::kTransposeMult, anInverseMatrix*flightVector );
  return m_chi2Value( 0,0 );
}


std::tuple<KFParticle, bool> KFParticle_Tools::buildMother( KFParticle vDaughters[], std::string daughterOrder[], bool isIntermediate, int intermediateNumber, int nTracks, bool constrainMass, float required_vertexID )
{
    KFParticle mother, inputTracks[ nTracks ];

    mother.SetConstructMethod( 2 );

    bool daughterMassCheck = true; 
    float unique_vertexID = 0;

    for ( int i = 0; i < nTracks; ++i )
    {
      float daughterMass = constrainMass ? particleMasses.find( daughterOrder[ i ].c_str() )->second : vDaughters[ i ].GetMass(); 
      inputTracks[ i ].Create( vDaughters[ i ].Parameters(),
                               vDaughters[ i ].CovarianceMatrix(),
                       (Int_t) vDaughters[ i ].GetQ(),
                               daughterMass );
      mother.AddDaughter( inputTracks[ i ] );
      if ( inputTracks[ i ].GetMass() == 0 ) daughterMassCheck = false;
      unique_vertexID += vDaughters[ i ].GetQ()*particleMasses.find( daughterOrder[ i ].c_str() )->second;
    }
   bool chargeCheck;
   if (m_get_charge_conjugate) chargeCheck = std::abs(unique_vertexID) == std::abs(required_vertexID) ? 1 : 0;
   else chargeCheck = unique_vertexID == required_vertexID ? 1 : 0;

   for ( int j = 0; j < nTracks; ++j ) inputTracks[ j ].SetProductionVertex( mother );

   float calculated_mass, calculated_mass_err;
   mother.GetMass( calculated_mass, calculated_mass_err );
   float calculated_pt = mother.GetPt();
   
   float min_mass = isIntermediate ? m_intermediate_mass_range[intermediateNumber].first : m_min_mass;
   float max_mass = isIntermediate ? m_intermediate_mass_range[intermediateNumber].second : m_max_mass;
   float min_pt = isIntermediate ? m_intermediate_min_pt[intermediateNumber] : m_mother_pt;

   bool goodCandidate = false;
   if ( calculated_mass > min_mass && calculated_mass < max_mass && 
        calculated_pt > min_pt && 
        daughterMassCheck && chargeCheck ) 
        goodCandidate = true;

   return std::make_tuple( mother, goodCandidate );
}


void KFParticle_Tools::constrainToVertex( KFParticle& particle, bool& goodCandidate, KFPVertex& vertex )
{
  KFParticle prod_vertex( vertex );
  //prod_vertex.AddDaughter( particle );
  //particle.SetProductionVertex( prod_vertex );
  //particle.SetAtProductionVertex( true );

  float calculated_fdchi2  = flightDistanceChi2( particle, vertex );
  float calculated_dira    = eventDIRA( particle, vertex );
  float calculated_ipchi2  = particle.GetDeviationFromVertex( prod_vertex );
  //float calculated_lifetime, calculated_lifetime_error;
  //mother.GetLifeTime( calculated_lifetime, calculated_lifetime_error );
  goodCandidate = false;
  if ( calculated_fdchi2 > m_fdchi2 && calculated_ipchi2 < m_mother_ipchi2 &&
       calculated_dira > m_dira_min && calculated_dira < m_dira_max )  goodCandidate = true;
}


std::tuple<KFParticle, bool> KFParticle_Tools::getCombination( KFParticle vDaughters[], std::string daughterOrder[], KFPVertex vertex, bool constrain_to_vertex, bool isIntermediate, int intermediateNumber, int nTracks, bool constrainMass, float required_vertexID )
{
   KFParticle candidate;
   bool isGoodCandidate;

   std::tie( candidate, isGoodCandidate ) = buildMother( vDaughters, daughterOrder, isIntermediate, intermediateNumber, nTracks, constrainMass, required_vertexID );   
   if ( constrain_to_vertex && isGoodCandidate && !isIntermediate ) constrainToVertex( candidate, isGoodCandidate, vertex );

   return std::make_tuple( candidate, isGoodCandidate );
}


std::vector<std::vector<std::string>> KFParticle_Tools::findUniqueDaughterCombinations( int start, int end )
{
  std::vector<int> vect_permutations;
  std::vector<std::vector<std::string>> uniqueCombinations;
  std::map<int, std::string> daughterMap;
  for( int i = start; i < end; i++)
  {
    daughterMap.insert( std::pair<int, std::string>( i, m_daughter_name[i].c_str() ) );
    vect_permutations.push_back(i);
  }
  int *permutations = &vect_permutations[0];

  do
  {
    std::vector<std::string> combination;
    for( int i = 0; i < (end - start); i++) combination.push_back( daughterMap.find(permutations[i])->second );
    uniqueCombinations.push_back( combination );
  } while ( std::next_permutation( permutations, permutations + vect_permutations.size() ) );      

  removeDuplicates(uniqueCombinations);

  return uniqueCombinations;
}


void KFParticle_Tools::removeDuplicates( std::vector<int> &v )
{ 
  auto end = v.end();
  for ( auto it = v.begin(); it != end; ++it ) 
  { 
    end = std::remove( it + 1, end, *it );
  }
  v.erase( end, v.end() );
}

void KFParticle_Tools::removeDuplicates( std::vector<std::vector<int>> &v )
{ 
  auto end = v.end();
  for ( auto it = v.begin(); it != end; ++it ) 
  { 
    end = std::remove( it + 1, end, *it );
  }
  v.erase( end, v.end() );
}


void KFParticle_Tools::removeDuplicates( std::vector<std::vector<std::string>> &v )
{ 
  auto end = v.end();
  for ( auto it = v.begin(); it != end; ++it ) 
  { 
    end = std::remove( it + 1, end, *it );
  }
  v.erase( end, v.end() );
}
