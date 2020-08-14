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

/** Mass Hypothesis Codes
 ** These are the codes used in the PDG to identify different particles
 ** that have also been coded into KFParticle
 ** PDG Code | Mass Index | Particle
 **         11 |  0 | Electron
 **         13 |  1 | Muon
 **         19 |  1 | Muon
 **        211 |  2 | Pion ( Charged )
 **        321 |  3 | Kaon ( Charged )
 **       2212 |  4 | Proton
 ** 1000010020 |  5 | Deuteron
 ** 1000010030 |  6 | Triton
 ** 1000020030 |  7 | Helium-3
 ** 1000020040 |  8 | Helium-4
 **       3112 |  9 | Sigma ( - )
 **       3222 | 10 | Sigma ( + )
 **       3312 | 11 | Xi
 **       3334 | 12 | Omega
 ** Use PDG codes in your analysis
 ** All other codes return a pion
 **/

#include "KFParticle_Tools.h"

//KFParticle stuff
#include "KFParticle.h"
#include "KFPTrack.h"
#include "KFPVertex.h"
#include "KFVertex.h"
#include "KFParticleDatabase.h"

/// Create necessary objects 
KFParticleDatabase kfpDatabase;

//Particle masses are in GeV
std::map<std::string, float> particleMasses = 
{ 
  { "electron", kfpDatabase.GetMass( 11 ) },
  { "muon",     kfpDatabase.GetMass( 13 ) },
  { "pion",     kfpDatabase.GetMass( 211 ) },
  { "kaon",     kfpDatabase.GetMass( 321 ) },
  { "proton",   kfpDatabase.GetMass( 2212 ) },
  //{ "photon" },       returnPDGMass( 22 ) },
  //{ "KS0",            returnPDGMass( 310 ) },
  //{ "Lambda0",        returnPDGMass( 3122 ) },
  { "pi0",      kfpDatabase.GetPi0Mass() },
  { "D0",       kfpDatabase.GetD0Mass() },
  { "Dplus",    kfpDatabase.GetDPlusMass() },
  { "J/psi",    3.09690 },
  { "phi",      1.019461 } 
 };


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
    m_comb_DCA( 0.2 ),
    m_vertex_chi2ndof( 15. ),
    m_fdchi2( 50. ),
    m_dira_min( 0.95 ),
    m_dira_max( 1.01 ),
    m_mother_pt( 0. ),
    m_constrain_to_vertex( true ), 
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
  std::vector<std::vector<int>> goodTracksThatMeet = findTwoProngs( daughterParticles, goodTrackIndex );
  if ( m_num_tracks > 2 ) goodTracksThatMeet = findThreeProngs( daughterParticles, goodTrackIndex, goodTracksThatMeet );
  if ( m_num_tracks > 3 ) goodTracksThatMeet = findFourProngs( daughterParticles, goodTrackIndex, goodTracksThatMeet );

  getCandidateDecay(selectedMother, selectedVertex, selectedDaughters, daughterParticles,
                    goodTracksThatMeet, primaryVertices, 0, m_num_tracks, false, 0);

}

/*
 *  This function is used to build a more complicated decay with intermediate particles such as D's or J/psi's
 */


void KFParticle_Tools::getCandidateDecay(std::vector<KFParticle>& selectedMother,
                                         std::vector<KFParticle>& selectedVertex,
                                         std::vector<std::vector<KFParticle>>& selectedDaughters,
                                         std::vector<KFParticle> daughterParticles,
                                         std::vector<std::vector<int>> goodTracksThatMeet,
                                         std::vector<KFPVertex> primaryVertices,
                                         int n_track_start, int n_track_stop, 
                                         bool isIntermediate, int intermediateNumber)
{
  std::vector<std::vector<std::string>> uniqueCombinations = findUniqueDaughterCombinations( n_track_start, n_track_stop );
  std::vector<KFParticle> goodCandidates, goodVertex, goodDaughters[ ( n_track_stop - n_track_start ) ];
  KFParticle candidate;
  bool isGood;

  for ( unsigned int i_comb = 0; i_comb < goodTracksThatMeet.size(); ++i_comb ) //Loop over all good track combinations
  {
    KFParticle daughterTracks[ ( n_track_stop - n_track_start ) ];

    for ( int i_track = 0; i_track < ( n_track_stop - n_track_start ); ++i_track ) { daughterTracks[ i_track ] = daughterParticles[ goodTracksThatMeet[ i_comb ][ i_track ] ]; } //Build array of the good tracks in that combination

      for ( unsigned int i_uc = 0; i_uc < uniqueCombinations.size(); ++i_uc ) //Loop over unique track PID assignments
      {
        for ( unsigned int i_pv = 0; i_pv < primaryVertices.size(); ++i_pv ) //Loop over all PVs in the event
        {
          std::string *names = &uniqueCombinations[ i_uc ][0];

          std::tie( candidate, isGood ) = getCombination( daughterTracks, names, primaryVertices[ i_pv ], m_constrain_to_vertex, isIntermediate, intermediateNumber );
          if (isGood)
          {
            goodCandidates.push_back( candidate );
            goodVertex.push_back( primaryVertices[ i_pv ] );
            for ( int i = 0; i < ( n_track_stop - n_track_start ); ++i )
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
        selectedVertex.push_back( goodVertex[ bestCombinationIndex ] );
        std::vector<KFParticle> particles;
        for ( int i = 0; i < ( n_track_stop - n_track_start ); ++i ) particles.push_back( goodDaughters[i][ bestCombinationIndex ] );
        selectedDaughters.push_back( particles );
      }
      goodCandidates.clear();
      goodVertex.clear();
      for (int j = 0; j < ( n_track_stop - n_track_start ); ++j ) goodDaughters[j].clear();
  }
}

KFPVertex KFParticle_Tools::makeVertex( PHCompositeNode *topNode )
{ 
   KFPVertex kfp_vertex;

   kfp_vertex.SetXYZ( m_dst_vertex->get_x(),
                      m_dst_vertex->get_y(),
                      m_dst_vertex->get_z() );

   kfp_vertex.SetCovarianceMatrix( m_dst_vertex->get_error( 0,0 ),
                                   m_dst_vertex->get_error( 1,0 ),
                                   m_dst_vertex->get_error( 1,1 ),
                                   m_dst_vertex->get_error( 2,0 ),
                                   m_dst_vertex->get_error( 2,1 ),
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


KFPTrack KFParticle_Tools::makeTrack( PHCompositeNode *topNode ) ///Return a KFPTrack from track vector and covariance matrix. No mass or vertex constraints
{ 
  KFPTrack kfp_track;

  float f_trackCovariance[21] = { m_dst_track->get_error( 0,0 ), 
                                  m_dst_track->get_error( 1,0 ), 
                                  m_dst_track->get_error( 1,1 ), 
                                  m_dst_track->get_error( 2,0 ), 
                                  m_dst_track->get_error( 2,1 ), 
                                  m_dst_track->get_error( 2,2 ),
                                  m_dst_track->get_error( 3,0 ), 
                                  m_dst_track->get_error( 3,1 ), 
                                  m_dst_track->get_error( 3,2 ), 
                                  m_dst_track->get_error( 3,3 ), 
                                  m_dst_track->get_error( 4,0 ), 
                                  m_dst_track->get_error( 4,1 ),
                                  m_dst_track->get_error( 4,2 ), 
                                  m_dst_track->get_error( 4,3 ),
                                  m_dst_track->get_error( 4,4 ), 
                                  m_dst_track->get_error( 5,0 ), 
                                  m_dst_track->get_error( 5,1 ), 
                                  m_dst_track->get_error( 5,2 ),
                                  m_dst_track->get_error( 5,3 ), 
                                  m_dst_track->get_error( 5,4 ), 
                                  m_dst_track->get_error( 5,5 ) }; 

  kfp_track.SetParameters( m_dst_track->get_x(),
                           m_dst_track->get_y(),
                           m_dst_track->get_z(),
                           m_dst_track->get_px(),
                           m_dst_track->get_py(),
                           m_dst_track->get_pz() );
  kfp_track.SetCovarianceMatrix( f_trackCovariance );
  kfp_track.SetNDF( m_dst_track->get_ndf() );
  kfp_track.SetChi2( m_dst_track->get_chisq() );
  kfp_track.SetCharge( m_dst_track->get_charge() );
  kfp_track.SetId( m_dst_track->get_id() );

  return kfp_track;
}


KFParticle KFParticle_Tools::makeParticle( PHCompositeNode *topNode, int massHypothesis ) ///Return a KFParticle from particle vector and covariance matrix. No mass or vertex constraints
{ 
  KFPTrack kfp_track = makeTrack( topNode );
  KFParticle kfp_particle( kfp_track, massHypothesis );

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
     daughterParticles.push_back( makeParticle( topNode, -1 ) ); ///Turn all dst tracks in KFP tracks
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
    if ( isGoodTrack( daughterParticles[ i_parts ], primaryVertices ) ) { goodTrackIndex.push_back( i_parts ); }
  }
  
  removeDuplicates( goodTrackIndex );

  return goodTrackIndex;
}


std::vector<std::vector<int>> KFParticle_Tools::findTwoProngs( std::vector<KFParticle> daughterParticles, std::vector<int> goodTrackIndex )
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
          if ( m_num_tracks == 2 && vertexchi2ndof < m_vertex_chi2ndof ) goodTracksThatMeet.push_back( combination );
          else goodTracksThatMeet.push_back( combination ); 
        }
      }
    }
  }
  
  return goodTracksThatMeet;
}


std::vector<std::vector<int>> KFParticle_Tools::findThreeProngs( std::vector<KFParticle> daughterParticles, std::vector<int> goodTrackIndex, std::vector<std::vector<int>> goodTracksThatMeet )
{ 
  unsigned int nGoodTwoProngs = goodTracksThatMeet.size();
  KFVertex threeParticleVertex;  

  for ( std::vector <int>::iterator i_it = goodTrackIndex.begin(); i_it != goodTrackIndex.end(); ++i_it )
  { 
    for ( unsigned int i_prongs = 0; i_prongs < nGoodTwoProngs; ++i_prongs )
    { 
      if( *i_it != goodTracksThatMeet[ i_prongs ][ 0 ] && *i_it != goodTracksThatMeet[ i_prongs ][ 1 ] )
      { 
        if( daughterParticles[ *i_it ].GetDistanceFromParticle( daughterParticles[ goodTracksThatMeet[ i_prongs ][ 0 ] ] ) < m_comb_DCA &&
            daughterParticles[ *i_it ].GetDistanceFromParticle( daughterParticles[ goodTracksThatMeet[ i_prongs ][ 1 ] ] ) < m_comb_DCA )
          { 
            threeParticleVertex += daughterParticles[*i_it ];
            threeParticleVertex += daughterParticles[ goodTracksThatMeet[ i_prongs ][ 0 ] ];
            threeParticleVertex += daughterParticles[ goodTracksThatMeet[ i_prongs ][ 1 ] ];
            float vertexchi2ndof = threeParticleVertex.GetChi2()/threeParticleVertex.GetNDF();
            std::vector<int> combination = { goodTracksThatMeet[ i_prongs ][ 0 ], goodTracksThatMeet[ i_prongs ][ 1 ],  *i_it, }; 
            if ( m_num_tracks == 3 &&  vertexchi2ndof < m_vertex_chi2ndof ) goodTracksThatMeet.push_back( combination );
            else goodTracksThatMeet.push_back( combination );
          }
      }
    }
  }

  goodTracksThatMeet.erase(  goodTracksThatMeet.begin(), goodTracksThatMeet.begin() + nGoodTwoProngs );   
  for( unsigned int i = 0; i < goodTracksThatMeet.size(); ++i ) std::sort( goodTracksThatMeet[i].begin(), goodTracksThatMeet[i].end() );
  removeDuplicates( goodTracksThatMeet );

  return goodTracksThatMeet;
}


std::vector<std::vector<int>>  KFParticle_Tools::findFourProngs( std::vector<KFParticle> daughterParticles, std::vector<int> goodTrackIndex, std::vector<std::vector<int>> goodTracksThatMeet )
{ 
  unsigned int nGoodThreeProngs = goodTracksThatMeet.size();
  KFVertex fourParticleVertex;  

  for ( std::vector <int>::iterator i_it = goodTrackIndex.begin(); i_it != goodTrackIndex.end(); ++i_it )
  { 
    for ( unsigned int i_prongs = 0; i_prongs < nGoodThreeProngs; ++i_prongs )
    { 
      if( *i_it != goodTracksThatMeet[ i_prongs ][ 0 ] && 
          *i_it != goodTracksThatMeet[ i_prongs ][ 1 ] &&
          *i_it != goodTracksThatMeet[ i_prongs ][ 2 ] )
      { 
        if( daughterParticles[ *i_it ].GetDistanceFromParticle( daughterParticles[ goodTracksThatMeet[ i_prongs ][ 0 ] ] ) < m_comb_DCA &&
            daughterParticles[ *i_it ].GetDistanceFromParticle( daughterParticles[ goodTracksThatMeet[ i_prongs ][ 1 ] ] ) < m_comb_DCA &&
            daughterParticles[ *i_it ].GetDistanceFromParticle( daughterParticles[ goodTracksThatMeet[ i_prongs ][ 2 ] ] ) < m_comb_DCA )
          { 
            fourParticleVertex += daughterParticles[*i_it ];
            fourParticleVertex += daughterParticles[ goodTracksThatMeet[ i_prongs ][ 0 ] ];
            fourParticleVertex += daughterParticles[ goodTracksThatMeet[ i_prongs ][ 1 ] ];
            fourParticleVertex += daughterParticles[ goodTracksThatMeet[ i_prongs ][ 2 ] ];
            float vertexchi2ndof = fourParticleVertex.GetChi2()/fourParticleVertex.GetNDF();
            std::vector<int> combination = { goodTracksThatMeet[ i_prongs ][ 0 ], 
                                             goodTracksThatMeet[ i_prongs ][ 1 ], 
                                             goodTracksThatMeet[ i_prongs ][ 2 ], 
                                             *i_it, }; 
            if ( vertexchi2ndof < m_vertex_chi2ndof ) goodTracksThatMeet.push_back( combination ); 
          } 
      }
    }
  }

  goodTracksThatMeet.erase(  goodTracksThatMeet.begin(), goodTracksThatMeet.begin() + nGoodThreeProngs );   
  for( unsigned int i = 0; i < goodTracksThatMeet.size(); ++i ) std::sort( goodTracksThatMeet[i].begin(), goodTracksThatMeet[i].end() );
  removeDuplicates( goodTracksThatMeet );

  return goodTracksThatMeet;
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

  flightVector( 0,0 ) = particle.GetX() -  kfp_vertex.GetX();
  flightVector( 1,0 ) = particle.GetY() -  kfp_vertex.GetY();
  flightVector( 2,0 ) = particle.GetZ() -  kfp_vertex.GetZ();  

  for ( int i = 0; i < 3; i++ ) { for ( int j = 0; j < 3; j++ ) { flightDistanceCovariance( i, j ) = particle.GetCovariance( i, j ) + kfp_vertex.GetCovariance( i, j ); } }

  TMatrixD anInverseMatrix( 3,3 );
  anInverseMatrix = flightDistanceCovariance.Invert();
  TMatrixD m_chi2Value( 1,1 );
  m_chi2Value = TMatrixD( flightVector, TMatrixD::kTransposeMult, anInverseMatrix*flightVector );
  return m_chi2Value( 0,0 );
}


std::tuple<KFParticle, bool> KFParticle_Tools::buildMother( KFParticle vDaughters[], std::string daughterOrder[], bool isIntermediate, int intermediateNumber )
{ 
    KFParticle mother, inputTracks[ m_num_tracks ];

    mother.SetConstructMethod( 2 );

    for ( int i = 0; i < m_num_tracks; ++i )
    { 
      inputTracks[ i ].Create( vDaughters[ i ].Parameters(),
                               vDaughters[ i ].CovarianceMatrix(),
                       (Int_t) vDaughters[ i ].GetQ(),
                               particleMasses.find( daughterOrder[ i ].c_str() )->second );
      mother.AddDaughter( inputTracks[ i ] ); 
    }

    for ( int j = 0; j < m_num_tracks; ++j ) inputTracks[ j ].SetProductionVertex( mother );

   float calculated_mass, calculated_mass_err;
   mother.GetMass( calculated_mass, calculated_mass_err );
   float calculated_pt = mother.GetPt();
   
   float min_mass = isIntermediate ? m_intermediate_mass_range[intermediateNumber].first : m_min_mass;
   float max_mass = isIntermediate ? m_intermediate_mass_range[intermediateNumber].second : m_max_mass;
   float min_pt = isIntermediate ? m_intermediate_min_pt[intermediateNumber] : m_mother_pt;

   bool goodCandidate = false;
   if ( calculated_mass > min_mass && calculated_mass < max_mass && 
        calculated_pt > min_pt && chargeChecker( inputTracks ) ) 
        goodCandidate = true;

   return std::make_tuple( mother, goodCandidate );
}


void KFParticle_Tools::constrainToVertex( KFParticle& particle, bool& goodCandidate, KFPVertex& vertex )
{
  KFParticle prod_vertex( vertex );

  prod_vertex.AddDaughter( particle );
  particle.SetProductionVertex( prod_vertex );
  particle.SetAtProductionVertex( true );

  float calculated_fdchi2 = flightDistanceChi2( particle, vertex );
  float calculated_dira = eventDIRA( particle, vertex );
   //float calculated_lifetime, calculated_lifetime_error;
   //mother.GetLifeTime( calculated_lifetime, calculated_lifetime_error );

  goodCandidate = false;
  if ( calculated_fdchi2 > m_fdchi2 && calculated_dira > m_dira_min && calculated_dira < m_dira_max )  goodCandidate = true;

}


std::tuple<KFParticle, bool> KFParticle_Tools::getCombination( KFParticle vDaughters[], std::string daughterOrder[], KFPVertex vertex, bool constrain_to_vertex, bool isIntermediate, int intermediateNumber )
{
   KFParticle candidate;
   bool isGoodCandidate;

   std::tie( candidate, isGoodCandidate ) = buildMother( vDaughters, daughterOrder, isIntermediate, intermediateNumber );   
   if ( constrain_to_vertex && isGoodCandidate ) constrainToVertex( candidate, isGoodCandidate, vertex );

   return std::make_tuple( candidate, isGoodCandidate );
}


bool KFParticle_Tools::chargeChecker( KFParticle vDaughters[] )
{
  double required_unique_ID, vertex_unique_ID;
  required_unique_ID = vertex_unique_ID = 0;

  for ( int i = 0; i < m_num_tracks; i++ ) 
  {
    vertex_unique_ID += vDaughters[i].GetQ()*vDaughters[i].GetMass();
    required_unique_ID += m_daughter_charge[i]*particleMasses.find( m_daughter_name[i].c_str() )->second ;
  }

  bool checkCharge = round(vertex_unique_ID*1e3) == std::abs(round(required_unique_ID*1e3));

  return checkCharge;
}


std::vector<std::vector<std::string>> KFParticle_Tools::findUniqueDaughterCombinations( int start, int end )
{
  std::vector<int> vect_permutations;
  std::vector<std::vector<std::string>> uniqueCombinations;
  std::map<int, std::string> daughterMap;

  for( int i = start; i < end; i++)
  {
    daughterMap.insert( std::pair<int, std::string>( i+1, m_daughter_name[i].c_str() ) );
    vect_permutations.push_back(i+1);
  }
  int *permutations = &vect_permutations[0];

  do
  {
    std::vector<std::string> combination;
    for( int i = 0; i < m_num_tracks; i++) combination.push_back( daughterMap.find(permutations[i])->second );
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


float KFParticle_Tools::returnPDGMass( const int pdgIndex ) ///Return mother masses from KFParticleDatabase ( pdg.lbl.gov/2019/reviews/rpp2019-rev-monte-carlo-numbering.pdf )
{ 
   float mass, width;
   kfpDatabase.GetMotherMass( pdgIndex, mass, width );
   return mass;
 }
