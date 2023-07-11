/**
 * @file g4tpc/EmbRecoMatchv1.cc
 * @author D. Stewart
 * @date September 2022
 * @brief Version 1 of EmbRecoMatch
 */
#include "EmbRecoMatchv1.h"

unsigned short EmbRecoMatchv1::add_match(
          unsigned short id_reco          
        , unsigned short nclusreco        
        , unsigned short nclusmatched     
        , unsigned short id_tpctrackseed  
        , unsigned short id_svtxtrackseed 
    ) 
{
  m_matches.push_back( { id_reco, id_tpctrackseed, id_svtxtrackseed, nclusreco, nclusmatched } );
  ++m_nMatches;
  return m_nMatches;
}
