/**
 * @file g4tpc/EmbRecoMatchv1.cc
 * @author D. Stewart
 * @date September 2022
 * @brief Version 1 of EmbRecoMatch
 */
#include "EmbRecoMatchv1.h"

EmbRecoMatchv1::EmbRecoMatchv1(
          unsigned short id_truth
        , unsigned short id_reco
        , unsigned short nclustruth
        , unsigned short nclusreco
        , unsigned short nclusmatched
        , unsigned short id_trackseed
        , unsigned short id_svtxtrackseed
    ) : m_idTruthTrack     { id_truth         }
      , m_idRecoTrack      { id_reco          }
      , m_nClustersTruth   { nclustruth       }
      , m_nClustersReco    { nclusreco        }
      , m_idTrackSeed      { id_trackseed     }
      , m_idSvtxTrackSeed  { id_svtxtrackseed }
      , m_nClustersMatched { nclusmatched     }
  {}
