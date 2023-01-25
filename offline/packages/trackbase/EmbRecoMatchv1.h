#ifndef TRACKBASE_EMBRECOMATCHV1_H
#define TRACKBASE_EMBRECOMATCHV1_H

#include "TrkrDefs.h"
#include "EmbRecoMatch.h"

#include <phool/PHObject.h>
#include <vector>
#include <array>

/* class VtxPoint; */
class EmbRecoMatchv1 : public EmbRecoMatch
{
  public:
    ~EmbRecoMatchv1() override = default;

    unsigned short idTruthTrack()   const override { return m_idTruthTrack;   };
    unsigned short nClustersTruth() const override { return m_nClustersTruth; };
    unsigned short nMatches()       const override { return m_nMatches;       };


  private:
    static constexpr int MATCH_idReco           = 0;
    static constexpr int MATCH_idTpcSeed        = 1;
    static constexpr int MATCH_idSvtxSeed       = 2;
    static constexpr int MATCH_nClustersReco    = 3;
    static constexpr int MATCH_nClustersMatched = 4;

  public:
    unsigned short idRecoTrack      (unsigned short i=0) const override { return m_matches[i][MATCH_idReco];           };
    unsigned short idTpcTrackSeed   (unsigned short i=0) const override { return m_matches[i][MATCH_idTpcSeed];        };
    unsigned short idSvtxTrackSeed  (unsigned short i=0) const override { return m_matches[i][MATCH_idSvtxSeed];       };

    unsigned short nClustersReco    (unsigned short i=0) const override { return m_matches[i][MATCH_nClustersReco];    };
    unsigned short nClustersMatched (unsigned short i=0) const override { return m_matches[i][MATCH_nClustersMatched]; };

    /* float meanClusterZDiff()   const override { return m_meanClusterZDiff; }; */
    /* float meanClusterPhiDiff() const override { return m_meanClusterPhiDiff; }; */

    EmbRecoMatchv1(
          unsigned short id_truth         = USHRT_MAX
        , unsigned short id_reco          = USHRT_MAX
        , unsigned short nclustruth       = USHRT_MAX
        , unsigned short nclusreco        = USHRT_MAX
        , unsigned short nclusmatched     = USHRT_MAX
        , unsigned short id_tpctrackseed  = USHRT_MAX
        , unsigned short id_svtxtrackseed = USHRT_MAX
    ) 
      : m_idTruthTrack   { id_truth }
      , m_nClustersTruth { nclustruth }
      , m_nMatches       { 1 }
      , m_matches {{ id_reco, id_tpctrackseed, id_svtxtrackseed, nclusreco, nclusmatched }}
    {};

    unsigned short add_match (
          unsigned short id_reco          = USHRT_MAX
        , unsigned short nclusreco        = USHRT_MAX
        , unsigned short nclusmatched     = USHRT_MAX
        , unsigned short id_tpctrackseed  = USHRT_MAX
        , unsigned short id_svtxtrackseed = USHRT_MAX
    );

  private:
    
    unsigned short m_idTruthTrack    {USHRT_MAX};
    unsigned short m_nClustersTruth  {USHRT_MAX};
    unsigned short m_nMatches        {USHRT_MAX};

    std::vector<std::array<unsigned short, 5>> m_matches;
    // each match has metrics defined by five unsigned ints: idReco, idTpcSee, idSvtxSeed, nClustersReco, nClusterMatched

  public:
    // PHObject virtual overloads
    void identify(std::ostream& os = std::cout) const override
    {
      os << "EmbRecoMatchv1 base class" << std::endl;
    };

 protected:
  ClassDefOverride(EmbRecoMatchv1, 1)
};

#endif // G4TPC_TruthTrack_h
