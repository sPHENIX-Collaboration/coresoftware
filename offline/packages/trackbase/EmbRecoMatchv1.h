#ifndef TRACKBASE_EMBRECOMATCHV1_H
#define TRACKBASE_EMBRECOMATCHV1_H

#include "TrkrDefs.h"
#include "EmbRecoMatch.h"

#include <phool/PHObject.h>
#include <vector>

/* class VtxPoint; */
class EmbRecoMatchv1 : public EmbRecoMatch
{
  public:
    ~EmbRecoMatchv1() override = default;
    unsigned short idTruthTrack()    const override { return m_idTruthTrack; };
    unsigned short idRecoTrack()     const override { return m_idRecoTrack;  };
    unsigned short idTrackSeed()     const override { return m_idTrackSeed;  };
    unsigned short idSvtxTrackSeed() const override { return m_idSvtxTrackSeed;  };

    unsigned short nClustersTruth()   const override { return m_nClustersTruth;   };
    unsigned short nClustersReco()    const override { return m_nClustersReco;    };
    unsigned short nClustersMatched() const override { return m_nClustersMatched; };

    /* float meanClusterZDiff()   const override { return m_meanClusterZDiff; }; */
    /* float meanClusterPhiDiff() const override { return m_meanClusterPhiDiff; }; */

    EmbRecoMatchv1() {};
    EmbRecoMatchv1(
          unsigned short id_truth 
        , unsigned short id_reco
        , unsigned short nclustruth   =0
        , unsigned short nclusreco    =0
        , unsigned short nclusmatched =0
        , unsigned short id_trackseed =0
        , unsigned short id_svtxtrackseed=0
        /* float meanclusZdiff, */
        /* float meanclusphidiff */
    );
    /* : */
    /*   m_idTruthTrack    { id_truth }, */
    /*   m_idRecoTrack     { id_reco  }, */
    /*   m_idTrackSeed     { id_trackseed }, */
    /*   m_idSvtxTrackSeed { id_svtxtrackseed }, */
    /*   m_nClustersTruth  { nclustruth }, */
    /*   m_nClustersReco   { nclusreco }, */
    /*   m_meanClusterZDiff { meanclusZdiff }, */
    /*   m_meanClusterPhiDiff { meanclusphidiff } */
  /* {}; */

  private:
    /* unsigned short m_idTruthTrack    {0}; */
    unsigned short m_idRecoTrack     {0};

    unsigned short m_nClustersTruth   {0};
    unsigned short m_nClustersReco    {0};

    unsigned short m_idTrackSeed     {0};
    unsigned short m_idSvtxTrackSeed {0};
    unsigned short m_nClustersMatched {0};

    /* float m_meanClusterZDiff   {0.}; */
    /* float m_meanClusterPhiDiff {0.}; */

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
