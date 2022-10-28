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

    EmbRecoMatchv1(
          unsigned short id_truth     = USHRT_MAX
        , unsigned short id_reco      = USHRT_MAX
        , unsigned short nclustruth   = USHRT_MAX
        , unsigned short nclusreco    = USHRT_MAX
        , unsigned short nclusmatched = USHRT_MAX
        , unsigned short id_trackseed = USHRT_MAX
        , unsigned short id_svtxtrackseed= USHRT_MAX
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
    unsigned short m_idTruthTrack    {USHRT_MAX};
    unsigned short m_idRecoTrack     {USHRT_MAX};

    unsigned short m_nClustersTruth   {USHRT_MAX};
    unsigned short m_nClustersReco    {USHRT_MAX};

    unsigned short m_idTrackSeed      {USHRT_MAX};
    unsigned short m_idSvtxTrackSeed  {USHRT_MAX};
    unsigned short m_nClustersMatched {USHRT_MAX};

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
