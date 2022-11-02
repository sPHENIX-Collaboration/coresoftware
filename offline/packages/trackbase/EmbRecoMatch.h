#ifndef TRACKBASE_EMBRECOMATCH_H
#define TRACKBASE_EMBRECOMATCH_H

#include "TrkrDefs.h"

#include <phool/PHObject.h>
#include <vector>
#include <climits>
#include <cmath>

/* class VtxPoint; */

class EmbRecoMatch : public PHObject
{
  public:
    EmbRecoMatch() {};
    ~EmbRecoMatch() override = default;

    virtual unsigned short idTruthTrack   () const { return USHRT_MAX; };
    virtual unsigned short nClustersTruth () const { return USHRT_MAX; };
    virtual unsigned short nMatches       () const { return 0;         };

    // the argument i is which match (ranked from 0 to nMatches-1, ordered by goodness of match)
    virtual unsigned short idRecoTrack      (unsigned short =0) const { return USHRT_MAX; };
    virtual unsigned short idTpcTrackSeed   (unsigned short =0) const { return USHRT_MAX; };
    virtual unsigned short idSvtxTrackSeed  (unsigned short =0) const { return USHRT_MAX; };

    virtual unsigned short nClustersReco    (unsigned short =0) const { return USHRT_MAX; };
    virtual unsigned short nClustersMatched (unsigned short =0) const { return USHRT_MAX; };

    virtual float nRatioMatched() const { return NAN; };

  struct Comp { 
    bool operator()(const unsigned int lhs, const EmbRecoMatch* rhs) const
    {return lhs < rhs->idTruthTrack();}
    bool operator()(const EmbRecoMatch* lhs, const unsigned int rhs) const
    {return lhs->idTruthTrack() < rhs;}
    bool operator()(const EmbRecoMatch* lhs, const EmbRecoMatch* rhs) const
    {return lhs->idTruthTrack() < rhs->idTruthTrack();}
  };

  protected:
  ClassDefOverride(EmbRecoMatch, 1)
};
// this and that
/* this is a comment, and I *really* care about it! */

#endif // G4TPC_TruthTrack_h
