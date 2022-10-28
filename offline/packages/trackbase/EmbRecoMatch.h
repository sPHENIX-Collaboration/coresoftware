#ifndef TRACKBASE_EMBRECOMATCH_H
#define TRACKBASE_EMBRECOMATCH_H

#include "TrkrDefs.h"

#include <phool/PHObject.h>
#include <vector>
#include <climits>

/* class VtxPoint; */

class EmbRecoMatch : public PHObject
{
  public:
    EmbRecoMatch() {};
    ~EmbRecoMatch() override = default;
    // The getters
    virtual unsigned short idTruthTrack()     const { return USHRT_MAX; };
    virtual unsigned short idRecoTrack()      const { return USHRT_MAX; };
    virtual unsigned short idTrackSeed()      const { return USHRT_MAX; };
    virtual unsigned short idSvtxTrackSeed()  const { return USHRT_MAX; };

    virtual unsigned short nClustersTruth()   const { return USHRT_MAX; };
    virtual unsigned short nClustersReco()    const { return USHRT_MAX; };
    virtual unsigned short nClustersMatched() const { return USHRT_MAX; };

    // The setters

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
