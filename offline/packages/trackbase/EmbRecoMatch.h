#ifndef TRACKBASE_EMBRECOMATCH_H
#define TRACKBASE_EMBRECOMATCH_H

#include "TrkrDefs.h"

#include <phool/PHObject.h>
#include <vector>

/* class VtxPoint; */

class EmbRecoMatch : public PHObject
{
  public:
    ~EmbRecoMatch() override = default;
    // The getters
    virtual unsigned short idTruthTrack() const = 0;
    virtual unsigned short idRecoTrack()  const = 0;
    virtual unsigned short idTrackSeed()  const = 0;
    virtual unsigned short idSvtxTrackSeed()  const = 0;

    virtual unsigned short nClustersTruth()   const = 0;
    virtual unsigned short nClustersReco()    const = 0;
    virtual unsigned short nClustersMatched() const = 0;

    /* virtual float meanClusterZDiff()   const = 0; // weighted difference track centers */
    /* virtual float meanClusterPhiDiff() const = 0; */

    // The setters

  struct Comp { 
    bool operator()(const unsigned int lhs, const EmbRecoMatch* rhs) const
    {return lhs < rhs->m_idTruthTrack;}
    bool operator()(const EmbRecoMatch* lhs, const unsigned int rhs) const
    {return lhs->m_idTruthTrack < rhs;}
    bool operator()(const EmbRecoMatch* lhs, const EmbRecoMatch* rhs) const
    {return lhs->m_idTruthTrack < rhs->m_idTruthTrack;}
  };


  protected:
  EmbRecoMatch(unsigned short _idTruthTrack) : 
    m_idTruthTrack { _idTruthTrack} {}; //
  unsigned short m_idTruthTrack;
  EmbRecoMatch() : m_idTruthTrack{0} {};
  ClassDefOverride(EmbRecoMatch, 1)
};
// this and that
/* this is a comment, and I *really* care about it! */

#endif // G4TPC_TruthTrack_h
