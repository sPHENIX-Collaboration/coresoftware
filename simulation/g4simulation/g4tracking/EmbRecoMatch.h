#ifndef TRACKBASE_EMBRECOMATCH_H
#define TRACKBASE_EMBRECOMATCH_H

#include <phool/PHObject.h>
#include <trackbase/TrkrDefs.h>
#include <climits>
#include <cmath>
#include <vector>

class EmbRecoMatch : public PHObject
{
 public:
  EmbRecoMatch(){};
  ~EmbRecoMatch() override = default;

  virtual unsigned short idTruthTrack() const { return std::numeric_limits<unsigned short>::max(); };
  virtual unsigned short nClustersTruth() const { return std::numeric_limits<unsigned short>::max(); };
  virtual unsigned short nMatches() const { return std::numeric_limits<unsigned short>::max(); };

  virtual unsigned short idRecoTrack(unsigned short = 0) const { return std::numeric_limits<unsigned short>::max(); };
  virtual unsigned short idTpcTrackSeed(unsigned short = 0) const { return std::numeric_limits<unsigned short>::max(); };
  virtual unsigned short idSvtxTrackSeed(unsigned short = 0) const { return std::numeric_limits<unsigned short>::max(); };

  virtual unsigned short nClustersReco(unsigned short = 0) const { return std::numeric_limits<unsigned short>::max(); };
  virtual unsigned short nClustersMatched(unsigned short = 0) const { return std::numeric_limits<unsigned short>::max(); };

  virtual float nRatioMatched() const { return std::numeric_limits<float>::quiet_NaN(); };

  struct Comp
  {
    bool operator()(const unsigned int lhs, const EmbRecoMatch* rhs) const
    {
      return lhs < rhs->idTruthTrack();
    }
    bool operator()(const EmbRecoMatch* lhs, const unsigned int rhs) const
    {
      return lhs->idTruthTrack() < rhs;
    }
    bool operator()(const EmbRecoMatch* lhs, const EmbRecoMatch* rhs) const
    {
      if (lhs->idTruthTrack() != rhs->idTruthTrack())
        return lhs->idTruthTrack() < rhs->idTruthTrack();
      else
        return lhs->idRecoTrack() < rhs->idRecoTrack();
    }
  };

 protected:
  ClassDefOverride(EmbRecoMatch, 1)
};

#endif  // G4TPC_TruthTrack_h
