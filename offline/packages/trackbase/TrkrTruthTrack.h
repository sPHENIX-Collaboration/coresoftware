#ifndef TRACKBASE_TRKRTRUTHTRACK_H
#define TRACKBASE_TRKRTRUTHTRACK_H

#include "TrkrDefs.h"

#include <phool/PHObject.h>
#include <vector>

/* class VtxPoint; */

class TrkrTruthTrack : public PHObject
{
 public:

  //! dtor
  ~TrkrTruthTrack() override = default;

  virtual unsigned int getTrackid () const = 0;

  virtual float getX0() const = 0;
  virtual float getY0() const = 0;
  virtual float getZ0() const = 0;

  virtual float getPseudoRapidity() const = 0;
  virtual float getPt() const = 0;
  virtual float getPhi() const = 0;

  virtual std::vector<TrkrDefs::cluskey>& getClusters() = 0;
  virtual void addCluster(TrkrDefs::cluskey) = 0;
  //std::map<unsigned int /*track id*/, std::vector<TrkrDefs::cluskey>

  virtual void setTrackid(unsigned int) = 0;
  virtual void setX0(float) = 0;
  virtual void setY0(float) = 0;
  virtual void setZ0(float) = 0;

  virtual void setPseudoRapity(float) = 0;
  virtual void setPt(float) = 0;
  virtual void setPhi(float) = 0;

  struct Comp { 
    bool operator()(const unsigned int lhs, const TrkrTruthTrack* rhs) const
    {return lhs < rhs->trackid;}
    bool operator()(const TrkrTruthTrack* lhs, const unsigned int& rhs) const
    {return lhs->trackid < rhs;}
    bool operator()(const TrkrTruthTrack* lhs, const TrkrTruthTrack* rhs) const
    {return lhs->trackid < rhs->trackid;}
  };

  // PHObject virtual overloads
  void identify(std::ostream& os = std::cout) const override
  {
    os << "TrkrTruthTrack base class" << std::endl;
  };

 protected:
  unsigned int trackid;
  TrkrTruthTrack() : trackid{0} {};
  ClassDefOverride(TrkrTruthTrack, 1)
};

#endif // G4TPC_TruthTrack_h

