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

  virtual unsigned int getTrackid () const { return 0.; };

  virtual float getX0() const { return 0.; };
  virtual float getY0() const { return 0.; };
  virtual float getZ0() const { return 0.; };

  virtual float getPseudoRapidity() const { return 0.; };
  virtual float getPt()             const { return 0.; };
  virtual float getPhi()            const { return 0.; };

  virtual std::vector<TrkrDefs::cluskey>& getClusters();
  virtual void addCluster(TrkrDefs::cluskey) {};
  //std::map<unsigned int /*track id*/, std::vector<TrkrDefs::cluskey>

  virtual void setTrackid(unsigned int) {};
  virtual void setX0(float) {};
  virtual void setY0(float) {};
  virtual void setZ0(float) {};

  virtual void setPseudoRapity(float) {};
  virtual void setPt(float)  {};
  virtual void setPhi(float) {};

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

