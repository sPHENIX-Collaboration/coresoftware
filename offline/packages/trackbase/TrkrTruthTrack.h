#ifndef TRACKBASE_TRKRTRUTHTRACK_H
#define TRACKBASE_TRKRTRUTHTRACK_H

#include "TrkrDefs.h"

#include <phool/PHObject.h>

#include <climits>
#include <cmath>
#include <vector>

class TrkrTruthTrack : public PHObject
{
 public:
  //! dtor
  ~TrkrTruthTrack() override = default;
  TrkrTruthTrack(){};

  virtual unsigned int getTrackid() const { return UINT_MAX; };

  virtual float getX0() const { return NAN; };
  virtual float getY0() const { return NAN; };
  virtual float getZ0() const { return NAN; };

  virtual float getPseudoRapidity() const { return NAN; };
  virtual float getPt() const { return NAN; };
  virtual float getPhi() const { return NAN; };

  virtual std::vector<TrkrDefs::cluskey>& getClusters();
  virtual void addCluster(TrkrDefs::cluskey){};
  // std::map<unsigned int /*track id*/, std::vector<TrkrDefs::cluskey>

  virtual void setTrackid(unsigned int){};
  virtual void setX0(float){};
  virtual void setY0(float){};
  virtual void setZ0(float){};

  virtual void setPseudoRapity(float){};
  virtual void setPt(float){};
  virtual void setPhi(float){};

  struct Comp
  {
    bool operator()(const unsigned int lhs, const TrkrTruthTrack* rhs) const
    {
      return lhs < rhs->getTrackid();
    }
    bool operator()(const TrkrTruthTrack* lhs, const unsigned int& rhs) const
    {
      return lhs->getTrackid() < rhs;
    }
    bool operator()(const TrkrTruthTrack* lhs, const TrkrTruthTrack* rhs) const
    {
      return lhs->getTrackid() < rhs->getTrackid();
    }
  };

  // PHObject virtual overloads
  void identify(std::ostream& os = std::cout) const override
  {
    os << "TrkrTruthTrack base class" << std::endl;
  };

 protected:
  ClassDefOverride(TrkrTruthTrack, 1)
};

#endif  // TRACKBASE_TRKRTRUTHTRACK_H
