#ifndef TRACKBASE_TRKRTRUTHTRACK_H
#define TRACKBASE_TRKRTRUTHTRACK_H

#include <phool/PHObject.h>
#include <trackbase/TrkrDefs.h>

#include <limits>
#include <vector>

class TrkrTruthTrack : public PHObject
{
 public:
  //! dtor
  ~TrkrTruthTrack() override = default;
  TrkrTruthTrack(){};

  virtual unsigned int getTrackid() const { return std::numeric_limits<unsigned int>::max(); };

  virtual float getX0() const { return std::numeric_limits<float>::quiet_NaN(); };
  virtual float getY0() const { return std::numeric_limits<float>::quiet_NaN(); };
  virtual float getZ0() const { return std::numeric_limits<float>::quiet_NaN(); };

  virtual float getPseudoRapidity() const { return std::numeric_limits<float>::quiet_NaN(); };
  virtual float getPt() const { return std::numeric_limits<float>::quiet_NaN(); };
  virtual float getPhi() const { return std::numeric_limits<float>::quiet_NaN(); };

  virtual std::vector<TrkrDefs::cluskey>& getClusters();
  virtual void addCluster(TrkrDefs::cluskey){};

  virtual bool has_hitsetkey(TrkrDefs::hitsetkey) const { return false; };
  virtual bool has_hitsetkey(TrkrDefs::cluskey) const { return false; };
  virtual std::pair<bool, TrkrDefs::cluskey> get_cluskey(TrkrDefs::hitsetkey) const { return {false, 0}; };  // bool is if there is the key, and if so, then hitsetket is the correspeonding key

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
