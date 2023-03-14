#ifndef G4TRACKING_TRKRTRUTHTRACKCONTAINER_H
#define G4TRACKING_TRKRTRUTHTRACKCONTAINER_H


#include <trackbase/TrkrDefs.h>
#include <phool/PHObject.h>
#include <map>
#include <vector>

class PHG4TruthInfoContainer;
class TrkrTruthTrack;

class TrkrTruthTrackContainer : public PHObject
{
 public:
  //using Vector = std::vector<TrkrTruthTrack*>;
  using Map = std::map<unsigned int, TrkrTruthTrack*>;
  using Iterator = Map::iterator; //Vector::iterator;
  using ConstIterator = Map::const_iterator;
  using Range = std::pair<Iterator, Iterator>;
  using ConstRange = std::pair<ConstIterator, ConstIterator>;

  //! reset method
  void Reset() override {}

  //! add a Track
  virtual void            addTruthTrack (TrkrTruthTrack*) { };
  virtual TrkrTruthTrack* getTruthTrack (unsigned int /*trackid*/) { return nullptr; };
  virtual TrkrTruthTrack* getTruthTrack (unsigned int /*trackid*/, PHG4TruthInfoContainer*) 
  { return nullptr; };
  virtual ConstRange      getTruthTrackRange () const;
  virtual bool            hasTrackid    (unsigned int /*trackid*/) const { return false;   };
  virtual Map&            getMap();
  int nhw() { return 6; };
  int nhw_cc();
  virtual int nhw_virt () { return 60; };
  
  // PHObject virtual overload
  void identify(std::ostream& os = std::cout) const override
  {
    os << "TrkrTruthTrackContainer base class" << std::endl;
  };

 protected:
  TrkrTruthTrackContainer() = default;
  ClassDefOverride(TrkrTruthTrackContainer, 1)
};

#endif  // G4TRACKING_TRUTHTRACKCONTAINER_H
