#ifndef TRACKBASE_TRKRTRUTHTRACKCONTAINER_H
#define TRACKBASE_TRKRTRUTHTRACKCONTAINER_H

#include "TrkrDefs.h"

#include <phool/PHObject.h>
#include <map>
#include <vector>

class TrkrTruthTrack;
class TrkrTruthTrackContainer : public PHObject
{
 public:
  using Vector = std::vector<TrkrTruthTrack*>;
  using Iterator = Vector::iterator;
  using ConstIterator = Vector::const_iterator;
  using Range = std::pair<Iterator, Iterator>;
  using ConstRange = std::pair<ConstIterator, ConstIterator>;

  //! reset method
  void Reset() override {}

  //! add a Track
  virtual void            addTruthTrack(TrkrTruthTrack*) {};
  virtual ConstRange      getTruthTrackRange() const;
  virtual bool            hasTrackid    (unsigned int /*trackid*/) const { return false;   };
  virtual TrkrTruthTrack* getTruthTrack (unsigned int /*trackid*/) const { return nullptr; };
  virtual Vector&         getTruthTracks();
  
  // PHObject virtual overload
  void identify(std::ostream& os = std::cout) const override
  {
    os << "TrkrTruthTrackContainer base class" << std::endl;
  };

 protected:
  TrkrTruthTrackContainer() = default;
  ClassDefOverride(TrkrTruthTrackContainer, 1)
};

#endif  // TRACKBASE_TRUTHTRACKCONTAINER_H
