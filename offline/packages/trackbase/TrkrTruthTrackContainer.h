#ifndef TRACKBASE_TRKRTRUTHTRACKCONTAINER_H
#define TRACKBASE_TRKRTRUTHTRACKCONTAINER_H

#include "TrkrDefs.h"

#include <phool/PHObject.h>
#include <vector>
#include <map>

class TrkrTruthTrack;
class TrkrTruthTrackContainer : public PHObject
{
  public:
  using Vector        = std::vector<TrkrTruthTrack*>;
  using Iterator      = Vector::iterator;
  using ConstIterator = Vector::const_iterator;
  using Range         = std::pair<Iterator, Iterator>;
  using ConstRange    = std::pair<ConstIterator, ConstIterator>;

  //! reset method
  void Reset() override {}

  //! add a Track
  virtual void       addTruthTrack(TrkrTruthTrack*) = 0;
  virtual ConstRange getTruthTrackRange() const =0;
  virtual bool       hasTrackid(unsigned int trackid) const =0;
  virtual TrkrTruthTrack* getTruthTrack(unsigned int trackid) const=0;
  virtual Vector&    getTruthTracks() =0;
  
  // PHObject virtual overload
  void identify(std::ostream& os = std::cout) const override
  {
    os << "TrkrTruthTrackContainer base class" << std::endl;
  };

 protected:
  TrkrTruthTrackContainer() = default;
  ClassDefOverride(TrkrTruthTrackContainer, 1)
};

#endif // TRACKBASE_TRUTHTRACKCONTAINER_H

