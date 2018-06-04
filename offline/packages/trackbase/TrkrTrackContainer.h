#ifndef TRACKBASE_TRKRTRACKCONTAINER_H
#define TRACKBASE_TRKRTRACKCONTAINER_H

#include "TrkrTrack.h"
#include "TrkrDefs.h"

#include <phool/PHObject.h>

#include <map>
#include <set>

class TrkrTrackContainer : public PHObject
{
 public:

  TrkrTrackContainer() {};

  virtual ~TrkrTrackContainer() {}
  void Reset() {};

  void identify(std::ostream &os = std::cout) const {};

 protected:

  ClassDef(TrkrTrackContainer, 1)
};

#endif //TRACKBASE_TRKRTRACKCONTAINER_H
