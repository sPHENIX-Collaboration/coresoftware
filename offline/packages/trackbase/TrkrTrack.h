#ifndef TRACKBASE_TRKRTRACK_H
#define TRACKBASE_TRKRTRACK_H

#include "TrkrDefs.h"

#include <phool/PHObject.h>

#include <iostream>


class TrkrTrack : public PHObject
{
 public:
  //! dtor
  virtual ~TrkrTrack() {}
  // PHObject virtual overloads
  virtual void identify(std::ostream& os = std::cout) const
  {
    os << "TrkrTrack base class" << std::endl;
  }
  virtual void Reset() {}
  virtual int isValid() const { return 0; }
 protected:
  TrkrTrack() {}
  ClassDef(TrkrTrack, 1);
};

#endif //TRACKBASE_TRKRTRACK_H
