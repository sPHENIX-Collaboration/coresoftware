#ifndef __TRACKERHIT_H__
#define __TRACKERHIT_H__

#include <PHG4Cell.h>
#include "TrackerDefs.h"

class TrackerHit: public PHG4Cell
{

public:

  virtual ~TrackerHit() {}

  virtual void identify(std::ostream& os = std::cout) const;
  virtual void Copy(TrackerHit const &hit);
  virtual void Reset();

  void set_hitid(TrackerDefs::keytype id) { hitid = id; }
  TrackerDefs::keytype get_hitid() { return hitid; }

protected:

private:

  TrackerDefs::keytype hitid; 

  ClassDef(TrackerHit, 2);

};

#endif