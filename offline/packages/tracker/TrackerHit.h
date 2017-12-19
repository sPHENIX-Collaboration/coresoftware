#ifndef __TRACKERHIT_H__
#define __TRACKERHIT_H__

#include <g4main/PHG4Hit.h>
#include <TObject.h>
#include "TrackerDefs.h"

#include <limits.h>
#include <map>

class TrackerHit : public TObject
{
 public:

  typedef std::map<PHG4HitDefs::keytype, float> EdepMap;
  typedef EdepMap::iterator EdepIterator;
  typedef EdepMap::const_iterator EdepConstIterator;
  typedef std::pair<EdepIterator, EdepIterator> EdepRange;
  typedef std::pair<EdepConstIterator, EdepConstIterator> EdepConstRange;

  virtual ~TrackerHit() {}
  virtual void identify(std::ostream& os = std::cout) const;
  virtual void Copy(TrackerHit const& hit);
  virtual void Reset();
  virtual void print() {};

  virtual void add_edep(const PHG4HitDefs::keytype g4hitid, const float edep) {return;}
  virtual EdepConstRange get_g4hits() {
    std::map <PHG4HitDefs::keytype, float> dummy;
    return std::make_pair(dummy.begin(), dummy.end());
  }

  virtual void set_hitid(const TrackerDefs::keytype id) { return; }
  virtual TrackerDefs::keytype get_hitid() const { return 0; }


 protected:
  TrackerHit() {}
 private:

  ClassDef(TrackerHit, 1);
};

#endif