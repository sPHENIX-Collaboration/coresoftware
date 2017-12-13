#ifndef __TRACKERHITV1_H__
#define __TRACKERHITV1_H__

#include "TrackerDefs.h"
#include "TrackerHit.h"

#ifdef __CINT__
#include <stdint.h>
#else
#include <cstdint>
#endif
#include <iostream>
#include <map>

class TrackerHitv1 : public TrackerHit
{
 public:
  TrackerHitv1();
  TrackerHitv1(TrackerDefs::keytype id);

  virtual ~TrackerHitv1();

  virtual void identify(std::ostream& os = std::cout) const;
  virtual void Reset();

  void set_hitid(const TrackerDefs::keytype id) { hitid = id; }
  TrackerDefs::keytype get_hitid() const { return hitid; }

  void add_edep(const PHG4HitDefs::keytype g4hitid, const float edep);

  EdepConstRange get_g4hits()
  {
    return std::make_pair(hitedeps.begin(), hitedeps.end());
  }

  void print() const;

  // D. McGlinchey - When hashing out schema for TrackerHit, included
  //                 prop_map here. However, the property map implementation
  //                 is split between PHG4Cell and PHG4Cellv1.
  //                 Is it truly necessary to reimplement here?
 protected:
 private:

  TrackerDefs::keytype hitid;
  EdepMap hitedeps;
  //! container for additional property
  //  there is a typedef as prop_map_t in PHG4Cellv1, but not included here
};

#endif