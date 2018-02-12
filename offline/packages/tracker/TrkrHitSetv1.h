#ifndef __TrkrHitSetV1_H__
#define __TrkrHitSetV1_H__

#include "TrackerDefs.h"
#include "TrkrHitSet.h"

#ifdef __CINT__
#include <stdint.h>
#else
#include <cstdint>
#endif
#include <iostream>
#include <map>

class TrkrHitSetv1 : public TrkrHitSet
{
 public:
  TrkrHitSetv1();

  virtual ~TrkrHitSetv1() {}

  virtual void identify(std::ostream& os = std::cout) const;
  virtual void Reset();
  void print() const;

  void set_hitid(const TrackerDefs::hitkeytype id) { hitid = id; }
  TrackerDefs::hitkeytype get_hitid() const { return hitid; }

  void set_truthid(const uint64_t id) { truthid = id; }
  uint64_t get_truthid() const { return truthid; }


 protected:
 private:

  TrackerDefs::hitkeytype hitid;
  uint64_t truthid;

};

#endif