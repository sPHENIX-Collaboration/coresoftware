#ifndef __TrkrHitSet_H__
#define __TrkrHitSet_H__

#include <g4main/PHG4Hit.h>
#include <TObject.h>
#include "TrackerDefs.h"

#ifdef __CINT__
#include <stdint.h>
#else
#include <cstdint>
#endif
#include <iostream>

class TrkrHitSet : public TObject
{
 public:

  TrkrHitSet() {}

  virtual ~TrkrHitSet() {}
  virtual void identify(std::ostream& os = std::cout) const;
  virtual void Copy(TrkrHitSet const& hit);
  virtual void Reset();
  virtual void print() {};

  virtual void set_hitid(const TrackerDefs::hitkeytype id) { return; }
  virtual TrackerDefs::hitkeytype get_hitid() const { return 0; }

  virtual void set_truthid(const uint64_t id) { return; }
  virtual uint64_t get_truthid() const { return NULL; }

 private:

  ClassDef(TrkrHitSet, 1);
};

#endif