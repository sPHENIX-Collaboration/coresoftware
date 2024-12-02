#ifndef FUN4ALLRAW_INTTRAWHITCONTAINER_H
#define FUN4ALLRAW_INTTRAWHITCONTAINER_H

#include <phool/PHObject.h>

class InttRawHit;

class InttRawHitContainer : public PHObject
{
 public:
  InttRawHitContainer() = default;
  virtual ~InttRawHitContainer() = default;

  virtual InttRawHit *AddHit() { return nullptr; }
  virtual InttRawHit *AddHit(InttRawHit *) { return nullptr; }
  virtual unsigned int get_nhits() { return 0; }
  virtual InttRawHit *get_hit(unsigned int) { return nullptr; }

 private:
  ClassDefOverride(InttRawHitContainer, 1)
};

#endif
