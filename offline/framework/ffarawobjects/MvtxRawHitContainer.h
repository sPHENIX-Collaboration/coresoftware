#ifndef FUN4ALLRAW_MVTXRAWHITCONTAINER_H
#define FUN4ALLRAW_MVTXRAWHITCONTAINER_H

#include <phool/PHObject.h>

class MvtxRawHit;

class MvtxRawHitContainer : public PHObject
{
 public:
  MvtxRawHitContainer() = default;
  virtual ~MvtxRawHitContainer() = default;

  virtual MvtxRawHit *AddHit() { return nullptr; }
  virtual MvtxRawHit *AddHit(MvtxRawHit *) { return nullptr; }
  virtual unsigned int get_nhits() { return 0; }
  virtual MvtxRawHit *get_hit(unsigned int) { return nullptr; }

 private:
  ClassDefOverride(MvtxRawHitContainer, 1)
};

#endif
