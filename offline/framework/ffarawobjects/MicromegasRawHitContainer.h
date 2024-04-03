#ifndef FUN4ALLRAW_MICROMEGASRAWHITCONTAINER_H
#define FUN4ALLRAW_MICROMEGASRAWHITCONTAINER_H

#include <phool/PHObject.h>

class MicromegasRawHit;

class MicromegasRawHitContainer : public PHObject
{
 public:
  MicromegasRawHitContainer() = default;
  virtual ~MicromegasRawHitContainer() = default;

  virtual MicromegasRawHit *AddHit() { return nullptr; }
  virtual MicromegasRawHit *AddHit(MicromegasRawHit *) { return nullptr; }
  virtual unsigned int get_nhits() { return 0; }
  virtual MicromegasRawHit *get_hit(unsigned int) { return nullptr; }

 private:
  ClassDefOverride(MicromegasRawHitContainer, 1)
};

#endif
