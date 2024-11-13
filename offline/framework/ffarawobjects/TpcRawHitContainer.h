#ifndef FUN4ALLRAW_TPCRAWHITCONTAINER_H
#define FUN4ALLRAW_TPCRAWHITCONTAINER_H

#include <phool/PHObject.h>

class TpcRawHit;

class TpcRawHitContainer : public PHObject
{
 public:
  TpcRawHitContainer() = default;
  virtual ~TpcRawHitContainer() = default;

  virtual TpcRawHit *AddHit() { return nullptr; }
  virtual TpcRawHit *AddHit(TpcRawHit *) { return nullptr; }
  virtual unsigned int get_nhits() { return 0; }
  virtual TpcRawHit *get_hit(unsigned int) { return nullptr; }
  virtual void setStatus(const unsigned int) { return; }
  virtual unsigned int getStatus() const { return 0; }
  virtual void setBco(const uint64_t) { return; }
  virtual uint64_t getBco() const { return 0; }

 private:
  ClassDefOverride(TpcRawHitContainer, 0)
};

#endif
