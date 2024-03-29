#ifndef FUN4ALLRAW_GL1RAWTHIT_H
#define FUN4ALLRAW_GL1RAWTHIT_H

#include <phool/PHObject.h>

#include <limits>

class Gl1RawHit : public PHObject
{
 public:
  Gl1RawHit() = default;
  virtual ~Gl1RawHit() = default;

  virtual uint64_t get_bco() const { return std::numeric_limits<uint64_t>::max(); }
  virtual void set_bco(const uint64_t) { return; }
  virtual int getEvtSequence() const { return std::numeric_limits<int>::min(); }
  virtual void setEvtSequence(const int) {return;}
  using PHObject::CopyFrom;
  virtual void CopyFrom(Gl1RawHit *gl1hit);

 private:
  ClassDefOverride(Gl1RawHit, 1)
};

#endif
