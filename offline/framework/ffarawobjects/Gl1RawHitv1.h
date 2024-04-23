#ifndef FUN4ALLRAW_GL1RAWHITV1_H
#define FUN4ALLRAW_GL1RAWHITV1_H

#include "Gl1RawHit.h"

#include <limits>

class Gl1RawHitv1 : public Gl1RawHit
{
 public:
  Gl1RawHitv1() = default;
  Gl1RawHitv1(Gl1RawHit *gl1hit);
  ~Gl1RawHitv1() override{};

  void Reset() override;
  /** identify Function from PHObject
      @param os Output Stream
   */
  void identify(std::ostream &os = std::cout) const override;
  uint64_t get_bco() const override { return bco; }
  // cppcheck-suppress virtualCallInConstructor
  void set_bco(const uint64_t val) override { bco = val; }

 protected:
  uint64_t bco = std::numeric_limits<uint64_t>::max();

  ClassDefOverride(Gl1RawHitv1, 1)
};

#endif
