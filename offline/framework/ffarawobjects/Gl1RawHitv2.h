#ifndef FUN4ALLRAW_GL1RAWHITV2_H
#define FUN4ALLRAW_GL1RAWHITV2_H

#include "Gl1RawHitv1.h"

#include <limits>

class Gl1RawHitv2 : public Gl1RawHitv1
{
 public:
  Gl1RawHitv2() = default;
  Gl1RawHitv2(Gl1RawHit *gl1hit);
  ~Gl1RawHitv2() override{};

  void Reset() override;
  /** identify Function from PHObject
      @param os Output Stream
   */
  void identify(std::ostream &os = std::cout) const override;
  // cppcheck-suppress virtualCallInConstructor
  int getEvtSequence() const override { return evtseq; }
  // cppcheck-suppress virtualCallInConstructor
  void setEvtSequence(const int i) override { evtseq = i; }

 protected:
  int evtseq{std::numeric_limits<int>::min()};

  ClassDefOverride(Gl1RawHitv2, 1)
};

#endif
