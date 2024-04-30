#ifndef FUN4ALLRAW_INTTRAWHITV2_H
#define FUN4ALLRAW_INTTRAWHITV2_H

#include "InttRawHitv1.h"

#include <limits>

class InttRawHitv2 : public InttRawHitv1
{
 public:
  InttRawHitv2() = default;
  InttRawHitv2(InttRawHit *intthit);
  ~InttRawHitv2() override = default;

  /** identify Function from PHObject
      @param os Output Stream
   */
  void identify(std::ostream &os = std::cout) const override;

  uint32_t get_event_counter() const override { return event_counter; }
  void set_event_counter(uint32_t val) override { event_counter = val; }

 protected:
  uint32_t event_counter = std::numeric_limits<uint32_t>::max();

  ClassDefOverride(InttRawHitv2, 1)
};

#endif
