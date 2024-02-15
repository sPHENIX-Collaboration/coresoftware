#ifndef FUN4ALLRAW_INTTRAWHITCONTAINERV1_H
#define FUN4ALLRAW_INTTRAWHITCONTAINERV1_H

#include "InttRawHitContainer.h"

class InttRawHit;
class TClonesArray;

class  InttRawHitContainerv1: public InttRawHitContainer
{
public:
  InttRawHitContainerv1();
  ~InttRawHitContainerv1() override;

  /// Clear Event
  void Reset() override;

  /** identify Function from PHObject
      @param os Output Stream
   */
  void identify(std::ostream &os = std::cout) const override;

  /// isValid returns non zero if object contains vailid data
  int isValid() const override;

  InttRawHit *AddHit() override;
  InttRawHit *AddHit(InttRawHit *intthit) override;
  unsigned int get_nhits() override;
  InttRawHit *get_hit(unsigned int index) override;

private:
    TClonesArray *InttRawHitsTCArray = nullptr;

  ClassDefOverride(InttRawHitContainerv1,1)
};

#endif
