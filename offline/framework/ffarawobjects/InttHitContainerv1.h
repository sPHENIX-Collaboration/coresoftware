#ifndef FUN4ALLRAW_INTTHITCONTAINERV1_H
#define FUN4ALLRAW_INTTHITCONTAINERV1_H

#include "InttHitContainer.h"

class InttHit;
class TClonesArray;

class  InttHitContainerv1: public InttHitContainer
{
public:
  InttHitContainerv1();
  ~InttHitContainerv1() override;

  /// Clear Event
  void Reset() override;

  /** identify Function from PHObject
      @param os Output Stream
   */
  void identify(std::ostream &os = std::cout) const override;

  /// isValid returns non zero if object contains vailid data
  int isValid() const override;

  InttHit *AddHit() override;
  unsigned int get_nhits() override;
  InttHit *get_hit(unsigned int) override {return nullptr;}

private:
    TClonesArray *InttHitsTCArray = nullptr;

  ClassDefOverride(InttHitContainerv1,1)    
};

#endif
