#ifndef FUN4ALLRAW_MVTXHITRAWCONTAINERV1_H
#define FUN4ALLRAW_MVTXHITRAWCONTAINERV1_H

#include "MvtxRawHitContainer.h"

class MvtxRawHit;
class TClonesArray;

class  MvtxRawHitContainerv1: public MvtxRawHitContainer
{
public:
  MvtxRawHitContainerv1();
  ~MvtxRawHitContainerv1() override;

  /// Clear Event
  void Reset() override;

  /** identify Function from PHObject
      @param os Output Stream
   */
  void identify(std::ostream &os = std::cout) const override;

  /// isValid returns non zero if object contains vailid data
  int isValid() const override;

  MvtxRawHit *AddHit() override;
  MvtxRawHit *AddHit(MvtxRawHit *mvtxhit) override;
  unsigned int get_nhits() override;
  MvtxRawHit *get_hit(unsigned int index) override;

private:
    TClonesArray *MvtxRawHitsTCArray = nullptr;

  ClassDefOverride(MvtxRawHitContainerv1,1)
};

#endif
