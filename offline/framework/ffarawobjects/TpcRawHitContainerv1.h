#ifndef FUN4ALLRAW_TPCHITRAWCONTAINERV1_H
#define FUN4ALLRAW_TPCHITRAWCONTAINERV1_H

#include "TpcRawHitContainer.h"

class TpcRawHit;
class TClonesArray;

class  TpcRawHitContainerv1: public TpcRawHitContainer
{
public:
  TpcRawHitContainerv1();
  ~TpcRawHitContainerv1() override;

  /// Clear Event
  void Reset() override;

  /** identify Function from PHObject
      @param os Output Stream
   */
  void identify(std::ostream &os = std::cout) const override;

  /// isValid returns non zero if object contains vailid data
  int isValid() const override;

  TpcRawHit *AddHit() override;
  TpcRawHit *AddHit(TpcRawHit *tpchit) override;
  unsigned int get_nhits() override;
  TpcRawHit *get_hit(unsigned int index) override;

private:
    TClonesArray *TpcRawHitsTCArray = nullptr;

  ClassDefOverride(TpcRawHitContainerv1,1)
};

#endif
