#ifndef FUN4ALLRAW_TPCHITRAWCONTAINERv3_H
#define FUN4ALLRAW_TPCHITRAWCONTAINERv3_H

#include "TpcRawHitContainer.h"

class TpcRawHit;
class TClonesArray;

// NOLINTNEXTLINE(hicpp-special-member-functions)
class TpcRawHitContainerv3 : public TpcRawHitContainer
{
 public:
  TpcRawHitContainerv3();
  ~TpcRawHitContainerv3() override;

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
  void setStatus(const unsigned int i) override { status = i; }
  unsigned int getStatus() const override { return status; }
  void setBco(const uint64_t i) override { bco = i; }
  uint64_t getBco() const override { return bco; }

 private:
  TClonesArray *TpcRawHitsTCArray{nullptr};
  uint64_t bco{0};
  unsigned int status{0};

  ClassDefOverride(TpcRawHitContainerv3, 1)
};

#endif
