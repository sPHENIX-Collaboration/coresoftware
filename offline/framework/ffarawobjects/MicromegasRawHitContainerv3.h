#ifndef FUN4ALLRAW_MICROMEGASRAWHITCONTAINERv3_H
#define FUN4ALLRAW_MICROMEGASRAWHITCONTAINERv3_H

#include "MicromegasRawHitContainer.h"

class MicromegasRawHit;
class TClonesArray;

// NOLINTNEXTLINE(hicpp-special-member-functions)
class MicromegasRawHitContainerv3 : public MicromegasRawHitContainer
{
 public:

  /// constructor
  explicit MicromegasRawHitContainerv3();

  /// destructor
  ~MicromegasRawHitContainerv3() override;

  /// Clear Event
  void Reset() override;

  /** identify Function from PHObject
      @param os Output Stream
   */
  void identify(std::ostream &os = std::cout) const override;

  /// isValid returns non zero if object contains vailid data
  int isValid() const override;

  MicromegasRawHit *AddHit() override;
  MicromegasRawHit *AddHit(MicromegasRawHit*) override;
  unsigned int get_nhits() override;
  MicromegasRawHit *get_hit(unsigned int) override;

 private:
  TClonesArray *MicromegasRawHitsTCArray{nullptr};
  ClassDefOverride(MicromegasRawHitContainerv3, 1)
};

#endif
