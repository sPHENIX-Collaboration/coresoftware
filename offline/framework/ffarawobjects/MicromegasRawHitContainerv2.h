#ifndef FUN4ALLRAW_MICROMEGASRAWHITCONTAINERV2_H
#define FUN4ALLRAW_MICROMEGASRAWHITCONTAINERV2_H
#include "MicromegasRawHitContainer.h"

class MicromegasRawHit;
class TClonesArray;

class MicromegasRawHitContainerv2 : public MicromegasRawHitContainer
{
  public:

  /// constructor
  explicit MicromegasRawHitContainerv2();

  /// destructor
  ~MicromegasRawHitContainerv2() override;

  /// Clear Event
  void Reset() override;

  /** identify Function from PHObject
  @param os Output Stream
  */
  void identify(std::ostream &os = std::cout) const override;

  /// isValid returns non zero if object contains vailid data
  int isValid() const override;

  MicromegasRawHit *AddHit() override;
  MicromegasRawHit *AddHit(MicromegasRawHit *) override;
  unsigned int get_nhits() override;
  MicromegasRawHit *get_hit(unsigned int) override;

 private:
  TClonesArray *MicromegasRawHitsTCArray = nullptr;

  ClassDefOverride(MicromegasRawHitContainerv2, 1)
};

#endif
