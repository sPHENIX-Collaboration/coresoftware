#ifndef __MBD_MBDPMTCONTAINERV1_H__
#define __MBD_MBDPMTCONTAINERV1_H__

#include "MbdPmtContainer.h"

#include <iostream>

#include <TClonesArray.h>

///
class MbdPmtContainerV1 : public MbdPmtContainer
{
public:
  /// ctor
  MbdPmtContainerV1();

  /// dtor
  virtual ~MbdPmtContainerV1();

  /// Clear Event
  void Reset() override;

  /** identify Function from PHObject
      @param os Output Stream
   */
  void identify(std::ostream &os = std::cout) const override;

  /// isValid returns non zero if object contains vailid data
  int isValid() const override;

  /** set number of pmts for Mbd
      @param ival Number of Mbd Pmt's
   */
  void set_npmt(const Short_t ival) override
  {
    npmt = ival;
    return;
  }

  /// get Number of Mbd Pmt's
  Short_t get_npmt() const override { return npmt; }

  /** get MbdHitPmt of Pmt iPmt in TClonesArray
      @param iPmt no of Pmt in TClonesArray
   */
  MbdPmtHit *get_pmt(const int iPmt) const override { return (MbdPmtHit*)MbdPmtHits->ConstructedAt(iPmt); }

private:
  TClonesArray *GetMbdPmtHits() const { return MbdPmtHits; }

  Short_t npmt = 0;
  TClonesArray *MbdPmtHits = nullptr;

  ClassDefOverride(MbdPmtContainerV1, 1)
};

#endif
