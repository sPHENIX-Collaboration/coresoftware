#ifndef __BBC_BBCPMTCONTAINERV1_H__
#define __BBC_BBCPMTCONTAINERV1_H__

#include "BbcPmtContainer.h"

#include <iostream>

#include <TClonesArray.h>

///
class BbcPmtContainerV1 : public BbcPmtContainer
{
public:
  /// ctor
  BbcPmtContainerV1();

  /// dtor
  virtual ~BbcPmtContainerV1();

  /// Clear Event
  void Reset() override;

  /** identify Function from PHObject
      @param os Output Stream
   */
  void identify(std::ostream &os = std::cout) const override;

  /// isValid returns non zero if object contains vailid data
  int isValid() const override;

  /** set number of pmts for Bbc
      @param ival Number of Bbc Pmt's
   */
  void set_npmt(const Short_t ival) override
  {
    npmt = ival;
    return;
  }

  /// get Number of Bbc Pmt's
  Short_t get_npmt() const override { return npmt; }

  /** get BbcHitPmt of Pmt iPmt in TClonesArray
      @param iPmt no of Pmt in TClonesArray
   */
  BbcPmtHit *get_pmt(const int iPmt) const override { return (BbcPmtHit*)BbcPmtHits->ConstructedAt(iPmt); }

private:
  TClonesArray *GetBbcPmtHits() const { return BbcPmtHits; }

  Short_t npmt = 0;
  TClonesArray *BbcPmtHits = nullptr;

  ClassDefOverride(BbcPmtContainerV1, 1)
};

#endif
