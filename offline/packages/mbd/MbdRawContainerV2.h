#ifndef MBD_MBDRAWCONTAINERV2_H__
#define MBD_MBDRAWCONTAINERV2_H__

#include "MbdRawContainer.h"

#include <iostream>

#include <TClonesArray.h>

///
class MbdRawContainerV2 : public MbdRawContainer
{
public:
  /// ctor
  MbdRawContainerV2();

  /// dtor
  virtual ~MbdRawContainerV2();

  /// Clear Event
  void Reset() override;

  /** identify Function from PHObject
      @param os Output Stream
   */
  void identify(std::ostream &out = std::cout) const override;

  /// isValid returns non zero if object contains vailid data
  int isValid() const override;

  /** Add Mbd data containing evt, clk, and femclk
      @param ievt   Event number
      @param iclk    XMIT clock
      @param ifemclk FEM clock
   */
  virtual void set_clocks(const Int_t ievt, const UShort_t iclk, const UShort_t ifemclk) override;

  /** get Event Number
   */
  virtual Int_t get_evt() const override;

  /** get XMIT Clock Counter
   */
  virtual UShort_t get_clock() const override;

  /** get FEM Clock Counter
   */
  virtual UShort_t get_femclock() const override;

  /** set number of pmts for Mbd
      @param ival Number of Mbd Pmt's
   */
  void set_npmt(const Short_t ival) override
  {
    if ( ival != MbdRawHits->GetEntries() )
    {
      std::cout << "ERROR, " << ival << " differs from " << MbdRawHits->GetEntries() << std::endl;
      std::cout << " Setting npmt to " << MbdRawHits->GetEntries() << std::endl; 
    }
    npmt = MbdRawHits->GetEntries();
    return;
  }

  /// get Number of Mbd Pmt's
  Short_t get_npmt() const override { return MbdRawHits->GetEntries(); }

  /** get MbdRawPmt of Pmt iPmt in TClonesArray
      @param iPmt no of Pmt in TClonesArray
   */
  MbdRawHit *get_pmt(const int iPmt) const override { return (MbdRawHit*)MbdRawHits->ConstructedAt(iPmt); }

private:
  TClonesArray *GetMbdRawHits() const { return MbdRawHits; }

  Int_t evt{-1};
  UShort_t clk{0};
  UShort_t femclk{0};
  Short_t npmt = 0;
  TClonesArray *MbdRawHits = nullptr;

  ClassDefOverride(MbdRawContainerV2, 1)
};

#endif
