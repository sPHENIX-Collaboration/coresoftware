#ifndef __MBD_MBDRAWHITV2_H__
#define __MBD_MBDRAWHITV2_H__

#include "MbdRawHit.h"

#include <cmath>
#include <iostream>
#include <limits>

class MbdRawHitV2 : public MbdRawHit
{
 public:
  MbdRawHitV2() = default;
  ~MbdRawHitV2() override = default;

  //! Just does a clear
  void Reset() override;

  //! Clear is used by TClonesArray to reset the tower to initial state without calling destructor/constructor
  void Clear(Option_t* = "") override;

  //! PMT number
  Short_t get_pmt() const override { return bpmt; }

  //! ADC
  Float_t get_adc() const override { return badc; }

  //! TDC from time channel
  Float_t get_ttdc() const override { return bttdc; }

  //! TDC from charge channel
  Float_t get_qtdc() const override { return bqtdc; }

  //! Chi2/NDF from charge channel waveform fit
  Float_t get_chi2ndf() const override { return (fitstat&0xfff)/100.; }

  //! Info about charge channel waveform fit
  UShort_t get_fitinfo() const override { return (fitstat>>12); }

  //! Set PMT data values
  void set_pmt(const Short_t pmt, const Float_t a, const Float_t tt, const Float_t tq) override
  {
    bpmt = pmt;
    badc = a;
    bttdc = tt;
    bqtdc = tq;
  }

  //! Store chi2/ndf (encoded in fitstat)
  void set_chi2ndf(const Double_t chi2ndf) override
  {
    UShort_t us_chi2ndf = 0;
    if (std::isfinite(chi2ndf) && chi2ndf > 0.)
    {
      const Double_t clipped = (chi2ndf > 40.95) ? 40.95 : chi2ndf;
      us_chi2ndf = static_cast<UShort_t>(clipped * 100.);
    }
    fitstat &= 0xf000;
    fitstat |= us_chi2ndf;
  }

  //! Store fitinfo (encoded in fitstat)
  void set_fitinfo(const UShort_t fitinfo) override
  {
    fitstat &= 0xfff;
    fitstat |= (fitinfo<<12);
  }

  //! Prints out exact identity of object
  void identify(std::ostream& out = std::cout) const override;

  //! isValid returns non zero if object contains valid data
  virtual int isValid() const override
  {
    if (std::isnan(get_ttdc())) return 0;
    return 1;
  }

 private:
  Short_t  bpmt;
  UShort_t fitstat; //waveform fit status
  Float_t  badc;
  Float_t  bttdc;
  Float_t  bqtdc;

  ClassDefOverride(MbdRawHitV2, 1)
};

#endif
