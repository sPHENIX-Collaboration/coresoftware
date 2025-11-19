#ifndef __MBD_MBDRAWHITV1_H__
#define __MBD_MBDRAWHITV1_H__

#include "MbdRawHit.h"

#include <cmath>
#include <iostream>
#include <limits>

class MbdRawHitV1 : public MbdRawHit
{
 public:
  MbdRawHitV1() = default;
  ~MbdRawHitV1() override = default;

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

  void set_pmt(const Short_t pmt, const Float_t a, const Float_t tt, const Float_t tq) override
  {
    bpmt = pmt;
    badc = a;
    bttdc = tt;
    bqtdc = tq;
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
  Short_t bpmt;
  Float_t badc;
  Float_t bttdc;
  Float_t bqtdc;

  ClassDefOverride(MbdRawHitV1, 1)
};

#endif
