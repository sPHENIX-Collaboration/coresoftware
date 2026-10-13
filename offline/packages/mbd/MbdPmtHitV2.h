#ifndef __MBD_MBDPMTHITV2_H__
#define __MBD_MBDPMTHITV2_H__

#include "MbdPmtHit.h"

#include <cmath>
#include <iostream>
#include <limits>

class MbdPmtHitV2 : public MbdPmtHit
{
 public:
  MbdPmtHitV2() = default;
  ~MbdPmtHitV2() override = default;

  //! Just does a clear
  void Reset() override;

  //! Clear is used by TClonesArray to reset the tower to initial state without calling destructor/constructor
  void Clear(Option_t* = "") override;

  //! PMT number
  Short_t get_pmt() const override { return bpmt; }

  //! Effective Nch in PMT
  Float_t get_q() const override { return bq; }

  //! Best time
  Float_t get_time() const override
  { 
    if ( !std::isnan(btq) && std::isnan(btt) )
    {
      return btq;
    }
    return btt;
  }

  //! Time from time channel
  Float_t get_tt() const override { return btt; }

  //! Time from charge channel
  Float_t get_tq() const override { return btq; }

  //! Chi2/NDF from charge channel waveform fit
  Float_t get_chi2ndf() const override { return (fitstat&0xfff)/100.; }

  //! Info about charge channel waveform fit
  UShort_t get_fitinfo() const override { return ((fitstat&0x7000)>>12); }

  //! Whether time channel was bad
  UShort_t get_badtdc() const override
  {
    if ( (fitstat&0x8000) != 0 )
    {
      return true;
    }
    return false;
  }

  //! Raw Info on Fits
  UShort_t get_fitstat() const override { return fitstat; }

  void set_pmt(const Short_t pmt, const Float_t q, const Float_t tt, const Float_t tq) override
  {
    bpmt = pmt;
    bq = q;
    btt = tt;
    btq = tq;
  }

  void set_fitstat(const UShort_t f) override { fitstat = f; }

  //! Prints out exact identity of object
  void identify(std::ostream& out = std::cout) const override;

  //! isValid returns non zero if object contains valid data
  virtual int isValid() const override
  {
    if (std::isnan(get_time())) return 0;
    return 1;
  }

 private:
  Short_t  bpmt;
  UShort_t fitstat;
  Float_t  bq;
  Float_t  btt;
  Float_t  btq;

  ClassDefOverride(MbdPmtHitV2, 1)
};

#endif
