// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef MBD_MBDPMTHIT_H
#define MBD_MBDPMTHIT_H

#include "MbdReturnCodes.h"

#include <phool/PHObject.h>
#include <phool/phool.h>

#include <iostream>

class MbdPmtHit : public PHObject
{
 public:
  MbdPmtHit() {}
  virtual ~MbdPmtHit() override = default;

  virtual Short_t get_pmt() const
  {
    PHOOL_VIRTUAL_WARNING;
    return -9999;
  }

  virtual Float_t get_q() const
  {
    PHOOL_VIRTUAL_WARNING;
    return MbdReturnCodes::MBD_INVALID_FLOAT;
  }

  virtual Float_t get_time() const
  {
    PHOOL_VIRTUAL_WARNING;
    return MbdReturnCodes::MBD_INVALID_FLOAT;
  }

  virtual Float_t get_tt() const
  {
    PHOOL_VIRTUAL_WARNING;
    return MbdReturnCodes::MBD_INVALID_FLOAT;
  }

  virtual Float_t get_tq() const
  {
    PHOOL_VIRTUAL_WARNING;
    return MbdReturnCodes::MBD_INVALID_FLOAT;
  }

  virtual Float_t get_npe() const
  {
    PHOOL_VIRTUAL_WARNING;
    return MbdReturnCodes::MBD_INVALID_FLOAT;
  }

  virtual Float_t get_chi2ndf() const
  {
    static int ctr = 0;
    if ( ctr<3 )
    {
      PHOOL_VIRTUAL_WARNING;
      ctr++;
    }
    return MbdReturnCodes::MBD_INVALID_FLOAT;
  }

  virtual UShort_t get_fitinfo() const
  {
    static int ctr = 0;
    if ( ctr<3 )
    {
      PHOOL_VIRTUAL_WARNING;
      ctr++;
    }
    return 0;
  }

  virtual UShort_t get_badtdc() const
  {
    static int ctr = 0;
    if ( ctr<3 )
    {
      PHOOL_VIRTUAL_WARNING;
      ctr++;
    }
    return 0;
  }

  virtual UShort_t get_fitstat() const
  {
    static int ctr = 0;
    if ( ctr<3 )
    {
      PHOOL_VIRTUAL_WARNING;
      ctr++;
    }
    return 0;
  }

  virtual void set_pmt(const Short_t /*pmt*/, const Float_t /*q*/, const Float_t /*tt*/, const Float_t /*tq*/)
  {
    PHOOL_VIRTUAL_WARNING;
  }

  virtual void set_simpmt(const Float_t /*npe*/)
  {
    PHOOL_VIRTUAL_WARNING;
  }

  virtual void set_npe(const Float_t /*npe*/)
  {
    PHOOL_VIRTUAL_WARNING;
  }

  virtual void set_fitstat(const UShort_t /*fitstat*/)
  {
    static int ctr = 0;
    if ( ctr<3 )
    {
      PHOOL_VIRTUAL_WARNING;
      ctr++;
    }
  }


  virtual void identify(std::ostream& out = std::cout) const override;

  virtual int isValid() const override { return 0; }

 private:
  ClassDefOverride(MbdPmtHit, 1)
};

#endif
