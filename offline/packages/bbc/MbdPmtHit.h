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

  virtual void set_pmt(const Short_t /*pmt*/, const Float_t /*q*/, const Float_t /*tt*/, const Float_t /*tq*/)
  {
    PHOOL_VIRTUAL_WARNING;
  }

  virtual void identify(std::ostream& os = std::cout) const override;

  virtual int isValid() const override { return 0; }

 private:
  ClassDefOverride(MbdPmtHit, 1)
};

#endif
