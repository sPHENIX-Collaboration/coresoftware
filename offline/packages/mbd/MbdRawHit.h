// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef MBD_MBDRAWHIT_H
#define MBD_MBDRAWHIT_H

#include "MbdReturnCodes.h"

#include <phool/PHObject.h>
#include <phool/phool.h>

#include <iostream>

class MbdRawHit : public PHObject
{
 public:
  MbdRawHit() {}
  virtual ~MbdRawHit() override = default;

  virtual Short_t get_pmt() const
  {
    PHOOL_VIRTUAL_WARNING;
    return -9999;
  }

  virtual Float_t get_adc() const
  {
    PHOOL_VIRTUAL_WARNING;
    return MbdReturnCodes::MBD_INVALID_FLOAT;
  }

  virtual Float_t get_ttdc() const
  {
    PHOOL_VIRTUAL_WARNING;
    return MbdReturnCodes::MBD_INVALID_FLOAT;
  }

  virtual Float_t get_qtdc() const
  {
    PHOOL_VIRTUAL_WARNING;
    return MbdReturnCodes::MBD_INVALID_FLOAT;
  }

  virtual void set_pmt(const Short_t /*pmt*/, const Float_t /*adc*/, const Float_t /*ttdc*/, const Float_t /*qtdc*/)
  {
    PHOOL_VIRTUAL_WARNING;
  }

  virtual void identify(std::ostream& out = std::cout) const override;

  virtual int isValid() const override { return 0; }

 private:
  ClassDefOverride(MbdRawHit, 1)
};

#endif
