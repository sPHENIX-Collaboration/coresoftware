// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef BBC_BBCPMTHIT_H
#define BBC_BBCPMTHIT_H

#include "BbcReturnCodes.h"

#include <phool/PHObject.h>
#include <phool/phool.h>

#include <iostream>

class BbcPmtHit : public PHObject
{
 public:
  BbcPmtHit() {}
  BbcPmtHit(const short /*pmt*/, const float /*adc*/, const float /*tdc0*/, const float /*tdc1*/) {}
  virtual ~BbcPmtHit() {}

  virtual short get_pmt() const
  {
    PHOOL_VIRTUAL_WARNING;
    return BbcReturnCodes::BBC_INVALID_SHORT;
  }
  virtual float get_adc() const
  {
    PHOOL_VIRTUAL_WARNING;
    return BbcReturnCodes::BBC_INVALID_FLOAT;
  }
  virtual float get_tdc0() const
  {
    PHOOL_VIRTUAL_WARNING;
    return BbcReturnCodes::BBC_INVALID_FLOAT;
  }
  virtual float get_tdc1() const
  {
    PHOOL_VIRTUAL_WARNING;
    return BbcReturnCodes::BBC_INVALID_FLOAT;
  }

  void identify(std::ostream& os = std::cout) const override;

 private:
  ClassDefOverride(BbcPmtHit, 1)
};

#endif
