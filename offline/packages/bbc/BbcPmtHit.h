// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef BBC_BBCPMTHIT_H
#define BBC_BBCPMTHIT_H

#include <iostream>

#include <phool/PHObject.h>
#include <phool/phool.h>


class BbcPmtHit : public PHObject
{
public:

  BbcPmtHit() {}
  BbcPmtHit(const Short_t /*pmt*/, const Float_t /*adc*/, const Float_t /*tdc0*/, const Float_t /*tdc1*/) {}
  virtual ~BbcPmtHit() { }

  virtual Short_t get_pmt()  const { PHOOL_VIRTUAL_WARNING; return -9999; }
  virtual Float_t get_adc()  const { PHOOL_VIRTUAL_WARNING; return -9999; }
  virtual Float_t get_tdc0() const { PHOOL_VIRTUAL_WARNING; return -9999; }
  virtual Float_t get_tdc1() const { PHOOL_VIRTUAL_WARNING; return -9999; }

  void identify(std::ostream& os = std::cout) const override;

private:
  ClassDefOverride(BbcPmtHit,1)
};

#endif
