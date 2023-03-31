// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef BBC_BBCNORTHSOUTH_H
#define BBC_BBCNORTHSOUTH_H

#include <phool/phool.h>
#include <phool/PHObject.h>


class BbcNorthSouth : public PHObject
{
public:
  BbcNorthSouth() = default;

  ~BbcNorthSouth() override = default;
  virtual void identify(std::ostream& os = std::cout) const override;

  virtual Short_t get_nPMT() const { PHOOL_VIRTUAL_WARNING; return -9999; }
  virtual float get_nCharge() const { PHOOL_VIRTUAL_WARNING; return -9999; }
  virtual float get_MeanTime() const { PHOOL_VIRTUAL_WARNING; return -9999; }

private:

  ClassDefOverride(BbcNorthSouth,1)

};

#endif
