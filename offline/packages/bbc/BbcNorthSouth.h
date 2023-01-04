#ifndef __BBCNORTHSOUTH_H__
#define __BBCNORTHSOUTH_H__

#include "phool/phool.h"
#include <TObject.h>


class BbcNorthSouth : public TObject
{
public:
  BbcNorthSouth() { }
  BbcNorthSouth(const Short_t /*npmts*/, const Float_t /*ncharge*/, const Float_t /*meantime*/) {}
  virtual ~BbcNorthSouth() { }
  virtual void identify(std::ostream& os = std::cout) const;

  virtual Short_t get_nPMT() const { PHOOL_VIRTUAL_WARNING; return -9999; }
  virtual Float_t get_nCharge() const { PHOOL_VIRTUAL_WARNING; return -9999; }
  virtual Float_t get_MeanTime() const { PHOOL_VIRTUAL_WARNING; return -9999; }

protected:

  virtual void Clear(Option_t * /*option*/ = "") override { PHOOL_VIRTUAL_WARNING; }

private:

  ClassDefOverride(BbcNorthSouth,1)

};

#endif
