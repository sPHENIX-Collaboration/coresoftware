#ifndef RAWTOWER_H_
#define RAWTOWER_H_

#include <phool/phool.h>
#include <phool/PHObject.h>
#include <iostream>
#include <map>

class RawTower : public PHObject {

 public:
  RawTower() {}
  virtual ~RawTower() {}

  virtual void Reset() { PHOOL_VIRTUAL_WARNING; }
  virtual int isValid() const { PHOOL_VIRTUAL_WARN("isValid()"); return 0; }
  virtual void identify(std::ostream& os=std::cout) const { PHOOL_VIRTUAL_WARN("identify()"); }

  virtual int get_bineta() const { PHOOL_VIRTUAL_WARN("get_ieta()"); return -1; }
  virtual int get_binphi() const { PHOOL_VIRTUAL_WARN("get_iphi()");; return -1; }
  virtual float get_energy() const { PHOOL_VIRTUAL_WARN("get_energy()"); return 0.0; }

  virtual bool is_adjacent(RawTower& tower) { PHOOL_VIRTUAL_WARNING; return false; }

  virtual std::pair< std::map<unsigned int,float>::const_iterator, std::map<unsigned int,float>::const_iterator > get_g4cells();
  virtual void add_ecell(const unsigned int g4cellid, const float ecell) {PHOOL_VIRTUAL_WARN("add_ecell(const unsigned int g4cellid, const float ecell)"); return;}

  ClassDef(RawTower,1)

};

#endif /* RAWTOWER_H_ */
