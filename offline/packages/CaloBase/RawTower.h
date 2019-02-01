#ifndef CALOBASE_RAWTOWER_H
#define CALOBASE_RAWTOWER_H

#include "RawTowerDefs.h"

#include <g4detectors/PHG4CellDefs.h>

#include <phool/PHObject.h>
#include <phool/phool.h>

#include <cmath>
#include <iostream>
#include <map>

class RawTower : public PHObject
{
 public:
  //! key type for cell map which should be consistent with CellKeyType
  typedef unsigned long long CellKeyType;

  typedef std::map<CellKeyType, float> CellMap;
  typedef CellMap::iterator CellIterator;
  typedef CellMap::const_iterator CellConstIterator;
  typedef std::pair<CellIterator, CellIterator> CellRange;
  typedef std::pair<CellConstIterator, CellConstIterator> CellConstRange;

  typedef std::map<int, float> ShowerMap;
  typedef ShowerMap::iterator ShowerIterator;
  typedef ShowerMap::const_iterator ShowerConstIterator;
  typedef std::pair<ShowerIterator, ShowerIterator> ShowerRange;
  typedef std::pair<ShowerConstIterator, ShowerConstIterator> ShowerConstRange;

  virtual ~RawTower() {}

  virtual void Reset() { PHOOL_VIRTUAL_WARNING; }
  virtual int isValid() const
  {
    PHOOL_VIRTUAL_WARN("isValid()");
    return 0;
  }
  virtual void identify(std::ostream& os = std::cout) const { PHOOL_VIRTUAL_WARN("identify()"); }

  virtual void set_id(RawTowerDefs::keytype id) { PHOOL_VIRTUAL_WARN("set_id()"); }
  virtual RawTowerDefs::keytype get_id() const
  {
    PHOOL_VIRTUAL_WARN("get_id()");
    return 0;
  }

  virtual void set_key(RawTowerDefs::keytype id) { set_id(id); }
  virtual RawTowerDefs::keytype get_key() const { return get_id(); }

  virtual int get_bineta() const
  {
    PHOOL_VIRTUAL_WARN("get_ieta()");
    return -1;
  }
  virtual int get_binphi() const
  {
    PHOOL_VIRTUAL_WARN("get_iphi()");
    return -1;
  }
  virtual int get_column() const
  {
    PHOOL_VIRTUAL_WARN("get_column()");
    return -1;
  }
  virtual int get_row() const
  {
    PHOOL_VIRTUAL_WARN("get_row()");
    return -1;
  }

  //! energy assigned to the tower. Depending on stage of process and DST node name, it could be energy deposition, light yield or calibrated energies
  virtual double get_energy() const
  {
    PHOOL_VIRTUAL_WARN("get_energy()");
    return 0.0;
  }
  //! energy assigned to the tower. Depending on stage of process and DST node name, it could be energy deposition, light yield or calibrated energies
  virtual void set_energy(const double)
  {
    PHOOL_VIRTUAL_WARN("set_energy()");
    return;
  }

  //! Time stamp assigned to the tower. Depending on the tower maker, it could be rise time or peak time.
  virtual float get_time() const
  {
    PHOOL_VIRTUAL_WARN("get_time()");
    return NAN;
  }
  //! Time stamp assigned to the tower. Depending on the tower maker, it could be rise time or peak time.
  virtual void set_time(const float)
  {
    PHOOL_VIRTUAL_WARN("set_time()");
    return;
  }

  virtual bool empty_g4cells() const { return true; }
  virtual size_t size_g4cells() const { return 0; }
  virtual CellConstRange get_g4cells() const
  {
    PHOOL_VIRTUAL_WARN("get_g4cells()");
    CellMap dummy;
    return make_pair(dummy.begin(), dummy.end());
  }
  virtual CellIterator find_g4cell(CellKeyType id) { return CellMap().end(); }
  virtual CellConstIterator find_g4cell(CellKeyType id) const { return CellMap().end(); }
  virtual void add_ecell(const CellKeyType g4cellid, const float ecell)
  {
    PHOOL_VIRTUAL_WARN("add_ecell(const CellKeyType g4cellid, const float ecell)");
    return;
  }
  virtual void clear_g4cells() {}

  virtual bool empty_g4showers() const { return true; }
  virtual size_t size_g4showers() const { return 0; }
  virtual ShowerConstRange get_g4showers() const
  {
    PHOOL_VIRTUAL_WARN("get_g4showers()");
    ShowerMap dummy;
    return make_pair(dummy.begin(), dummy.end());
  }
  virtual ShowerIterator find_g4shower(int id) { return ShowerMap().end(); }
  virtual ShowerConstIterator find_g4shower(int id) const { return ShowerMap().end(); }
  virtual void add_eshower(const int g4showerid, const float eshower)
  {
    PHOOL_VIRTUAL_WARN("add_eshower(const unsigned int g4showerid, const float eshower)");
    return;
  }
  virtual void clear_g4showers() {}

 protected:
  RawTower() {}

  ClassDef(RawTower, 1)
};

#endif
