#ifndef CALOBASE_RAWTOWER_H
#define CALOBASE_RAWTOWER_H

#include "RawTowerDefs.h"

#include <phool/PHObject.h>
#include <phool/phool.h>

#include <climits>
#include <cmath>
#include <cstddef>  // for size_t
#include <iostream>
#include <map>
#include <string>  // for string
#include <utility>

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

  ~RawTower() override {}

  void Reset() override { PHOOL_VIRTUAL_WARNING; }
  int isValid() const override
  {
    PHOOL_VIRTUAL_WARN("isValid()");
    return 0;
  }
  void identify(std::ostream& /*os*/ = std::cout) const override { PHOOL_VIRTUAL_WARN("identify()"); }

  virtual void set_id(RawTowerDefs::keytype) { PHOOL_VIRTUAL_WARN("set_id()"); }
  virtual RawTowerDefs::keytype get_id() const
  {
    PHOOL_VIRTUAL_WARN("get_id()");
    return 0;
  }

  virtual void set_key(RawTowerDefs::keytype id) { set_id(id); }
  virtual RawTowerDefs::keytype get_key() const { return get_id(); }

  virtual int get_bineta() const
  {
    PHOOL_VIRTUAL_WARN("get_bineta()");
    return -1;
  }

  virtual int get_bintheta() const
  {
    return get_bineta();
  }

  virtual int get_binphi() const
  {
    PHOOL_VIRTUAL_WARN("get_binphi()");
    return -1;
  }

  //! return layer ID assigned to tower
  virtual int get_binl() const
  {
    PHOOL_VIRTUAL_WARN("get_binl()");
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

  //! energy assigned to the tower. Depending on stage of process and DST node name, it could be energy deposition, light yield or calibrated energies
  virtual double get_scint_gammas() const
  {
    PHOOL_VIRTUAL_WARN("get_scint_gammas()");
    return 0.0;
  }
  //! scint_gammas assigned to the tower. Depending on stage of process and DST node name, it could be scint_gammas deposition, light yield or calibrated energies
  virtual void set_scint_gammas(const double)
  {
    PHOOL_VIRTUAL_WARN("set_scint_gammas()");
    return;
  }

  //! energy assigned to the tower. Depending on stage of process and DST node name, it could be energy deposition, light yield or calibrated energies
  virtual double get_cerenkov_gammas() const
  {
    PHOOL_VIRTUAL_WARN("get_cerenkov_gammas()");
    return 0.0;
  }
  //! cerenkov_gammas assigned to the tower. Depending on stage of process and DST node name, it could be cerenkov_gammas deposition, light yield or calibrated energies
  virtual void set_cerenkov_gammas(const double)
  {
    PHOOL_VIRTUAL_WARN("set_cerenkov_gammas()");
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
  virtual CellConstRange get_g4cells() const;
  virtual CellIterator find_g4cell(CellKeyType id);
  virtual CellConstIterator find_g4cell(CellKeyType id) const;
  virtual void add_ecell(const CellKeyType /*g4cellid*/, const float /*ecell*/)
  {
    PHOOL_VIRTUAL_WARN("add_ecell(const CellKeyType g4cellid, const float ecell)");
    return;
  }
  virtual void clear_g4cells() {}

  virtual bool empty_g4showers() const { return true; }
  virtual size_t size_g4showers() const { return 0; }
  virtual ShowerConstRange get_g4showers() const;
  virtual ShowerIterator find_g4shower(int /*id*/);
  virtual ShowerConstIterator find_g4shower(int /*id*/) const;
  virtual void add_eshower(const int /*g4showerid*/, const float /*eshower*/)
  {
    PHOOL_VIRTUAL_WARN("add_eshower(const unsigned int g4showerid, const float eshower)");
    return;
  }
  virtual void clear_g4showers() {}

  //! Procedure to add a new PROPERTY tag:
  //! 1.add new tag below with unique value,
  //! 2.add a short name to RawTower::get_property_info
  enum PROPERTY
  {  //

    //! Scintillation photon count or energy
    prop_scint_gammas = 1,

    //! Cherenkov photon count or energy
    prop_cerenkov_gammas = 2,

    //! max limit in order to fit into 8 bit unsigned number
    prop_MAX_NUMBER = UCHAR_MAX
  };

  virtual bool has_property(const PROPERTY /*prop_id*/) const { return false; }
  virtual double get_property(const PROPERTY /*prop_id*/) const { return NAN; }
  virtual void set_property(const PROPERTY /*prop_id*/, const double /*value*/) { return; }
  static const std::string get_property_info(PROPERTY prop_id);

 protected:
  RawTower() {}

  virtual unsigned int get_property_nocheck(const PROPERTY /*prop_id*/) const { return UINT_MAX; }
  virtual void set_property_nocheck(const PROPERTY /*prop_id*/, const unsigned int) { return; }

  ClassDefOverride(RawTower, 1)
};

#endif
